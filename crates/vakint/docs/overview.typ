// Overview section of the benchmark document (generated data only):
//   crates/vakint/docs/benchmark/tools/overview.py --data-root DATA_ROOT --tables tables.json [--overlay OVERLAY] --json overview.json
// `overview-section(o)`: what vakint computes today, symbolically or numerically and with what precision, from one to five loops:
// (a) counts per loop and Laurent order, (c) the result of any L-loop input, (d) the named constants, (e) how the forms were
// obtained and checked; PREVIEW mode adds the pending entries. (b), the per-master table, is precision table 1 (precision_tables.typ).

#let ov-sans = ("DejaVu Sans", "Roboto", "Noto Sans", "Libertinus Serif")
#let ov-navy = rgb("173f65")
#let ov-pale = rgb("eaf2fa")
#let ov-style = (
  "S": (fill: rgb("cdebc5"), ink: rgb("1b5e20")),
  "Z": (fill: rgb("e4efe1"), ink: rgb("4a6b45")),
  "S*": (fill: rgb("fde0b4"), ink: rgb("8a4b00")),
  "N": (fill: rgb("f8c9c4"), ink: rgb("8e1b12")),
  "P": (fill: rgb("e3dcf7"), ink: rgb("4b2a91")),
)
#let ov-na = rgb("f3f5f8")
#let ov-eps(o) = [$epsilon^(#o)$]
#let ov-key(o) = if o < 0 { "-" + str(-o) } else { str(o) }
#let ov-dag = text(fill: ov-style.P.ink, weight: "bold")[†]

#let ov-banner(o) = {
  let pre = o.mode == "preview"
  let tr = o.tree
  block(fill: if pre { ov-style.P.fill } else { ov-pale }, inset: 7pt, radius: 3pt, width: 100%, text(size: 8pt, font: ov-sans)[
    #text(weight: "bold", fill: if pre { ov-style.P.ink } else { ov-navy }, size: 9pt)[#if pre [PREVIEW] else [RELEASED]]
    #h(4pt) data tree #if tr.label != none [#tr.label] else [#raw(tr.data_root) (#tr.version)] · 5L basis #tr.five_loop_basis ·
    tables #raw(tr.tables_content_sha256.slice(0, 16)) (#tr.tables_status) · generated #o.generated_utc
    #if pre [
      #linebreak()
      Overlay #raw(o.overlay.path.split("/").last()) (#raw(o.overlay.sha256.slice(0, 12)), #o.overlay.created_utc): #o.overlay.title.
      Overlay results are *pending* (#ov-dag, violet): derived and checked numerically, but in no released data tree; they are
      counted apart and never as released (#o.overlay.states.pairs().map(((k, v)) => [#v #k]).join(", ")).
      #if o.overlay.basis != none [
        5L on #o.overlay.basis.name: #o.overlay.basis.added.map(x => [#x.name #x.design_name (replaces #x.replaces)]).join(", ");
        #o.overlay.basis.status.
      ]
    ] else [
      #linebreak()
      Only what this data tree ships. Results derived since the release (R operation, new basis masters, adopted PSLQ forms) are
      in the PREVIEW edition of this document.
    ]
  ])
}

#let ov-legend(pre) = {
  set text(7.4pt, font: ov-sans)
  let sw(k, body) = box(baseline: 2pt, rect(height: 9pt, inset: (x: 3pt, y: 1.5pt), radius: 1.5pt,
    fill: ov-style.at(k).fill, text(size: 6.8pt, weight: "bold", fill: ov-style.at(k).ink, body)))
  grid(columns: 3, column-gutter: 14pt, row-gutter: 4pt, align: left + horizon,
    [#sw("S")[S] exact, periods at any precision],
    [#sw("Z")[0] structural zero],
    [#sw("S*")[S\*(X)] exact, periods to X digits],
    [#sw("N")[N(X)] numerical only, X digits¹],
    ..if pre { ([#sw("P")[†] pending (preview)],) } else { () },
    [#box(baseline: 2pt, rect(height: 9pt, inset: (x: 3pt, y: 1.5pt), radius: 1.5pt, fill: ov-na, text(size: 6.8pt, fill: luma(130))[—])) not available],
  )
  v(-0.2em)
  text(size: 6.8pt, fill: luma(80))[¹ Every published digit figure of a numerical value (N(X), numerical payloads, finite literals)
    is its validated claim minus one digit (user ruling 2026-11-11), which absorbs discrepancies of 0.3–1 digit between a symbolic
    form and the stored claim. The lifted constants behind S\*(X) keep their own digits.]
}

// (a) counts per loop and order: rows = loop x class (classes absent at a loop are left out), columns = Laurent orders
#let ov-summary(o) = {
  set text(font: ov-sans, size: 7pt)
  set par(justify: false, leading: 0.4em)
  let pre = o.mode == "preview"
  let ords = o.orders1
  let classes = (("S", [S]), ("Z", [0]), ("S*", [S\*(X)]), ("N", [N(X)]))
  let hdr(b) = table.cell(fill: ov-navy, text(fill: white, weight: "bold", b))
  table(
    columns: (2.2cm, 1.05cm) + ords.map(_ => 1fr) + (1.05cm,),
    align: (x, y) => if x <= 1 { left + horizon } else { center + horizon },
    inset: (x: 1.8pt, y: 2.2pt),
    stroke: (x, y) => (bottom: 0.4pt + luma(228)),
    table.header(hdr[loops], hdr[class], ..ords.map(k => hdr(ov-eps(k))), hdr[total]),
    ..for s in o.summary {
      let rel = s.kind == "released"
      let get(k, c) = { let x = s.orders.at(ov-key(k), default: none); if x == none { none } else { x } }
      let present = classes.filter(((c, _)) => ords.any(k => { let x = get(k, c); x != none and (x.at(c) > 0 or (not rel and x.pending.at(c) > 0)) }))
      let n = present.len() + 1
      let out = ()
      out.push(table.cell(rowspan: n, fill: if rel { ov-pale } else { ov-style.P.fill },
        text(weight: "bold", fill: if rel { ov-navy } else { ov-style.P.ink })[#s.label \ #text(weight: "regular", size: 6.2pt)[#s.masters masters#if rel and s.loop == 5 and s.basis != none [ \ basis #s.basis]#if not rel [ \ incl. pending #ov-dag]]]))
      // held row
      out.push(text(size: 6.2pt, fill: luma(110))[held])
      for k in ords { let x = get(k, "S"); out.push(if x == none { text(fill: luma(200))[—] } else { text(size: 6.2pt, fill: luma(110))[#x.held] }) }
      out.push([])
      for (c, name) in present {
        let st = ov-style.at(c)
        out.push(table.cell(fill: st.fill.lighten(40%), text(fill: st.ink, weight: "semibold", size: 6.6pt, name)))
        let tot = 0
        let ptot = 0
        for k in ords {
          let x = get(k, c)
          if x == none { out.push([]); continue }
          let v = x.at(c)
          let p = if rel { 0 } else { x.pending.at(c) }
          tot += v
          ptot += p
          let extra = if c == "S*" and x.sstar_digits != none {
            let (lo, hi) = (x.sstar_digits.at(0), x.sstar_digits.at(1))
            [#lo#if hi != lo [–#hi] d]
          } else if c == "N" and x.n_digits_min != none [≥#x.n_digits_min, med #calc.round(x.n_digits_median, digits: 1)]
          out.push(if v == 0 and p == 0 { text(fill: luma(195))[·] } else {
            table.cell(fill: st.fill.lighten(55%), {
              text(weight: "bold", fill: st.ink)[#v]
              if p > 0 { text(fill: ov-style.P.ink, weight: "bold")[ +#p#ov-dag] }
              if extra != none { linebreak(); text(size: 5.2pt, fill: st.ink, extra) }
            })
          })
        }
        out.push([#text(weight: "bold", fill: st.ink)[#tot]#if ptot > 0 [#text(fill: ov-style.P.ink, weight: "bold")[ +#ptot]]])
      }
      out
    },
  )
}

// (c) result of any L-loop input; letters = the notes of precision table 2
#let ov-marks(c) = {
  let m = c.at("marks", default: ())
  if m.len() > 0 { text(size: 5.6pt, weight: "regular", fill: luma(90), super(m.join(","))) }
}
#let ov-rcell(c) = {
  let st = c.state
  if st == "below" { return table.cell(fill: white, []) }
  if st == "none" { return table.cell(fill: ov-na, text(size: 8pt, fill: luma(170))[—]) }
  if st == "partial" {
    return table.cell(fill: ov-na, text(size: 6.4pt, fill: luma(105))[not every input#ov-marks(c) \ #text(size: 5.6pt)[others #c.returned]])
  }
  if st == "numerical" {
    return table.cell(fill: ov-style.N.fill, text(size: 6.4pt, fill: ov-style.N.ink)[*numerical* \ #text(size: 5.4pt)[v0.2.x #c.v02x; via the basis: not measured]])
  }
  if st == "exact" {
    let s = if c.at("sstar_X", default: none) != none { ov-style.at("S*") } else { ov-style.S }
    return table.cell(fill: s.fill, stroke: 1pt + ov-style.P.ink, text(size: 8pt, weight: "bold", fill: s.ink)[exact#ov-dag])
  }
  if st == "fails-closed" {
    return table.cell(fill: ov-style.P.fill, text(size: 6.4pt, fill: ov-style.P.ink)[*fails closed* \ #text(size: 5.6pt)[#c.open_pending / #c.n_contributors entering coeff. numerical#if c.open_released != c.open_pending [ (#c.open_released w/o pending)]]])
  }
  let s = ov-style.at(c.class)
  let pend = c.at("pending", default: false)
  table.cell(fill: s.fill, stroke: if pend { 1pt + ov-style.P.ink } else { none }, {
    text(size: 8pt, weight: "bold", fill: s.ink)[#if c.class == "S" [exact] else if c.class == "S*" [exact S\*] else [N]#if pend { ov-dag }]
    if c.at("guaranteed", default: false) { linebreak(); text(size: 6pt, weight: "bold", fill: s.ink)[guaranteed²] }
    if c.class != "S" { linebreak(); text(size: 6.6pt, fill: s.ink)[#if c.class == "N" [≥ ]#c.X d] }
    if c.at("others", default: none) != none { linebreak(); text(size: 5.4pt, fill: s.ink)[others #c.others] }
    ov-marks(c)
  })
}

#let ov-result(o) = {
  set text(font: ov-sans)
  set par(justify: false, leading: 0.4em)
  let o2 = o.orders2
  let rows = o.result + o.former
  table(
    columns: (2.2cm,) + o2.map(_ => 1fr),
    align: center + horizon,
    inset: (x: 1.5pt, y: 4pt),
    stroke: 0.6pt + white,
    table.header(
      table.cell(fill: ov-navy, text(size: 7.4pt, fill: white, weight: "bold")[input]),
      ..o2.map(k => table.cell(fill: ov-navy, text(size: 7.4pt, fill: white, weight: "bold", ov-eps(k)))),
    ),
    ..for r in rows {
      let pre = r.at("kind", default: "former") != "released"
      (table.cell(fill: if pre { ov-style.P.fill } else { ov-pale },
        text(size: 7.4pt, weight: "bold", fill: if pre { ov-style.P.ink } else { ov-navy })[#r.label#if r.at("replaced_by", default: none) != none [#linebreak()#text(size: 5.8pt, weight: "regular")[replaced by #r.replaced_by]]]),
       ..o2.map(k => ov-rcell(r.cells.at(ov-key(k)))))
    },
  )
}

// (d) named constants
#let ov-constants(o) = {
  set text(font: ov-sans, size: 6.2pt)
  set par(justify: false, leading: 0.35em)
  show regex("\^\{[^}]*\}|\^\([^)]*\)|\^[-−]?[0-9a-zA-Zεδ]+"): it => super(it.text.slice(1).trim("{").trim("}").trim("(").trim(")"))
  show raw: set text(size: 6pt)
  let hdr(b) = table.cell(fill: ov-navy, text(fill: white, weight: "bold", b))
  table(
    columns: (3.1cm, 1fr, 1.05cm, 1.05cm, 2.5cm, 4.1cm),
    align: (x, y) => if x in (2, 3, 4) { center + horizon } else { left + horizon },
    inset: (x: 2.2pt, y: 1.7pt),
    stroke: (x, y) => (bottom: 0.4pt + luma(230)),
    table.header(hdr[constant], hdr[definition], hdr[digits], hdr[any precision], hdr[coefficients (released)], hdr[traditional form (verified, not yet adopted)]),
    ..for c in o.constants {
      let fin = c.kind == "finite"
      (
        [#if c.display != none [#c.display #text(fill: luma(120))[(#raw(c.code))]] else [#raw(c.code)]],
        [#c.definition],
        if fin { text(fill: ov-style.at("S*").ink, weight: "semibold")[#c.digits] } else { text(fill: ov-style.S.ink)[any] },
        if fin { text(fill: ov-style.at("S*").ink)[no] } else { text(fill: ov-style.S.ink)[yes] },
        [#c.used_total#if c.used_in.len() > 0 [ #text(fill: luma(110))[(#c.used_in.pairs().map(((k, v)) => [#k #v]).join(", "))]]#if c.at("pending_uses", default: 0) > 0 [ #text(fill: ov-style.P.ink)[+#c.pending_uses#ov-dag]]],
        [#if c.traditional != none [#c.traditional#if c.at("traditional_check", default: none) != none [ #text(fill: luma(120))[(#c.traditional_check)]]] else [#text(fill: luma(170))[—]]],
      )
    },
  )
}

#let ov-pending(o) = {
  set text(font: ov-sans, size: 6.1pt)
  set par(justify: false, leading: 0.35em)
  let hdr(b) = table.cell(fill: ov-style.P.ink, text(fill: white, weight: "bold", b))
  table(
    columns: (2cm, 1.4cm, 4.6cm, 1fr, 1.9cm),
    align: (x, y) => if x in (1, 4) { center + horizon } else { left + horizon },
    inset: (x: 2.2pt, y: 1.6pt),
    stroke: (x, y) => (bottom: 0.4pt + luma(230)),
    table.header(hdr[coefficient], hdr[pending class], hdr[review state], hdr[how obtained (evidence)], hdr[agreement / claimed d]),
    ..for p in o.pending {
      let low = p.below_claim
      (
        [#p.loop\L/#p.master #ov-eps(p.order)],
        text(fill: ov-style.P.ink, weight: "semibold")[#p.label],
        [#p.state#if p.at("evidence_class", default: none) != none [ #text(fill: ov-style.P.ink, weight: "bold")[#{("exact": "exact derivation", "fit": "fitted counterterm", "R*": "confirmed by R*").at(p.evidence_class)}]]],
        [#if p.method != none [#p.method] else if p.kind == "new-master" [value of the new basis master] else [#p.kind]#if p.at("fit", default: none) != none and p.method == none [ · #p.fit]#if p.flags.len() > 0 [ · #text(fill: rgb("8e1b12"))[#p.flags.join("; ")]]],
        [#if p.agreement_digits != none [#text(fill: if low { rgb("8e1b12") } else { ov-style.S.ink }, weight: "semibold")[#p.agreement_digits]] else [—] / #if p.claimed_digits != none [#p.claimed_digits] else [—]#if low [ #text(fill: rgb("8e1b12"))[below claim]]],
      )
    },
  )
}

#let overview-section(o) = {
  [= What vakint computes today <overview>]
  v(-0.3em)
  ov-banner(o)
  v(0.15em)
  text(size: 9.4pt)[For every single-scale vacuum integral from one to five loops: which Laurent coefficients vakint returns
    *exactly* (S: in constants it evaluates at any precision; S\*(X): in constants known to X digits), which only *numerically*
    (N(X): X validated digits), and which are *structural zeros*. (a) counts the master coefficients, (c) states what a user gets
    for an arbitrary $L$-loop input, (d) lists the constants, (e) says how the forms were obtained and checked. The per-master
    table (b) is precision table 1 (section "Precision of the masters and of the results").]
  v(0.1em)
  ov-legend(o.mode == "preview")
  v(0.3em)
  text(9pt, weight: "bold", fill: ov-navy)[(a) Master coefficients per loop and Laurent order]
  text(7.2pt, fill: luma(80))[ — number of master coefficients in each class (held: masters holding the order); S\*: range of the digits of the limiting constants;
    N: fewest and median validated digits (rounded down). #if o.mode == "preview" [Violet rows: every coefficient counted once, in its
    best class; +n#ov-dag = how many of them come from the overlay (pending), apart from the released count before it.]]
  ov-summary(o)
  v(0.4em)
  text(9pt, weight: "bold", fill: ov-navy)[(c) Result of any single-scale vacuum input with $L$ loops]
  text(7.2pt, fill: luma(80))[ — worst case over the admitted inputs (RustRed backend), spurious poles of the reduction and
    lower-loop factors of 5L products included (letters: the notes of precision table 2). "Not every input": some inputs need master orders
    the data do not hold and get MissingMasterOrders (no value, never a wrong one). #if o.former.len() > 0 [Former masters: reduced
    onto the new basis; an order that was exact comes back exact once every numerical coefficient entering it is exact, otherwise
    vakint fails closed at that order (no value).]]
  ov-result(o)
  for n in o.at("r_operation_notes", default: ()) {
    text(size: 6.1pt, fill: luma(60), font: ov-sans)[² *#n.cell guaranteed exact* (R-operation): the pole coefficients of every
      single-scale vacuum integral are free of the 4L-only constants; verified exactly on vakint's own reduction data:
      #n.evidence. So the #n.constants.join(", ") terms that
      spurious poles bring in through #n.cancelled_orders.join(", ") cancel identically, and vakint already returns this order
      free of them. Proof: #n.proof. #linebreak()]
  }
  if o.result.any(r => r.kind == "preview") {
    let r = o.result.filter(r => r.kind == "preview").first()
    text(size: 6.6pt, fill: luma(80), font: ov-sans)[#r.label: #r.note.]
  }
  v(0.4em)
  text(9pt, weight: "bold", fill: ov-navy)[(e) How the symbolic forms were obtained and checked]
  block(width: 100%, text(size: 7pt)[
    #set par(justify: true, leading: 0.45em, spacing: 0.5em)
    #for p in o.method [*#p.title.* #p.text #parbreak()]
  ])
  pagebreak(weak: true)
  text(9pt, weight: "bold", fill: ov-navy)[(d) Named constants]
  text(7.2pt, fill: luma(80))[ — every constant in the exact coefficients of the data tree. Built-ins are evaluated by vakint at any
    requested precision; finite records (4L/constants.json) are known to the digits shown, which bound every S\* coefficient that
    contains them. #if o.names != none [Traditional forms (periods/l6 dictionary): #o.names.status; full formulas and the
    explicit integral definitions in the appendix below and in #raw("docs/benchmark/data/CONSTANTS-TRADITIONAL.md"). D5bar3 is
    undecided.]]
  ov-constants(o)
  if o.names != none {
    pagebreak(weak: true)
    text(9pt, weight: "bold", fill: ov-navy)[Appendix: constants, traditional forms and definitions]
    text(7.2pt, fill: luma(80))[ — #o.names.status. Formulas with more than 20 terms (D5bar2, D6bar2) are in
      #raw("docs/benchmark/data/CONSTANTS-TRADITIONAL.md"). Notation: w = exp(iπ/3); Li\_{a1..ar}(x1..xr) = Σ\_{k1>..>kr≥1} Π\_i x\_i^k\_i / k\_i^a\_i;
      Ls\_j^(k)(t) = −∫\_0^t φ^k ln^(j−1−k) |2 sin(φ/2)| dφ. #if o.names.convention != none [#o.names.convention]]
    {
      set text(font: ov-sans, size: 5.6pt)
      set par(justify: false, leading: 0.3em)
      show regex("\^\{[^}]*\}|\^\([^)]*\)|\^[-−]?[0-9a-zA-Zεδ]+"): it => super(it.text.slice(1).trim("{").trim("}").trim("(").trim(")"))
      table(columns: (3.2cm, 1fr), inset: (x: 2pt, y: 1.2pt), stroke: (x, y) => (bottom: 0.4pt + luma(232)),
        ..for c in o.constants.filter(c => c.at("name_kind", default: none) != none) {
          ([#raw(c.code)#if c.name_kind == "undecided" [ #text(fill: rgb("8e1b12"), weight: "bold")[undecided]]],
           [#if c.name_kind == "reduces" {
              if c.formula != none [#c.code = #c.formula] else [#c.formula_terms terms: see CONSTANTS-TRADITIONAL.md]
              if c.traditional_check != none [ #text(fill: luma(110))[(#c.traditional_check)]]
              if c.at("definition_full", default: none) != none [ #linebreak()#text(fill: luma(90))[Definition: #c.definition_full]]
            } else [#c.definition_full]])
        })
    }
  }
  if o.pending.len() > 0 or o.at("listed", default: ()).len() > 0 {
    pagebreak(weak: true)
    text(9pt, weight: "bold", fill: ov-style.P.ink)[Pending entries of this preview (#o.pending.len())]
    text(7.2pt, fill: luma(80))[ — agreement = digits shared with the stored numerical value it would replace (claimed digits as
      measured, before the −1 publication offset); "below claim": more than one digit below that claim.]
    ov-pending(o)
    let big = o.pending.filter(p => p.below_claim)
    text(size: 6.6pt, fill: luma(70), font: ov-sans)[Agreement more than one digit below the claim (#big.len()):
      #if big.len() == 0 [none] else [#big.map(p => [#p.master ε#super(str(p.order)) #p.agreement_digits / #p.claimed_digits]).join("; ")].
      Discrepancies up to one digit are absorbed by the −1 publication offset (user ruling 2026-11-11).]
  }
  if o.at("listed", default: ()).len() > 0 {
    v(0.4em)
    text(9pt, weight: "bold", fill: ov-style.P.ink)[Listed, not applied (#o.listed.len())]
    text(7.2pt, fill: luma(80))[ — numerical 5L poles whose candidate is marginal (to be confirmed), whose overall counterterm the
      fit did not identify, or which have no candidate yet; they stay numerical in every count.]
    {
      set text(font: ov-sans, size: 6.1pt)
      let hdr(b) = table.cell(fill: luma(110), text(fill: white, weight: "bold", b))
      table(columns: (2.3cm, 3.6cm, 1fr, 1.9cm), inset: (x: 2.2pt, y: 1.6pt), stroke: (x, y) => (bottom: 0.4pt + luma(230)),
        table.header(hdr[coefficient], hdr[state], hdr[evidence / method], hdr[agreement / claimed d]),
        ..for p in o.listed {
          ([#p.loop\L/#p.master #ov-eps(p.order)], [#p.state], [#if p.evidence != none [#p.evidence] else [#p.method]],
           [#if p.agreement_digits != none [#p.agreement_digits] else [—] / #if p.claimed_digits != none [#p.claimed_digits] else [—]])
        })
    }
  }
}
