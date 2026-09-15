#import "../lib.typ": *

#set page(width: auto, height: auto, margin: 12pt)

#let cases = (
  (name: "mink", rep: mink(4)),
  (name: "euc", rep: euc(4)),
  (name: "lor", rep: lor(4)),
  (name: "bis", rep: bis(4)),
  (name: "spf", rep: spf(4)),
  (name: "cof", rep: cof(3)),
  (name: "coad", rep: coad(8)),
  (name: "cos", rep: cos(6)),
)

#let palette-rows = cases.map(case => {
  let T = tensor("T_" + case.name, namespace: "index_palette_example")
  let expression = T(
    slot(case.rep, 1),
    slot(case.rep, 2),
    slot(case.rep, 5),
  )
  ([#case.name], $#to-typst(expression)$)
}).flatten()

#grid(
  columns: (auto, 1fr),
  column-gutter: 1em,
  row-gutter: 0.45em,
  ..palette-rows,
)

#let p = vector("p")
#let manual = p(slot(mink(4), $zeta_7$))

#v(0.8em)
#grid(
  columns: (auto, 1fr),
  column-gutter: 1em,
  row-gutter: 0.45em,
  [first automatic index], $#to-typst(p(slot(mink(4), 1)))$,
  [after one palette wrap], $#to-typst(p(slot(mink(4), 5)))$,
  [manual index], $#to-typst(manual)$,
)
