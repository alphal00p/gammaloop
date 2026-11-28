#import "macros.typ": *

#set document(title: "One-loop renormalisation of the Standard Model in the unbroken phase with a non-anticommuting gamma5", author: "gamma5 working notes")
#set page(paper: "a4", margin: (x: 2.2cm, y: 2.4cm), numbering: "1")
#set text(font: "New Computer Modern", size: 10.5pt, lang: "en")
#show math.equation: set text(font: "New Computer Modern Math")
#set par(justify: true, leading: 0.62em)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)", supplement: [Eq.])
#show heading.where(level: 1): it => { pagebreak(weak: true); v(0.4em); it; v(0.6em) }
#show heading.where(level: 2): it => { v(0.5em); it; v(0.3em) }
#show link: set text(fill: rgb("#1f4e79"))
#show ref: set text(fill: rgb("#1f4e79"))
#show figure.caption: set text(size: 9.5pt)
#set table(inset: 5pt)

#align(center)[
  #text(size: 17pt, weight: "bold")[One-loop renormalisation of the Standard Model \ in the unbroken phase with a non-anticommuting $g5$]
  #v(0.5em)
  #text(size: 11.5pt)[A diagrammatic BMHV calculation with Slavnov--Taylor-restoring counterterms]
  #v(0.8em)
  #text(size: 10pt)[Working report, #datetime(year: 2026, month: 11, day: 27).display("[day] [month repr:long] [year]")]
]
#v(2.5em)
#align(center)[#block(width: 86%)[
  #set par(justify: true)
  #set align(left)
  #text(weight: "bold")[Abstract.] We carry out the complete one-loop renormalisation of the Standard Model with one fermion family in the unbroken phase, in dimensional regularisation with the Breitenlohner--Maison--'t Hooft--Veltman treatment of $g5$, in the Feynman gauge and the $overline("MS")$ scheme. Every chiral fermion receives a sterile partner, so that the only breaking of gauge invariance by the regularisation is the evanescent part of the fermion kinetic term. The calculation is purely diagrammatic. The pole parts of eighteen classes of one-particle-irreducible Green's functions give the divergent counterterms, including their evanescent parts. The one-loop graphs with one insertion of the breaking vertex give the finite counterterms that restore the Slavnov--Taylor identity. All counterterms, their Feynman rules, and the renormalised action in the unbroken and in the broken phase are given explicitly. Each step is stated as an operation on Feynman graphs and as the one-loop case of a construction defined at any loop order, built on the $Rop$-operation. Every class was computed twice independently, the results agree with the literature where it exists, and the restored identity was checked against loop functions with complete finite parts. The last section lists what the present tools lack to automate the recipe.
]]
#v(1em)

#outline(depth: 2, indent: auto)

#include "sections/01_introduction.typ"
#include "sections/02_method.typ"
#include "sections/03_setup.typ"
#include "sections/03b_operator_bases.typ"
#include "sections/04_classes.typ"
#include "sections/05_conclusion.typ"
#include "sections/06_automation.typ"
#include "sections/07_supporting_material.typ"

#bibliography("refs.bib", style: "american-physics-society", title: "References")
#context [#metadata(query(heading).map(h => (h.location().page(), h.level, if h.has("label") { str(h.label) } else { "" }))) <toc-data>]
