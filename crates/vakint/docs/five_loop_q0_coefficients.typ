#set document(title: "Saved five-loop zero-momentum coefficients", author: "ValentinHirschi")
#set page(paper: "a4", margin: (x: 17mm, top: 16mm, bottom: 16mm),
  footer: context align(right, text(size: 8pt, fill: rgb("687584"))[
    Vakint · Saved zero-momentum result · #counter(page).display("1")
  ]))
#set text(font: "Libertinus Serif", size: 10.3pt, fill: rgb("243343"))
#set par(justify: false, leading: 0.57em, spacing: 0.65em)
#set heading(numbering: none)
#show heading.where(level: 1): set text(size: 23pt, weight: "bold", fill: rgb("173f65"))
#show heading.where(level: 2): set text(size: 14.5pt, weight: "bold", fill: rgb("173f65"))
#show raw: set text(font: "DejaVu Sans Mono", size: 7.7pt)
#show link: set text(fill: rgb("225c92"))
#let ink = rgb("173f65")
#let blue = rgb("eaf2fa")
#let gray = rgb("f1f4f7")
#let green = rgb("edf6f0")
#let panel(body, fill: blue) = block(width: 100%, inset: 9pt, radius: 4pt, fill: fill, body)
#let small(body) = text(size: 8.3pt, fill: rgb("5b6876"), body)
#let metric(value, label) = [
  #text(size: 23pt, weight: "bold", fill: ink)[#value]
  #linebreak()#text(size: 9pt)[#label]
]

= Saved five-loop zero-momentum coefficients

#text(size: 11.5pt)[Inspect the completed reduction and its scalar coefficients without repeating IBPs.]

#v(4pt)
#table(columns: (1fr, 1fr, 1fr), inset: (x: 11pt, y: 12pt),
  stroke: 0.5pt + white, fill: blue,
  metric("55", "individual master identities"),
  metric("58", "formal master products"),
  metric("449", "real Laurent coefficients"),
)

The saved expression is the complete zero-momentum contraction of the four-rung
quark-loop fixture, with its native factors and unaveraged metric/color
contraction. Its 55 identities comprise 46 five-loop masters and nine lower-loop
masters: one at 1L, one at 2L, two at 3L and five at 4L. All remain formal,
including repeated factors and products. No five-loop master value is computed.

== What is expanded

Write the unit-mass result as $g^10 sum_(alpha) P_(alpha)(epsilon) C_(alpha)(epsilon)$,
where $g$ is the unchanged #raw("UFO::GC_11") and $P_(alpha)$ is a canonical
master product with its full multiplicities. The supplied scalar jets obey

$ C_(alpha)(epsilon) = sum_(k=v_(alpha))^5 c_(alpha,k) epsilon^k + cal(O)(epsilon^6). $

The original output already contains $d=4-2 epsilon$ and the normalization
$C_(alpha)(epsilon) = R_(alpha)(epsilon) exp(-5 gamma_("E") epsilon)$. Neither is applied again.
Only the derived coefficient artifacts set $M_("UV")=mu_("r")=1$; the scale symbol
is absent in this saved input and no scale factor is inserted. The original
retained artifact keeps its mass dependence and normalization unchanged.

#panel(fill: green)[
  *The remainder is coefficientwise.* $cal(O)(epsilon^6)$ describes each scalar
  $C_(alpha)$, with the masters left unexpanded. The masters themselves depend on
  $epsilon$ and can have poles, so this is not a fully evaluated physical
  amplitude through $epsilon^5$.
]

== Load or reproduce

Use the existing coherent Symbolica Community installation containing Vakint:

```sh
python crates/vakint/examples/five_loop_q0_coefficients.py verify
python crates/vakint/examples/five_loop_q0_coefficients.py regenerate \
  --input PATH/result.symbolica --output NEW
```

Verification reads the saved artifacts and their precision/provenance records.
Regeneration starts from the retained expression, groups full master products,
and evaluates only their scalar jets. Both workflows avoid graph preparation and
IBP reduction. Regeneration requires a new output directory.

#small[Artifacts: #raw("examples/data/five_loop_q0/") relative to the Vakint crate.
The companion uses ordinary public APIs in one Symbolica kernel, with no private
shared-library path, implicit download, or Mathematica requirement.]

#pagebreak()
== Precision and faithful interchange

Each combined scalar Laurent coefficient receives one tracked evaluation with
300 requested decimal digits and the fixed 32-digit guard allowance. There is
no automatic retry. The minimum measured nonzero precision in this artifact is
330 decimal digits; final rounding is checked separately.

#table(columns: (1.25fr, 1fr, 1.7fr), inset: 7pt,
  stroke: 0.4pt + rgb("d5dee7"),
  fill: (x, y) => if y == 0 { blue } else if calc.odd(y) { gray } else { white },
  table.header([*Component*], [*Stored precision*], [*Interpretation*]),
  [Nonzero real coefficient], [998 binary bits], [At least 300 decimal digits after export rounding],
  [Exact imaginary zero], [1104 binary bits], [Exact zero; storage width is not a finite accuracy limit],
)

The sidecar retains componentwise achieved precision, exact rational centers and
absolute error radii, error quality and warnings. Radii include final binary
rounding. An independent development oracle at 380 digits checked all 449 real
coefficients against those saved radii; it is not an end-user dependency.

#table(columns: (1.3fr, 2.05fr), inset: 7pt,
  stroke: 0.4pt + rgb("d5dee7"),
  table.header([*Artifact*], [*Role*]),
  [#raw("q0-retained.atom / .expr")], [Original complete symbolic reduction],
  [#raw("q0-coefficients-exact-…")], [Exact scalar jets through $epsilon^5$, masters formal],
  [#raw("q0-coefficients-300digits-…atom")], [Checked native Float values with component storage precision],
  [#raw("coefficient-table.json")], [Lossless exact-ratio/explicit-bit reconstruction plus uncertainty metadata],
  [#raw("…300digits-…expr")], [Human decimal display; not a lossless Float round trip],
  [#raw("provenance.json")], [Checksums, source bindings, validation and preserved failure history],
)

The native Atom and exact-ratio JSON reconstruct identically. Canonical decimal
printing can lose a binary bit when parsed again; numerical #raw(".expr") files
therefore remain displays. The original exact text and exact symbolic jets have
their separate checks.

== Scientific and postprocessing provenance

#panel[
  *Original reduction:* executable #raw("490fd6492677"), external data
  #raw("d4384149b507"), retained expression #raw("da909dc0b8f1").
  *Coefficient postprocessor:* actual runtime source #raw("4111049372c0"),
  host executable #raw("67937cb5291d"). The literal production-source comparison
  to #raw("03e80d538c6e") covers 676 compiled production inputs; it does not
  relabel the executable as a rebuild of that later commit. Full identities are
  recorded in the provenance sidecar.
]

The 58 exact jets and tracked passes completed in *3.775 s*. A strict export
check then exposed a storage mismatch: native exact-zero imaginary components
kept 1104 bits, while JSON had recorded 998. A *2.205 s* metadata-only completion
repaired those 449 exact-zero precision records and performed the independent validation.
The failed receipt remains preserved. Every native and expression file stayed
byte-for-byte unchanged; no coefficient was reevaluated.

#small[These are development #emph[postprocessing] observations, not IBP reduction
benchmarks or cache speedups. The read-only repair codec is pinned specifically
to Symbolica #raw("9cd57678d8817ec33c05cacba4925162f1c86b0b"); it makes no claim
about future native formats. See
#link("benchmark/benchmark.typ")[the reduction benchmark] for reduction measurements.]
