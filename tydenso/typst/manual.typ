#import "@preview/tidy:0.4.3"
#import "lib.typ" as tydenso

#let manifest = toml("typst.toml")
#let package-version = manifest.package.version
#let repository = "https://github.com/alphal00p/gammaloop"
#let spenso-guide = "https://symbolica.io/docs/python_api/community/spenso/index.html"
#let idenso-guide = "https://symbolica.io/docs/python_api/community/idenso/index.html"
#let accent = rgb("#704b7c")
#let pale-accent = rgb("#f5eef8")
#let muted = rgb("#61616b")

#set document(title: "Tydenso Manual", author: "Tydenso contributors")
#set page(
  paper: "a4",
  margin: (x: 20mm, top: 19mm, bottom: 18mm),
  header: context {
    if counter(page).get().first() > 1 {
      set text(size: 8pt, fill: muted)
      grid(columns: (1fr, auto), [Tydenso], [Version #package-version])
      line(length: 100%, stroke: 0.35pt + rgb("#d7cfda"))
    }
  },
  footer: context {
    if counter(page).get().first() > 1 {
      align(center, text(size: 8pt, fill: muted, counter(page).display("1")))
    }
  },
)
#set text(font: "Libertinus Serif", size: 10.5pt)
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.1.1")
#show heading: set text(font: "Libertinus Serif", fill: accent)
#show heading.where(level: 4): set heading(numbering: none)
#show link: set text(fill: accent)
#show raw: set text(font: "DejaVu Sans Mono")

#let callout(title, body) = block(
  width: 100%,
  breakable: true,
  inset: 9pt,
  radius: 3pt,
  fill: pale-accent,
  stroke: (left: 2.2pt + accent),
)[
  #text(weight: "bold", fill: accent)[#title]
  #body
]

#let example-preamble = "#import tydenso: *\n"
#let docs = tidy.parse-module(
  read("lib.typ"),
  name: "tydenso",
  scope: (tydenso: tydenso),
  preamble: example-preamble,
)
#let reference-style = {
  let style = tidy.utilities.get-style-functions(tidy.styles.default)
  style.show-example = (..args) => block(
    tidy.styles.default.show-example(..args),
    breakable: false,
  )
  style
}
#let worked-example(code) = tidy.styles.default.show-example(
  raw(code.text, lang: "typ", block: true),
  scope: (tydenso: tydenso),
  preamble: example-preamble,
  mode: "markup",
  dir: ttb,
  scale-preview: 100%,
)
#show raw.where(lang: "worked"): worked-example

#align(center)[
  #v(25mm)
  #text(size: 35pt, weight: "bold", fill: accent)[Tydenso]
  #v(5mm)
  #text(size: 16pt, fill: muted)[Symbolic tensor algebra inside Typst]
  #v(13mm)
  #text(size: 11pt)[User manual · version #package-version]
  #v(23mm)
  #block(width: 76%, inset: 14pt, radius: 5pt, fill: pale-accent)[
    Build indexed tensors as readable Typst values, transform them with Idenso,
    and print them with Spenso's own notation.
  ]
  #v(34mm)
  #link(repository)[Repository] ·
  #link(repository + "/blob/main/tydenso/LICENSE")[MIT license] ·
  #link(spenso-guide)[Spenso Python API] ·
  #link(idenso-guide)[Idenso Python API]
]

#pagebreak()
#text(size: 22pt, weight: "bold", fill: accent)[Contents]
#v(5mm)
#outline(title: none, depth: 2, indent: auto)
#pagebreak()

= Start with a contraction

Tydenso is the tensor-focused companion to Tymbolica. It has its own package,
manual, notation, and WebAssembly engine. You do not need Tymbolica to construct,
simplify, inspect, or display a tensor expression.

The first example contracts a Minkowski metric with a vector. A representation
is an ordinary Typst dictionary; `slot` turns an index label into another
dictionary; `vector` creates a callable rank-one tensor name.

```worked
#let V = mink(4)
#let mu = slot(V, 1)
#let nu = slot(V, 2)
#let p = vector("p")

#let before = math($#metric(V, mu, nu) #p(nu)$)
#let after = simplify-metrics(before)

$
  #to-typst(before)
  quad arrow.r.long
  #to-typst(after)
$
```

Here `mink` turns the integer slots `1` and `2` into $mu$ and $nu$. The
constructor calls already produce readable math. Their hidden, versioned
metadata contains the exact Symbolica Atom, so `math` does not have to infer a
tensor from its appearance. It returns Atom bytes for algebra; `to-typst`
renders transformed results.

== Installation

Until Tydenso is published separately, install the `tydenso/typst` directory
from this repository as its own local package. On Linux, place or symlink it at:

#raw(
  "~/.local/share/typst/packages/local/tydenso/" + package-version,
  lang: "text",
  block: true,
)

Then import it independently:

#raw(
  "#import \"@local/tydenso:" + package-version + "\": *",
  lang: "typ",
  block: true,
)

From a checkout, `#import "path/to/gammaloop/tydenso/typst/lib.typ": *` works as
well. The package directory already contains its compressed engine and loader.

== How this corresponds to the Python API

The public concepts follow Spenso's Python surface, adapted to what is natural
in Typst. Python objects can overload calls; Typst dictionaries cannot, so
index construction is the explicit `slot(V, 1)`. Tensor names are functions,
which keeps the pleasant `p(mu)` spelling.

#table(
  columns: (1fr, 1fr),
  inset: 7pt,
  stroke: 0.4pt + rgb("#d7cfda"),
  table.header([*Python*], [*Typst*]),
  [```python
V = Representation.mink(4)
mu = V("mu")
nu = V("nu")
p = TensorName("p")
expr = TensorName.g(mu, nu) * p(nu)
```],
  [```typst
#let V = mink(4)
#let mu = slot(V, 1)
#let nu = slot(V, 2)
#let p = vector("p")
#let expr = math($#metric(V, mu, nu) #p(nu)$)
```],
)

Representations and slots stay transparent dictionaries. Symbols and tensor
calls are visible Typst math with exact Atom metadata; completed algebraic
expressions are Atom bytes.

= Build tensor expressions

== Representations and slots are data

Built-in constructors cover the representations initialized by Spenso and
Idenso: `mink`, `euc`, `lor`, `bis`, `spf`, `cof`, `coad`, and `cos`. Use
`representation` when a model introduces another representation.

Each built-in owns a conventional palette for integer slots:

#table(
  columns: (auto, 1fr, 1.25fr, auto),
  inset: 5pt,
  stroke: 0.4pt + rgb("#d7cfda"),
  table.header([*Constructor*], [*Class*], [*Palette*], [*Base row*]),
  [`mink`], [inline metric], [$mu, nu, rho, sigma$], [top],
  [`euc`], [self-dual], [$i, j, k, l$], [top],
  [`lor`], [dualizable], [$mu, nu, rho, sigma$], [top],
  [`bis`], [self-dual], [$a, b, c, d$], [bottom],
  [`spf`], [dualizable], [$alpha, beta, gamma, delta$], [top],
  [`cof`], [dualizable], [$i, j, k, l$], [top],
  [`coad`], [self-dual], [$a, b, c, d$], [top],
  [`cos`], [dualizable], [$I, J, K, L$], [top],
)

Every palette starts at 1. After its fourth entry it repeats with a numeric
subscript, so `slot(mink(4), 5)` displays as $mu_1$. A dualizable
representation uses the same palette for its dual; only the script row flips.

A custom representation can name its automatic indices. The palette repeats;
each pass adds a subscript. Manually written indices use the same notation but
keep their own symbolic identity.

```worked
#let M = representation(
  "M",
  4,
  namespace: "example",
  self-dual: true,
  indices: ($mu$, $nu$),
)
#let T = tensor("T", namespace: "example")
#let expression = T(
  slot(M, 1),       // mu
  slot(M, 2),       // nu
  slot(M, 3),       // mu_1
  slot(M, $rho_2$), // exactly the index written here
)

#let detailed = notation(with-dim: true)
$ #to-typst(expression, notation: detailed) $
```

Here the qualified form distinguishes the indices as members of $M_4$. A
custom representation without `indices` keeps numeric display. A string is a
Symbolica identifier and uses Symbolica's native quoted Typst output; Typst
math such as `$rho_2$` instead creates an opaque index symbol carrying exactly
that display notation. Neither form consults the automatic palette.

```worked
#let color = cof(3)
#let i = slot(color, 1)
#let j-lower = slot(color, 2, dual: true)

#table(
  columns: (auto, 1fr),
  inset: 5pt,
  [representation], [#color.name],
  [dimension], [#color.dimension],
  [first index], [#i.index],
  [second index is dual], [#j-lower.dual],
)
```

Because those values are dictionaries, document code can validate dimensions,
generate index families in a loop, or attach its own metadata before asking the
plugin to construct an Atom. `dual-representation` returns the paired
representation; `metric`, `identity-tensor`, and `flat-tensor` accept either
index labels or already constructed slots.

== Tensor names carry symmetry

`tensor` mirrors the useful part of Spenso's `TensorName`: a name plus optional
symmetric, antisymmetric, cyclic, or linear behavior. Symbolica applies those
attributes while it constructs the function, not as a later cosmetic step.

```worked
#let V = euc(3)
#let mu-slot = slot(V, $mu$)
#let nu-slot = slot(V, $nu$)
#let F = tensor("F", antisymmetric: true)

#let cancellation = math($#F(mu-slot, nu-slot) + #F(nu-slot, mu-slot)$)
$ F^(mu nu) + F^(nu mu) = #to-typst(cancellation) $
```

Only one of `symmetric`, `antisymmetric`, and `cycle-symmetric` may be true for
one tensor name. The `linear` flag is independent.

== Write tensor algebra as math

Tensor and vector constructors are Typst functions. Interpolate their calls in
math with `#F(...)`; interpolate scalar symbols with `#mass`. Parsely reads the
ordinary arithmetic, while each annotated value contributes its exact Atom.

```worked
#let V = mink("D")
#let mu = slot(V, 1)
#let nu = slot(V, 2)
#let p = vector("p")
#let mass = symbol("m", namespace: "model")

#let shell = math($#metric(V, mu, nu) #p(mu) #p(nu) - #mass^2$)
$ #to-typst(shell) $
```

The structural `add`, `mul`, `neg`, `sub`, `div`, and `pow` functions remain
useful for generated expressions. They accept Atom bytes, annotated content,
numbers, slots, and representation dictionaries. Both paths construct the
same Atom; neither reconstructs tensors from printed subscripts.

== Give vectors mathematical labels

Raw Typst math in a tensor or vector argument describes how one atomic label
is written. It can contain arithmetic, attachments, fractions, calls, accents,
and the usual structural math elements. This is useful for momentum labels
that should stay visually rich without becoming part of the tensor algebra.

```worked
#let V = mink(4)
#let p = vector("p")

#let routed = p($arrow(x + y)$, V)
#let algebraic = p(math($x + y$), V)

#table(
  columns: (auto, 1fr),
  inset: 5pt,
  [one display-defined label], [$ #to-typst(routed) $],
  [a genuine Symbolica sum], [$ #to-typst(algebraic) $],
)
```

The first call stores one opaque Symbolica symbol together with a portable,
structured Typst display tree. The second call stores the actual sum $x + y$,
so Symbolica may inspect or transform it. Tymbolica and Rubi preserve the
display attachment without needing to understand it; Tydenso restores it when
the payload returns.

== Build open and closed Spenso chains

The chain helpers construct Spenso's actual Atom heads; they are not drawing
commands. Compact endpoints are ordinary rank-one tensors carrying a
representation, while `gamma(mu)` supplies Spenso's `in` and `out`
placeholders for a chain factor.

```worked
#let M = mink(4)
#let B = bis(4)
#let mu = slot(M, 1)
#let nu = slot(M, 2)

#let p = vector("p")
#let q = vector("q")
#let u = vector("u")
#let v = vector("v")

#let scalar = dot(p(1, M), q(2, M))
#let open = chain(
  u(1, B),
  v(2, B),
  gamma(mu),
  gamma(p(1, M)),
  gamma(nu),
)
#let closed = trace(
  B,
  cyclic(gamma(mu), gamma(p(1, M)), gamma(nu)),
)

$ #to-typst(scalar) $
$ #to-typst(open) $
$ #to-typst(closed) $
```

Use `gamma(mu, a, b)` when the bispinor slots are explicit. `chain` keeps an
ordered open sequence; `cyclic` marks the factor list of a closed sequence;
and `trace` combines that cycle with its representation. The
#link(repository + "/blob/main/tydenso/typst/examples/spenso-notation.typ")[standalone
notation example] also inspects the constructed trees and checks their exact
Spenso shapes.

Explicit slots are also valid chain endpoints; they render as endpoint scripts
on the chain body, while compact vectors produce the named bra and ket ends.

A string and Typst math deliberately mean different things. `slot(B, "mu")`
uses the semantic Symbolica identifier `mu`; Symbolica prints its
multi-character name as the quoted Typst text `"mu"`. `slot(B,
$mu$)` creates a different, opaque index whose attached notation prints as the
Greek $mu$. They are not interchangeable contraction labels.

```worked
#let B = bis(4)
#let q = vector("q")
#let identifier = slot(B, "mu")
#let notation = slot(B, $mu$)

#table(
  columns: (auto, 1fr),
  [semantic identifier `"mu"`], [$ #to-typst(q(identifier)) $],
  [Typst notation `$mu$`], [$ #to-typst(q(notation)) $],
)
```

Use palette integers for a representation's standard indices, and use math
content when the written notation itself defines a manual index. Plain `$mu$`
works when the name is unbound. After `#let mu = ...`, Typst resolves `$mu$` to
that binding; write `$std.sym.mu$` when a literal Greek symbol must remain
unambiguous.

= Transform tensors

Idenso provides the domain-specific transformations. The most common ones
simplify metrics, gamma matrices, and color structures; selective expanders and
index-wrapping functions are available for lower-level workflows.

== Contract a metric chain

The next expression contains two metrics and one vector. Simplifying metrics
eliminates both contracted dummy indices in one pass.

```worked
#let V = mink(4)
#let mu = slot(V, 1)
#let nu = slot(V, 2)
#let rho = slot(V, 3)
#let p = vector("p")

#let chain = mul(
  metric(V, mu, nu),
  metric(V, nu, rho),
  p(rho),
)
#let reduced = simplify-metrics(chain)

$
  #to-typst(chain)
  quad arrow.r.long
  #to-typst(reduced)
$
```

`list-dangling` returns the free indices as Atom payloads. Pass an element to
`to-typst`, `to-string`, or `inspect` just like any other expression.

== Know which family to call

#table(
  columns: (1.1fr, 2.2fr),
  inset: 6pt,
  stroke: 0.4pt + rgb("#d7cfda"),
  table.header([*Family*], [*Purpose*]),
  [`simplify-metrics`, `expand-metrics`, `expand-mink`, `expand-bis`],
  [Metric, Minkowski, and bispinor structure.],
  [`simplify-gamma`, `dirac-adjoint`],
  [Dirac chains and conjugation.],
  [`simplify-color`, `expand-color`],
  [Color generators and structure constants.],
  [`cook-function`, `cook-indices`, `wrap-dummies`, `wrap-indices`],
  [Canonical index organization and explicit wrappers.],
  [`to-dots`],
  [Rewrite supported contractions as dot products.],
)

The mathematical conventions and supported identities track
#link(idenso-guide)[Idenso's API].

= Display and inspect

== Control tensor notation in Typst

`to-typst` returns math content. Tydenso lays out slots, tensor arguments,
chains, dots, and traces with ordinary Typst functions, so notation can follow
the document's own visual language. `to-string` remains available when a
compact Spenso expression is useful for logs or debugging.

```worked
#let V = mink(4)
#let p = vector("p")
#let expression = p(1, slot(V, 1))

#let detailed = notation(with-dim: true, commas: true)

Rendered: $ #to-typst(expression, notation: detailed) $

Compact Spenso form:
#raw(to-string(expression), block: true)
```

Upper and lower indices align in matching columns. `bis` uses the conventional
bottom row; the other built-in base orientations use the top row. A
dualizable representation puts its dual orientation on the opposite row,
while a self-dual representation does not flip.

Compact vectors nested inside another tensor have three layouts. The default
`"ports"` layout keeps their bra or ket outside the tensor and marks the
contracted argument position with a hollow port. `"schoonschip"` puts the bold
vector label directly into that aligned column. `"call"` groups upper entries
before a semicolon and lower entries after it.

```worked
#let L = lor(4)
#let p = vector("p", namespace: "notation_layout")
#let A = tensor("A", namespace: "notation_layout")
#let expression = A(
  slot(L, 1),
  p(1, L),
  slot(L, 2),
  slot(L, 3, dual: true),
)

#let scripts = notation(tensor-layout: "schoonschip")
#let call = notation(tensor-layout: "call")

#grid(
  columns: 2,
  gutter: 0.8em,
  [aligned scripts], $ #to-typst(expression, notation: scripts) $,
  [row-grouped call], $ #to-typst(expression, notation: call) $,
)
```

The script layout preserves every original argument column. The call layout
preserves order within each row; the exact Atom annotation retains the complete
original order. A call containing only lower entries starts with a semicolon,
as in $A(; rho)$. Commas are on by default in this mode and can be disabled
with `commas: false`.

The most common layout switches live directly on `notation`: `tensor-layout`,
`with-dim`, `parens`, `commas`, `symbol-scripts`, `index-gap`, and `factor-gap`.
`print-settings` configures the separate compact string returned by
`to-string`.

Exact namespaced heads and complete calls can be restyled without changing the
Atom. A call renderer receives its arguments as both structured nodes and
already rendered math. Tydenso attaches the exact call after the renderer
returns, so the result can still be interpolated into `math` losslessly.

```worked
#let x = symbol("x", namespace: "notation_example")
#let f = function("f", namespace: "notation_example")
#let expression = f(add(x, 1))

#let display = notation(
  heads: ("notation_example::x": $xi$),
  calls: (
    "notation_example::f": ctx => {
      let (argument,) = ctx.visual-arguments
      $cal(F)[#argument]$
    },
  ),
)

$ #to-typst(expression, notation: display) $
```

Use `tags` for a named family and `classes` for Symbolica attributes such as
`"antisymmetric"`. Tags are canonical namespaced identities—for example
`"model::kinematic"`—so they remain unambiguous when expressions move between
plugins. Document overrides win over Tydenso's defaults and over display labels
carried by a payload.

== Inspect the Atom as CBOR data

Atom bytes remain the lossless transformation and interchange format. Before a
transformation, annotated constructor content carries those same bytes in
metadata. `inspect` accepts either form and decodes it to recursive Typst data.
Every node has a `kind`; function nodes also expose their full and short names,
arguments, and symmetry flags.

```worked
#let V = euc(3)
#let A = tensor("A", symmetric: true)
#let expression = A(slot(V, 1), slot(V, 2))
#let tree = inspect(expression)

#table(
  columns: (auto, 1fr),
  inset: 5pt,
  [node kind], [#tree.kind],
  [function], [#tree.short-name],
  [arguments], [#tree.arguments.len()],
  [symmetric], [#tree.symmetric],
)
```

#callout(
  [Why keep both forms?],
  [
    CBOR is convenient for layouts, diagnostics, and package-level APIs.
    Reconstructing algebra from that tree would lose Symbolica state and exact
    coefficient details, so transformations continue to consume Atom bytes.
  ],
)

= Work alongside Tymbolica

Tydenso and Tymbolica share the same Atom payload. A package can construct and
simplify tensors with Tydenso, then pass the result to Tymbolica for a general
algebraic operation. Custom representations keep the information Tydenso needs
to use them again, including their index palette.

```typst
#import "@local/tydenso:0.1.0" as tensors
#import "@local/tymbolica:0.1.0" as algebra

#let V = tensors.representation(
  "M",
  4,
  namespace: "model",
  self-dual: true,
  indices: ($mu$, $nu$),
)
#let p = tensors.vector("p")
#let expression = p(tensors.slot(V, 1))

#let expanded = algebra.expand(expression)
#tensors.to-typst(expanded)
```

Keep package versions aligned. Matrix payloads from Tymbolica are a different
format and are not Tydenso inputs.

= API reference

The groups below cover the imported top-level API. `init` returns the same
surface as a dictionary bound to a selected plugin module.

#let reference-groups = (
  (
    title: [Tensor construction],
    names: (
      "tensor", "vector", "symbol", "function", "representation", "mink", "euc", "lor", "bis",
      "spf", "cof", "coad", "cos", "slot", "metric", "identity-tensor",
      "flat-tensor", "dual-representation",
    ),
  ),
  (
    title: [Expression construction],
    names: ("math", "atom", "add", "mul", "neg", "sub", "div", "pow"),
  ),
  (
    title: [Spenso products and chains],
    names: ("dot", "gamma", "gamma0", "gamma5", "projp", "projm", "chain", "cyclic", "trace"),
  ),
  (
    title: [Printing and inspection],
    names: ("notation", "merge-notation", "to-typst", "print-settings", "to-string", "inspect"),
  ),
  (
    title: [Metric and index transformations],
    names: (
      "cook-function", "cook-indices", "expand-bis", "expand-metrics",
      "expand-mink", "expand-mink-bis", "list-dangling", "simplify-metrics",
      "to-dots", "wrap-dummies", "wrap-indices",
    ),
  ),
  (
    title: [Dirac and color transformations],
    names: (
      "dirac-adjoint", "expand-color", "simplify-color", "simplify-gamma",
    ),
  ),
  (
    title: [Advanced configuration],
    names: ("init",),
  ),
)

#let public-names = docs.functions.filter(
  doc => doc.name.slice(0, 1) != "_",
).map(doc => doc.name)
#let grouped-names = reference-groups.map(group => group.names).flatten()
#assert.eq(
  grouped-names.sorted(),
  public-names.sorted(),
  message: "Every public API function must appear in exactly one reference group.",
)

#for group in reference-groups [
  == #group.title

  #{
    let subset = docs
    subset.functions = group.names.map(name =>
      docs.functions.find(doc => doc.name == name)
    )
    subset.variables = ()
    tidy.show-module(
      subset,
      style: reference-style,
      show-module-name: false,
      show-outline: false,
      sort-functions: none,
      omit-private-definitions: true,
      first-heading-level: 2,
      break-param-descriptions: true,
    )
  }
]

= Compatibility and licensing

This manual describes Tydenso #package-version and Typst 0.15 or newer.
Tydenso's original source is released under the
#link(repository + "/blob/main/tydenso/LICENSE")[MIT License]. The bundled engine also
contains Spenso, Idenso, and Symbolica; their own upstream terms continue to
apply. In particular, the MIT license for this interface does not relicense
Symbolica. Read #link("https://symbolica.io/license/")[Symbolica's license]
before redistributing or deploying the WebAssembly bundle.

= Acknowledgements

Tydenso would not exist without #link("https://symbolica.io/")[Symbolica], the
exact algebra engine beneath its Atom representation. Thank you to Symbolica's
contributors, and to the GammaLoop contributors for
#link(spenso-guide)[Spenso] and #link(idenso-guide)[Idenso], whose tensor model,
printers, and identities define this package's mathematics.

Thanks also to #link("https://typst.app/universe/package/tidy/")[Tidy], which
powers the worked examples and generated reference, and to
#link("https://typst.app/universe/package/parsely/")[Parsely], which parses
Typst math while preserving the semantic Atom metadata.
