// Shared presentation and routing helpers for every product manual.
//
// Data ownership is deliberately split three ways:
// - registry.toml owns product identity, component versions, feature gates, and relationships;
// - generated Rust/Python catalogs own signatures and item-level API metadata;
// - the files below docs/products/<product>/content own narrative guidance and examples.
//
// Product links climb from products/<current>/<channel>/ to the shared products root.

#let docs-channel = sys.inputs.at("channel", default: "latest")
#let docs-snapshot-tag = sys.inputs.at("snapshot-tag", default: "")
#let source-revision = sys.inputs.at("git-commit", default: "main")

#let ink = rgb("#172033")
#let muted = rgb("#586174")
#let accent = rgb("#5b4bdb")
#let accent-soft = rgb("#efedff")
#let rule = rgb("#d9dce5")
#let paper = rgb("#fbfbfd")

#let product-link(id, label: none) = {
  let shown = if label == none { id } else { label }
  let route = if docs-channel == "snapshot" {
    "../../../" + id + "/snapshots/" + docs-snapshot-tag + "/"
  } else {
    "../../" + id + "/latest/"
  }
  link(route)[#shown]
}

#let source-link(path, label: none) = {
  let shown = if label == none { path } else { label }
  link("https://github.com/alphal00p/gammaloop/blob/" + source-revision + "/" + path)[#shown]
}

#let callout(title, body) = block(
  width: 100%,
  fill: accent-soft,
  stroke: (left: 3pt + accent),
  inset: (x: 11pt, y: 8pt),
  radius: 2pt,
)[
  #text(fill: ink, weight: "bold")[#title]
  #v(3pt)
  #body
]

#let boundary(title, body) = block(
  width: 100%,
  fill: paper,
  stroke: .7pt + rule,
  inset: 11pt,
  radius: 3pt,
)[
  #text(fill: ink, weight: "bold")[#title]
  #v(4pt)
  #body
]

#let catalog-contract(rust-scope: none, python-scope: none) = {
  let rust-line = if rust-scope == none {
    [No Rust catalog is mounted for this page.]
  } else {
    [Rust scope: #raw(rust-scope)]
  }
  let python-line = if python-scope == none {
    [No Python catalog is mounted for this page.]
  } else {
    [Python scope: #raw(python-scope)]
  }

  callout("Generated-reference contract", [
    Signatures, fields, variants, defaults, feature gates, and source locations come from
    ordered generated descriptors. Hand-written prose explains intent and workflow; it must
    not duplicate a signature as a second source of truth. Rust documentation explicitly
    marked as Typst markup may be evaluated as markup. Python docstrings and all other
    external metadata are escaped as plain text.

    #rust-line  #linebreak() #python-line
  ])
}

#let release-note(body) = callout("Release-history policy", body)

#let product-document(
  title: none,
  tagline: none,
  version: none,
  owner: none,
  body: [],
) = {
  set document(title: title)
  set page(paper: "a4", margin: (x: 20mm, y: 18mm))
  set text(size: 10.5pt, fill: ink)
  set par(justify: true, leading: .66em)
  set heading(numbering: "1.1")
  show heading: it => {
    v(10pt)
    set text(fill: ink)
    it
    v(2pt)
  }
  show link: set text(fill: accent)
  show raw.where(block: true): it => block(
    width: 100%,
    fill: rgb("#f1f3f7"),
    inset: 9pt,
    radius: 2pt,
  )[#it]

  align(center)[
    #text(size: 25pt, weight: "bold", fill: ink)[#title]
    #v(5pt)
    #text(size: 12pt, fill: muted)[#tagline]
    #v(8pt)
    #text(size: 8.5pt, fill: muted)[Version source: #raw(version) #h(1em) Owner: #owner]
  ]

  v(15pt)
  line(length: 100%, stroke: .8pt + rule)
  v(8pt)
  body
}
