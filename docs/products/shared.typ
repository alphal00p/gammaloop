// Shared presentation and routing helpers for every project manual.

#let docs-channel = sys.inputs.at("channel", default: "latest")
#let docs-snapshot-tag = sys.inputs.at("snapshot-tag", default: "")
#let source-revision = sys.inputs.at("git-commit", default: "main")

#let ink = rgb("#3d2645")
#let muted = rgb("#625b63")
#let accent = rgb("#745080")
#let accent-soft = rgb("#eadff0")
#let rule = rgb("#d9cfd9")
#let paper = rgb("#ffffff")
#let canvas = rgb("#f9f6f0")

#let product-link(id, label: none, path: "") = {
  let shown = if label == none { id } else { label }
  let route = if docs-channel == "snapshot" {
    "../../../" + id + "/snapshots/" + docs-snapshot-tag + "/"
  } else {
    "../../" + id + "/latest/"
  }
  link(route + path)[#shown]
}

#let source-link(path, label: none) = {
  let shown = if label == none { path } else { label }
  link("https://github.com/alphal00p/gammaloop/blob/" + source-revision + "/" + path)[#shown]
}

#let developer-link(id, source, label) = {
  if docs-channel == "snapshot" {
    source-link("docs/architecture/" + source, label: label + " at this revision")
  } else {
    link("../../../developers/architecture/" + id + "/")[#label]
  }
}

#let callout(title, body) = block(
  width: 100%,
  fill: accent-soft,
  stroke: (left: 3pt + accent),
  inset: (x: 11pt, y: 8pt),
  radius: 2pt,
)[
  #text(fill: ink, weight: "bold")[#title]
  #linebreak()
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
  #linebreak()
  #body
]

#let catalog-contract(rust-scope: none, python-scope: none) = {
  let rust-line = if rust-scope == none {
    []
  } else {
    [Rust packages: #raw(rust-scope)]
  }
  let python-line = if python-scope == none {
    []
  } else {
    [Python modules: #raw(python-scope)]
  }

  callout("API reference", [
    The reference for this version lists public items, signatures, fields, defaults, feature
    requirements, and links to their source definitions.

    #rust-line  #linebreak() #python-line
  ])
}

#let release-note(body) = callout("History coverage", body)

#let product-document(
  title: none,
  tagline: none,
  version: none,
  owner: none,
  body: [],
) = {
  set document(title: title)
  set page(paper: "a4", margin: (x: 20mm, y: 18mm), fill: canvas)
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
    fill: paper,
    stroke: .7pt + ink,
    inset: 9pt,
    radius: 3pt,
  )[#it]

  align(center)[
    #text(size: 25pt, weight: "bold", fill: ink)[#title]
    #v(5pt)
    #text(size: 12pt, fill: muted)[#tagline]
  ]

  v(15pt)
  line(length: 100%, stroke: .8pt + rule)
  v(8pt)
  body
}
