// Default V1 GammaLoop grid template.
#import "fig-index.typ": tree

#let _alignment(value) = if type(value) != str {
  value
} else if value == "center+horizon" {
  center + horizon
} else if value == "left+horizon" {
  left + horizon
} else if value == "right+horizon" {
  right + horizon
} else if value == "center+top" {
  center + top
} else if value == "center+bottom" {
  center + bottom
} else {
  panic("unsupported grid alignment " + repr(value))
}

#let _render(config) = {
  let page-format = config.at("page_format", default: "portrait")
  let grid-cols = config.at("columns", default: 3)
  let grid-rows = config.at("rows", default: 5)
  let grid-row-h = config.at("row_height", default: 1fr)
  let grid-align = _alignment(config.at("align", default: "center+horizon"))
  let top-margin = config.at("top_margin", default: 12mm)
  let bottom-margin = config.at("bottom_margin", default: 30mm)
  let left-margin = config.at("left_margin", default: 10mm)
  let right-margin = config.at("right_margin", default: 10mm)
  let cells-per-page = grid-cols * grid-rows
  let show-page-numbers = config.at("page_numbers", default: true)
  let graph-label-position = config.at("graph_label_position", default: "top")
  let show-graph-labels = graph-label-position != "none"
  let graph-label-size = config.at("graph_label_size", default: 8pt)
  let image-height = if show-graph-labels { 80% } else { 100% }

  set page(
    paper: "a4",
    margin: (
      top: top-margin,
      bottom: bottom-margin,
      left: left-margin,
      right: right-margin,
    ),
    flipped: page-format == "landscape",
    numbering: if show-page-numbers { "1 / 1" } else { none },
    number-align: right + top,
  )

  let render-card(fig, index) = box(
    inset: 0pt,
    width: 100%,
    height: 100%,
  )[
    #if show-graph-labels and graph-label-position == "top" [
      #text(size: graph-label-size)[Graph \##index]
      #v(3pt)
    ]

    #image(fig.path, width: 100%, height: image-height, fit: "contain")

    #if show-graph-labels and graph-label-position == "bottom" [
      #v(3pt)
      #text(size: graph-label-size)[Graph \##index]
    ]
  ]
  let empty-cell = box(inset: 0pt, width: 100%, height: 100%)[]
  let heading-level(depth) = if depth + 2 > 6 { 6 } else { depth + 2 }
  let render-node(node, path, depth) = {
    if path.len() > 0 {
      heading(level: heading-level(depth))[#path.join(" / ")]
      v(6pt)
    }
    if node.figures.len() > 0 {
      let cards = node.figures.enumerate(start: 1).map(((index, fig)) => render-card(fig, index))
      for start in range(0, cards.len(), step: cells-per-page) {
        let end = calc.min(start + cells-per-page, cards.len())
        let slice = cards.slice(start, end)
        let padded = if slice.len() < cells-per-page {
          slice + ((cells-per-page - slice.len()) * (empty-cell,))
        } else {
          slice
        }
        grid(
          columns: grid-cols,
          rows: grid-rows * (grid-row-h,),
          align: grid-align,
          gutter: 0pt,
          ..padded,
        )
        if end < cards.len() {
          pagebreak()
        }
      }
    }
    for name in node.order {
      let child = node.folders.at(name)
      render-node(child, path + (name,), depth + 1)
    }
  }

  if tree.folders.processes.folders.at("amplitudes", default: none) != none {
    [= Amplitudes]
    for name in tree.folders.processes.folders.amplitudes.order {
      [#name]
      render-node(tree.folders.processes.folders.amplitudes.folders.at(name), (), 0)
    }
  }
  if tree.folders.processes.folders.at("cross_sections", default: none) != none {
    [= Cross-section]
    for name in tree.folders.processes.folders.cross_sections.order {
      [#name]
      render-node(tree.folders.processes.folders.cross_sections.folders.at(name), (), 0)
    }
  }
}

#let render(config) = _render(config)
