// Default GammaLoop grid template.
#import "fig-index.typ": tree

#let page_format = sys.inputs.at("page_format", default: "portrait")
#let grid_cols = int(sys.inputs.at("columns", default: "3"))
#let grid_rows = int(sys.inputs.at("rows", default: "5"))
#let grid_row_h = eval(sys.inputs.at("row_height", default: "1fr"))
#let grid_align = eval(sys.inputs.at("align", default: "center+horizon"))
#let top_margin = eval(sys.inputs.at("top_margin", default: "12mm"))
#let bottom_margin = eval(sys.inputs.at("bottom_margin", default: "30mm"))
#let left_margin = eval(sys.inputs.at("left_margin", default: "10mm"))
#let right_margin = eval(sys.inputs.at("right_margin", default: "10mm"))
#let cells_per_page = grid_cols * grid_rows

#let page_numbers = sys.inputs.at("page_numbers", default: "on")
#let show_page_numbers = if type(page_numbers) == str {
  let value = page_numbers
  not (
    value == "off"
      or value == "OFF"
      or value == "none"
      or value == "NONE"
      or value == "false"
      or value == "FALSE"
  )
} else {
  true
}

#let graph_label_position = sys.inputs.at("graph_label_position", default: "top")
#let show_graph_labels = graph_label_position != "none"
#let graph_label_size = eval(sys.inputs.at("graph_label_size", default: "8pt"))
#let image_height = if show_graph_labels { 80% } else { 100% }

#set page(
  paper: "a4",
  margin: (
    top: top_margin,
    bottom: bottom_margin,
    left: left_margin,
    right: right_margin,
  ),
  flipped: page_format == "landscape",
  numbering: if show_page_numbers { "1 / 1" } else { none },
  number-align: right + top,
)

#let render-card(fig, index) = box(
  inset: 0pt,
  width: 100%,
  height: 100%,
)[
  #if show_graph_labels and graph_label_position == "top" [
    #text(size: graph_label_size)[Graph \##index]
    #v(3pt)
  ]

  #image(
    fig.path,
    width: 100%,
    height: image_height,
    fit: "contain",
  )

  #if show_graph_labels and graph_label_position == "bottom" [
    #v(3pt)
    #text(size: graph_label_size)[Graph \##index]
  ]
]

#let empty-cell = box(
  inset: 0pt,
  width: 100%,
  height: 100%,
)[]

#let heading-level(depth) = if depth + 2 > 6 { 6 } else { depth + 2 }

#let render-node(node, path, depth) = {
  if path.len() > 0 {
    heading(level: heading-level(depth))[#path.join(" / ")]
    v(6pt)
  }

  if node.figures.len() > 0 {
    let cards = node.figures.enumerate(start: 1).map(((index, fig)) => render-card(fig, index))

    for start in range(0, cards.len(), step: cells_per_page) {
      let end = calc.min(start + cells_per_page, cards.len())
      let slice = cards.slice(start, end)
      let padded = if slice.len() < cells_per_page {
        slice + ((cells_per_page - slice.len()) * (empty-cell,))
      } else {
        slice
      }

      grid(
        columns: grid_cols,
        rows: grid_rows * (grid_row_h,),
        align: grid_align,
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

#if tree.folders.processes.folders.at("amplitudes", default: none) != none {
  [= Amplitudes]
  for name in tree.folders.processes.folders.amplitudes.order {
    [#name]
    render-node(tree.folders.processes.folders.amplitudes.folders.at(name), (), 0)
  }
}

#if tree.folders.processes.folders.at("cross_sections", default: none) != none {
  [= Cross-section]
  for name in tree.folders.processes.folders.cross_sections.order {
    [#name]
    render-node(tree.folders.processes.folders.cross_sections.folders.at(name), (), 0)
  }
}
