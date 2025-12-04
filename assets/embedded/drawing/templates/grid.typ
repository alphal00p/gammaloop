// Default grid template that mirrors the folder hierarchy of generated figures.
#import "fig-index.typ": tree

// --- Global layout configuration from inputs ---------------------------------

#let page_format    = sys.inputs.at("page_format", default: "portrait")
#let grid_cols      = int(sys.inputs.at("columns", default: 3))
#let grid_rows      = int(sys.inputs.at("rows", default: 5))
#let grid_row_h     = eval(sys.inputs.at("row_height", default: "1fr"))
#let grid_align     = eval(sys.inputs.at("align", default: "center+horizon"))
#let cells_per_page = grid_cols * grid_rows

// Page-number option: "on" (default) or "off"/"none"/"false"
#let page_numbers      = sys.inputs.at("page_numbers", default: "on")
#let show_page_numbers = if type(page_numbers) == "str" {
    let v = page_numbers.lower()
    not (v == "off" or v == "none" or v == "false")
  } else {
    true
  }

// Graph label position: "top", "bottom", or "none" (default: "bottom")
#let graph_label_position = sys.inputs.at("graph_label_position", default: "top")

#let show_graph_labels = graph_label_position != "none"

// Graph label font size (e.g. "8pt", "10pt")
#let graph_label_size = eval(sys.inputs.at("graph_label_size", default: "8pt"))

// Image height inside cell (reserve space if label is shown)
#let image_height = if show_graph_labels { 80% } else { 100% }

#set page(
  paper: "a4",
  margin: 12mm,
  flipped: page_format == "landscape",
  // "1 / 1" pattern = current page / total pages
  numbering: if show_page_numbers { "1 / 1" } else { none },
  // top-right
  number-align: right + top,
)

// --- Card / heading helpers --------------------------------------------------

// Each card fills the grid cell. When labels are enabled, the image uses
// only part of the cell height so there is room for the label inside the cell.
#let render-card(fig, idx) = box(
  inset: 0pt,
  width: 100%,
  height: 100%,
)[
  // Label on top, if requested
  #if show_graph_labels and graph_label_position == "top" [
    #text(size: graph_label_size)[Graph \##idx]
    #v(3pt)
  ]

  // Image
  #image(
    fig.path,
    width: 100%,
    height: image_height,
    fit: "contain",
  )

  // Label on bottom, if requested
  #if show_graph_labels and graph_label_position == "bottom" [
    #v(3pt)
    #text(size: graph_label_size)[Graph \##idx]
  ]
]

// empty cell used for padding so layout is always Ncol × Nrow
#let empty-cell = box(
  inset: 0pt,
  width: 100%,
  height: 100%,
)[
]

#let heading-level(depth) = if depth + 2 > 6 { 6 } else { depth + 2 }

// --- Recursive rendering with paging ----------------------------------------

#let render-node(node, path, depth) = {
  if path.len() > 0 {
    heading(level: heading-level(depth))[#path.join(" / ")]
    v(6pt)
  }

  if node.figures.len() > 0 {
    // Build cards with indices starting from 1
    let cards = node.figures.enumerate(start: 1).map(((i, fig)) => render-card(fig, i))

    // Chunk into pages of at most grid_cols * grid_rows items
    for i in range(0, cards.len(), step: cells_per_page) {
      let end = if i + cells_per_page > cards.len() {
        cards.len()
      } else {
        i + cells_per_page
      }

      let slice = cards.slice(i, end)

      // Pad with empty cells so we always have exactly cells_per_page items
      let padded = if slice.len() < cells_per_page {
        let pad = cells_per_page - slice.len()
        slice + (pad * (empty-cell,))
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

// --- Top-level sections ------------------------------------------------------

#if tree.folders.processes.folders.at("amplitudes", default: none) != none {
  [= Amplitudes]
  for f in tree.folders.processes.folders.amplitudes.order {
    [#f]
    render-node(tree.folders.processes.folders.amplitudes.folders.at(f), (), 0)
  }
}

#if tree.folders.processes.folders.at("cross_sections", default: none) != none {
  [= Cross-Section]
  for f in tree.folders.processes.folders.cross_sections.order {
    [#f]
    render-node(
      tree.folders.processes.folders.cross_sections.folders.at(f),
      (),
      0,
    )
  }
}
