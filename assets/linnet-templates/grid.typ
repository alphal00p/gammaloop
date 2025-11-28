// Default grid template that mirrors the folder hierarchy of generated figures.
#import "fig-index.typ": cols, tree

#set page(margin: 12mm)

#let render-card(fig) = box(
  inset: 0pt,
  width: 100%,
)[
  #image(fig.path, width: 100%)
]

#let heading-level(depth) = if depth + 2 > 6 { 6 } else { depth + 2 }

#let render-node(node, path, depth) = {
  if path.len() > 0 {
    heading(level: heading-level(depth))[#path.join(" / ")]
    v(6pt)
  }

  if node.figures.len() > 0 {
    grid(
      columns: 5,
      gutter: 0pt,
      ..node.figures.map(render-card)
    )
  }

  for name in node.order {
    let child = node.folders.at(name)
    render-node(child, path + (name,), depth + 1)
  }
}


#if tree.folders.processes.folders.amplitudes != none {
  [= Amplitudes]
  for f in tree.folders.processes.folders.amplitudes.order {
    [#f]
    render-node(tree.folders.processes.folders.amplitudes.folders.at(f), (), 0)
  }
}

#if tree.folders.processes.folders.at("cross-section", default: none) != none {
  [= Cross-Section]
  for f in tree.folders.processes.folders.cross-section.order {
    [#f]
    render-node(
      tree.folders.processes.folders.cross-section.folders.at(f),
      (),
      0,
    )
  }
}
