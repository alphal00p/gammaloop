// Default V1 grid template that mirrors the folder hierarchy of generated figures.
#import "fig-index.typ": tree

#let render-card(fig) = box(
  inset: 10pt,
  radius: 6pt,
  stroke: 0.6pt + gray,
  fill: luma(98%),
  width: 100%,
)[
  #image(fig.path, width: 100%)
  #v(6pt)
  #text(size: 9pt, fill: gray)[#fig.relative]
]

#let heading-level(depth) = if depth + 2 > 6 { 6 } else { depth + 2 }

#let render-node(node, path, depth, columns) = {
  if path.len() > 0 {
    heading(level: heading-level(depth))[#path.join(" / ")]
    v(6pt)
  }

  if node.figures.len() > 0 {
    grid(columns: columns, gutter: 12pt, ..node.figures.map(render-card))
    v(12pt)
  }

  for name in node.order {
    let child = node.folders.at(name)
    render-node(child, path + (name,), depth + 1, columns)
  }
}

#let render(config) = {
  set page(margin: 12mm)
  let columns = config.at("columns", default: 3)
  if tree.figures.len() == 0 and tree.order.len() == 0 {
    text(fill: gray)[No figures were generated.]
  } else {
    render-node(tree, (), 0, columns)
  }
}
