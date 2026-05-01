// Default grid template that mirrors the folder hierarchy of generated figures.
#import "fig-index.typ": tree

#set page(margin: 12mm)

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

#let render-node(node, path, depth) = {
  if path.len() > 0 {
    heading(level: heading-level(depth))[#path.join(" / ")]
    v(6pt)
  }

  if node.figures.len() > 0 {
    grid(
      columns: eval(sys.inputs.at("columns", default: "3")),
      gutter: 12pt,
      ..node.figures.map(render-card)
    )
    v(12pt)
  }

  for name in node.order {
    let child = node.folders.at(name)
    render-node(child, path + (name,), depth + 1)
  }
}

#if tree.figures.len() == 0 and tree.order.len() == 0 {
  text(fill: gray)[No figures were generated.]
} else {
  render-node(tree, (), 0)
}
