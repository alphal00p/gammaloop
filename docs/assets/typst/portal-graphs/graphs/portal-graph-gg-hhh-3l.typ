#import "../figure.typ": render

// Saved process graph: gg -> hhh, three loops.
#render(
  read("../../../../../examples/cli/gg_hhh/3L/3L_graph.dot"),
  amplitude-mode: true,
)
