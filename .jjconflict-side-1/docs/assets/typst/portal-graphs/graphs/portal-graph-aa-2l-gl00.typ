#import "../figure.typ": render

// Saved process graph: aa -> aa, two loops, GL00.
#render(
  read("../../../../../examples/cli/aa_aa/2L/graphs/GL00.dot"),
  amplitude-mode: true,
)
