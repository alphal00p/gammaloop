#import "../figure.typ": render

// Saved process graph: aa -> aa, two loops, GL08.
#render(
  read("../../../../../examples/cli/aa_aa/2L/graphs/GL08.dot"),
  amplitude-mode: true,
  momentum-arrows: true,
)
