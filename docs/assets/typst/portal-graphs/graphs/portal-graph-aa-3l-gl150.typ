#import "../figure.typ": render

// Generated process graph: aa -> aa, three loops, GL150.
#render(
  read(
    "../../../../../examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL150.dot",
  ),
  amplitude-mode: true,
  momentum-arrows: true,
)
