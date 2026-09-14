#import "../figure.typ": render

// GammaLoop graph fixture: gg -> hhh.
#render(
  read("../../../../../tests/resources/graphs/gghhh.dot"),
  amplitude-mode: true,
)
