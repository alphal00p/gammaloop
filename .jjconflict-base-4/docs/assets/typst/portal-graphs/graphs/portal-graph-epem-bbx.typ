#import "../figure.typ": render

// UV-renormalization graph fixture: e+ e- scattering with a d-quark loop.
#render(
  read("../../../../../tests/resources/graphs/uv_tests/epem_a_bbx.dot"),
  amplitude-mode: true,
)
