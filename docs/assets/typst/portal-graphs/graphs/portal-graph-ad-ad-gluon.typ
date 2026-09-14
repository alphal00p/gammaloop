#import "../figure.typ": render

// UV-renormalization graph fixture: ad -> ad with a gluon loop.
#render(
  read("../../../../../tests/resources/graphs/uv_tests/ad_ad_1L_gluon.dot"),
  amplitude-mode: true,
)
