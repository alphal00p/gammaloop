#import "../figure.typ": render

// Cross-section fixture: e+ e- -> t tbar with a physical cut.
#render(
  read("../../../../../tests/resources/graphs/epemttbar.dot"),
  cross-section-mode: true,
  momentum-arrows: true,
)
