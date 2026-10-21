#import "../figure.typ": render

// Current e+ e- -> t tbar fixture. The portal adapter opens its sewn
// initial-state connections for this drawing; it does not apply a final-state cut.
#render(
  read("../../../../../tests/resources/graphs/epemttbar.dot"),
  cross-section-mode: true,
  momentum-arrows: true,
)
