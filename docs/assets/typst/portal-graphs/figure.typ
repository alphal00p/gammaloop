#import "layout.typ": layout, layout-with-cut-curves
#import "../theme.typ": palette
#import "../../../../assets/embedded/drawing/templates/physics-edge-style.typ": (
  momentum-arrow-defaults,
)

// Portal variant of GammaLoop's generated figure template. The graph source,
// family and optional momentum-arrow treatment live in each asset's `.typ`
// file; layout and particle drawing remain Linnest's responsibility.
#let render(
  input,
  amplitude-mode: false,
  cross-section-mode: false,
  momentum-arrows: false,
  cut-curves: false,
) = {
  set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm), fill: none)
  set text(fill: palette.ink)

  let draw-layout = if cut-curves { layout-with-cut-curves } else { layout }

  context draw-layout(
    input,
    columns: 1,
    unit: 1.5,
    typst-fields: "plain",
    edge-style-options: (
      momentum-arrows: momentum-arrows,
      show-edge-index: false,
      show-momentum: false,
      show-particle: auto,
      label-fill: palette.ink,
      momentum-arrow-stroke: momentum-arrow-defaults.stroke
        + (paint: palette.ink),
    ),
    style-options: (node-label: none),
    amplitude-mode: amplitude-mode,
    cross-section-mode: cross-section-mode,
    additional-data: (
      seed: 42,
    ),
  )
}
