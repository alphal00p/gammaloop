#import "layout.typ": layout, layout-with-cut-curves
#import "../theme.typ": palette

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
      show-particle: false,
      momentum-arrow-stroke: (
        paint: palette.ink,
        thickness: 1.0pt,
        cap: "round",
      ),
      momentum-arrow-mark: (end: "straight", scale: 1.1),
    ),
    amplitude-mode: amplitude-mode,
    cross-section-mode: cross-section-mode,
    additional-data: (
      steps: 1200,
      seed: 42,
    ),
  )
}
