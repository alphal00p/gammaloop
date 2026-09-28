#import "../theme.typ": palette
#import "../../../../crates/kurvst/typst/examples/knot-logo.typ": construction

#set page(height: auto, width: auto, margin: 10mm, fill: none)

// How the GammaLoop mark is built: each letter's outline, the cut strips around
// the over stroke, and the centerlines, for the developer documentation.
#construction(
  colors: (palette.cut-blue, palette.accent, palette.soft, palette.cut-blue),
  cut: palette.cut-red,
  line: palette.ink,
  tint: 50%,
)
