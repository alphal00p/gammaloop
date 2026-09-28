#import "../theme.typ": palette
#import "../../../../crates/kurvst/typst/examples/knot-logo.typ": logo

#set page(height: auto, width: auto, margin: 10mm, fill: none)

// GammaLoop mark: the Kurvst knot logo, whose over/under cuts are geometry
// rather than strokes painted in the page colour.
#logo(palette.ink)
