// Shared website palette. All SVG sources read the same `theme` input so the
// light and dark assets cannot drift independently.
#let theme = str(sys.inputs.at("theme", default: "light")).trim("\"")

#let palette = if theme == "light" {
  (
    surface: rgb("#f9f6f0"),
    ink: rgb("#3d2645"),
    accent: rgb("#6f4d85"),
    soft: rgb("#d8b9e3"),
    cut-blue: rgb("#0072b2"),
    cut-red: rgb("#d94f3d"),
  )
} else if theme == "dark" {
  (
    surface: rgb("#211a23"),
    ink: rgb("#f8effa"),
    accent: rgb("#d8b9e3"),
    soft: rgb("#6f4d85"),
    cut-blue: rgb("#62b7e8"),
    cut-red: rgb("#ff8372"),
  )
} else {
  panic("website SVG theme must be light or dark")
}
