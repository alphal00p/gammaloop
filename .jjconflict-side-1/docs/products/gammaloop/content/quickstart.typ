#import "../../shared.typ": boundary

#let quickstart = [
= Your first GammaLoop calculation

Generate the same scalar one-loop bubble through the interface that fits your work. Every route
uses GammaLoop's built-in scalar model, disables integrated UV generation, and evaluates the same
three-dimensional loop-momentum point. The result is a compact software check rather than a
normalized collider prediction.

#boundary("Command line", [
  Use the executable for run cards, interactive commands, persistent states, and complete
  calculations. This is the shortest route for most GammaLoop users.

  #link("quickstart/cli/")[Use GammaLoop from the command line →]
])

#boundary("Python", [
  Use the stateful Python API from notebooks, scripts, or workflow automation. It drives the same
  command session and returns structured evaluation records.

  #link("quickstart/python/")[Use GammaLoop from Python →]
])

#boundary("Rust", [
  Embed GammaLoop when a Rust application must own the state lifecycle, issue commands, or retain
  typed evaluation results without crossing a Python boundary.

  #link("quickstart/rust/")[Use GammaLoop from Rust →]
])

If you are unsure, begin with the #link("quickstart/cli/")[command-line guide]. The
#link("reference/interfaces/")[interface guide] explains the longer-term ownership and precision
boundaries once the first calculation works.
]
