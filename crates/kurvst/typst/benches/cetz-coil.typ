#import "@preview/cetz:0.3.4" as cetz
#import "common.typ": count, turns, width, spacing, amplitude, stroke

#set page(width: 120mm, height: auto, margin: 8mm)

#cetz.canvas({
  for index in range(count) {
    let y = -index * spacing
    cetz.decorations.coil(
      cetz.draw.line((0, y), (width, y), stroke: none, mark: none),
      segments: turns,
      amplitude: amplitude,
      factor: 150%,
      stroke: stroke,
      mark: none,
    )
  }
})
