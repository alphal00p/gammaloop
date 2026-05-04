#let input-int(name, default) = int(sys.inputs.at(name, default: str(default)))
#let input-float(name, default) = float(sys.inputs.at(name, default: str(default)))

#let count = input-int("count", 250)
#let turns = input-int("turns", 18)
#let samples-per-period = input-int("samples-per-period", 16)
#let width = input-float("width", 9.0)
#let spacing = input-float("spacing", 0.24)
#let amplitude = input-float("amplitude", 0.12)
#let wavelength = width / turns
#let stroke = black + 0.35pt

#let straight-segment(y) = (
  start: (x: 0, y: y),
  ctrl_a: (x: width / 3, y: y),
  ctrl_b: (x: width * 2 / 3, y: y),
  end: (x: width, y: y),
)
