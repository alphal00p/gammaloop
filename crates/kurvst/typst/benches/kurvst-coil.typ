#import "@preview/cetz:0.3.4" as cetz
#import "../src/lib.typ" as kurvst
#import "common.typ": count, samples-per-period, spacing, amplitude, wavelength, stroke, straight-segment

#set page(width: 120mm, height: auto, margin: 8mm)

#cetz.canvas({
  for index in range(count) {
    let path = kurvst.pattern-segment(
      straight-segment(-index * spacing),
      pattern: kurvst.coil(samples-per-period: samples-per-period, longitudinal-scale: 1.25),
      amplitude: amplitude,
      wavelength: wavelength,
      samples-per-period: samples-per-period,
      coil-longitudinal-scale: 1.25,
    )
    kurvst.cetz-pattern(path, stroke: stroke)
  }
})
