#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, cetz,edge,hide
#import "@preview/mitex:0.2.6": *

#let massive = 1mm
#let massless = 0.5mm
#let source_stroke(c:black, thickness:0.5mm,dash:none) = (stroke:(paint:c,thickness:thickness)+dash)
#let sink_stroke(c:black, thickness:0.5mm,dash:none) = (stroke:source_stroke(c:c.lighten(50%), thickness:thickness,dash:dash).stroke)
#let wave = (decorations:cetz.decorations.wave.with(amplitude: 4pt,segment-length:0.2))
#let double = (extrude:(-0.5mm, 0.5mm))
#let arrow = (marks:((inherit:"solid",rev:false,pos:1.1,scale:50%),))
#let antiarrow = (marks:((inherit:"solid",rev:true,pos:1.1,scale:50%),))
#let arrowmap = orientation => if orientation == "Default"{
  arrow
} else if orientation == "Reversed"{
  antiarrow
} else{
  (:)
}
#let coil = (decorations:cetz.decorations.coil.with(amplitude: 4pt,segment-length:0.2))
#let zigzag = (decorations:cetz.decorations.zigzag.with(amplitude: 4pt,segment-length:0.2))
#let dashed = (dash: (0.1em, 0.5em))
#let dotted = (dash: (0.01em, 0.3em))
// Auto-generated particle styles from model (computed in Rust)
#let map = (
  "a": (source:source_stroke(c: black, thickness: massless) + wave, sink:sink_stroke(c: black, thickness: massless) + wave, label:mi(`{\gamma}`)),
  "Z": (source:source_stroke(c: black, thickness: massive) + wave, sink:sink_stroke(c: black, thickness: massive) + wave, label:mi(`{Z}`)),
  "W+": (source:source_stroke(c: blue, thickness: massive) + zigzag, sink:sink_stroke(c: blue, thickness: massive) + zigzag, label:mi(`{W^+}`)),
  "W-": (source:source_stroke(c: blue, thickness: massive) + zigzag, sink:sink_stroke(c: blue, thickness: massive) + zigzag, label:mi(`{W^-}`)),
  "g": (source:source_stroke(c: black, thickness: massless) + coil, sink:sink_stroke(c: black, thickness: massless) + coil, label:mi(`{g}`)),
  "ghA": (source:source_stroke(c: black, thickness: massless,dash: dotted), sink:sink_stroke(c: black, thickness: massless,dash: dotted), label:mi(`\tilde{\gamma}`)),
  "ghA~": (source:source_stroke(c: black, thickness: massless,dash: dotted), sink:sink_stroke(c: black, thickness: massless,dash: dotted), label:mi(`\overline{\tilde{\gamma}}`)),
  "ghZ": (source:source_stroke(c: black, thickness: massive,dash: dotted), sink:sink_stroke(c: black, thickness: massive,dash: dotted), label:mi(`\tilde{Z}`)),
  "ghZ~": (source:source_stroke(c: black, thickness: massive,dash: dotted), sink:sink_stroke(c: black, thickness: massive,dash: dotted), label:mi(`\bar{\tilde{Z}}`)),
  "ghWp": (source:source_stroke(c: blue, thickness: massive,dash: dotted), sink:sink_stroke(c: blue, thickness: massive,dash: dotted), label:mi(`{\tilde{W}^+}`)),
  "ghWp~": (source:source_stroke(c: blue, thickness: massive,dash: dotted), sink:sink_stroke(c: blue, thickness: massive,dash: dotted), label:mi(`{\bar{\tilde{W}}^+}`)),
  "ghWm": (source:source_stroke(c: blue, thickness: massive,dash: dotted), sink:sink_stroke(c: blue, thickness: massive,dash: dotted), label:mi(`{\tilde{W}^-}`)),
  "ghWm~": (source:source_stroke(c: blue, thickness: massive,dash: dotted), sink:sink_stroke(c: blue, thickness: massive,dash: dotted), label:mi(`{\bar{\tilde{W}}^-}`)),
  "ghG": (source:source_stroke(c: black, thickness: massless,dash: dotted), sink:sink_stroke(c: black, thickness: massless,dash: dotted), label:mi(`{\tilde{g}}`)),
  "ghG~": (source:source_stroke(c: black, thickness: massless,dash: dotted), sink:sink_stroke(c: black, thickness: massless,dash: dotted), label:mi(`{\bar{\tilde{g}}}`)),
  "ve": (source:source_stroke(c: black, thickness: massless), sink:sink_stroke(c: black, thickness: massless), label:mi(`{\nu_e}`)),
  "ve~": (source:source_stroke(c: black, thickness: massless), sink:sink_stroke(c: black, thickness: massless), label:mi(`{\overline{\nu}_e}`)),
  "vm": (source:source_stroke(c: black, thickness: massless), sink:sink_stroke(c: black, thickness: massless), label:mi(`{\nu_\mu}`)),
  "vm~": (source:source_stroke(c: black, thickness: massless), sink:sink_stroke(c: black, thickness: massless), label:mi(`{\overline{\nu}_\mu}`)),
  "vt": (source:source_stroke(c: black, thickness: massless), sink:sink_stroke(c: black, thickness: massless), label:mi(`{\nu_\tau}`)),
  "vt~": (source:source_stroke(c: black, thickness: massless), sink:sink_stroke(c: black, thickness: massless), label:mi(`{\overline{\nu}_\tau}`)),
  "u": (source:source_stroke(c: blue, thickness: massless), sink:sink_stroke(c: blue, thickness: massless), label:mi(`{u}`)),
  "u~": (source:source_stroke(c: blue, thickness: massless), sink:sink_stroke(c: blue, thickness: massless), label:mi(`{\overline{u}}`)),
  "c": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{c}`)),
  "c~": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{\overline{c}}`)),
  "t": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{t}`)),
  "t~": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{\overline{t}}`)),
  "d": (source:source_stroke(c: blue, thickness: massless), sink:sink_stroke(c: blue, thickness: massless), label:mi(`{d}`)),
  "d~": (source:source_stroke(c: blue, thickness: massless), sink:sink_stroke(c: blue, thickness: massless), label:mi(`{\overline{d}}`)),
  "s": (source:source_stroke(c: blue, thickness: massless), sink:sink_stroke(c: blue, thickness: massless), label:mi(`{s}`)),
  "s~": (source:source_stroke(c: blue, thickness: massless), sink:sink_stroke(c: blue, thickness: massless), label:mi(`{\overline{s}}`)),
  "b": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{b}`)),
  "b~": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{\overline{b}}`)),
  "H": (source:source_stroke(c: black, thickness: massive,dash: dashed), sink:sink_stroke(c: black, thickness: massive,dash: dashed), label:mi(`{H}`)),
  "G0": (source:source_stroke(c: black, thickness: massive,dash: dashed), sink:sink_stroke(c: black, thickness: massive,dash: dashed), label:mi(`{G_0}`)),
  "G+": (source:source_stroke(c: blue, thickness: massive,dash: dashed), sink:sink_stroke(c: blue, thickness: massive,dash: dashed), label:mi(`{G^+}`)),
  "G-": (source:source_stroke(c: blue, thickness: massive,dash: dashed), sink:sink_stroke(c: blue, thickness: massive,dash: dashed), label:mi(`{G^-}`)),
  "e-": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{e^-}`)),
  "e+": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{e^+}`)),
  "mu-": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{\mu^-}`)),
  "mu+": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{\mu^+}`)),
  "ta-": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{\tau^-}`)),
  "ta+": (source:source_stroke(c: blue, thickness: massive), sink:sink_stroke(c: blue, thickness: massive), label:mi(`{\tau^+}`)),
)
