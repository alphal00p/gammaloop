#import "physics-edge-style.typ": mi, palette, massive, massless, dashed, dotted, stroke-style, source-stroke, sink-stroke, fermion-flow, wave, coil, zigzag, default-edge, default-map, style as physics-style

// Auto-generated particle styles from the canonical model. The reusable
// drawing callbacks live in physics-edge-style.typ.
#let generated-map = (
  "a": (source:source-stroke(c: palette.ink, thickness: massless) + wave, sink:sink-stroke(c: palette.ink, thickness: massless) + wave, label:mi(`{\gamma}`)),
  "Z": (source:source-stroke(c: palette.ink, thickness: massive) + wave, sink:sink-stroke(c: palette.ink, thickness: massive) + wave, label:mi(`{Z}`)),
  "W+": (source:source-stroke(c: palette.accent, thickness: massive) + zigzag, sink:sink-stroke(c: palette.accent, thickness: massive) + zigzag, label:mi(`{W^+}`)),
  "W-": (source:source-stroke(c: palette.accent, thickness: massive) + zigzag, sink:sink-stroke(c: palette.accent, thickness: massive) + zigzag, label:mi(`{W^-}`)),
  "g": (source:source-stroke(c: palette.ink, thickness: massless) + coil, sink:sink-stroke(c: palette.ink, thickness: massless) + coil, label:mi(`{g}`)),
  "ghA": (source:source-stroke(c: palette.ink, thickness: massless, dash: dotted), sink:sink-stroke(c: palette.ink, thickness: massless, dash: dotted), label:mi(`\tilde{\gamma}`)),
  "ghA~": (source:source-stroke(c: palette.ink, thickness: massless, dash: dotted), sink:sink-stroke(c: palette.ink, thickness: massless, dash: dotted), label:mi(`\overline{\tilde{\gamma}}`)),
  "ghZ": (source:source-stroke(c: palette.ink, thickness: massive, dash: dotted), sink:sink-stroke(c: palette.ink, thickness: massive, dash: dotted), label:mi(`\tilde{Z}`)),
  "ghZ~": (source:source-stroke(c: palette.ink, thickness: massive, dash: dotted), sink:sink-stroke(c: palette.ink, thickness: massive, dash: dotted), label:mi(`\bar{\tilde{Z}}`)),
  "ghWp": (source:source-stroke(c: palette.accent, thickness: massive, dash: dotted), sink:sink-stroke(c: palette.accent, thickness: massive, dash: dotted), label:mi(`{\tilde{W}^+}`)),
  "ghWp~": (source:source-stroke(c: palette.accent, thickness: massive, dash: dotted), sink:sink-stroke(c: palette.accent, thickness: massive, dash: dotted), label:mi(`{\bar{\tilde{W}}^+}`)),
  "ghWm": (source:source-stroke(c: palette.accent, thickness: massive, dash: dotted), sink:sink-stroke(c: palette.accent, thickness: massive, dash: dotted), label:mi(`{\tilde{W}^-}`)),
  "ghWm~": (source:source-stroke(c: palette.accent, thickness: massive, dash: dotted), sink:sink-stroke(c: palette.accent, thickness: massive, dash: dotted), label:mi(`{\bar{\tilde{W}}^-}`)),
  "ghG": (source:source-stroke(c: palette.ink, thickness: massless, dash: dotted), sink:sink-stroke(c: palette.ink, thickness: massless, dash: dotted), label:mi(`{\tilde{g}}`)),
  "ghG~": (source:source-stroke(c: palette.ink, thickness: massless, dash: dotted), sink:sink-stroke(c: palette.ink, thickness: massless, dash: dotted), label:mi(`{\bar{\tilde{g}}}`)),
  "ve": (source:source-stroke(c: palette.ink, thickness: massless), sink:sink-stroke(c: palette.ink, thickness: massless), label:mi(`{\nu_e}`)) + fermion-flow,
  "ve~": (source:source-stroke(c: palette.ink, thickness: massless), sink:sink-stroke(c: palette.ink, thickness: massless), label:mi(`{\overline{\nu}_e}`)) + fermion-flow,
  "vm": (source:source-stroke(c: palette.ink, thickness: massless), sink:sink-stroke(c: palette.ink, thickness: massless), label:mi(`{\nu_\mu}`)) + fermion-flow,
  "vm~": (source:source-stroke(c: palette.ink, thickness: massless), sink:sink-stroke(c: palette.ink, thickness: massless), label:mi(`{\overline{\nu}_\mu}`)) + fermion-flow,
  "vt": (source:source-stroke(c: palette.ink, thickness: massless), sink:sink-stroke(c: palette.ink, thickness: massless), label:mi(`{\nu_\tau}`)) + fermion-flow,
  "vt~": (source:source-stroke(c: palette.ink, thickness: massless), sink:sink-stroke(c: palette.ink, thickness: massless), label:mi(`{\overline{\nu}_\tau}`)) + fermion-flow,
  "u": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{u}`)) + fermion-flow,
  "u~": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{\overline{u}}`)) + fermion-flow,
  "c": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{c}`)) + fermion-flow,
  "c~": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{\overline{c}}`)) + fermion-flow,
  "t": (source:source-stroke(c: palette.accent, thickness: massive), sink:sink-stroke(c: palette.accent, thickness: massive), label:mi(`{t}`)) + fermion-flow,
  "t~": (source:source-stroke(c: palette.accent, thickness: massive), sink:sink-stroke(c: palette.accent, thickness: massive), label:mi(`{\overline{t}}`)) + fermion-flow,
  "d": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{d}`)) + fermion-flow,
  "d~": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{\overline{d}}`)) + fermion-flow,
  "s": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{s}`)) + fermion-flow,
  "s~": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{\overline{s}}`)) + fermion-flow,
  "b": (source:source-stroke(c: palette.accent, thickness: massive), sink:sink-stroke(c: palette.accent, thickness: massive), label:mi(`{b}`)) + fermion-flow,
  "b~": (source:source-stroke(c: palette.accent, thickness: massive), sink:sink-stroke(c: palette.accent, thickness: massive), label:mi(`{\overline{b}}`)) + fermion-flow,
  "H": (source:source-stroke(c: palette.ink, thickness: massive, dash: dashed), sink:sink-stroke(c: palette.ink, thickness: massive, dash: dashed), label:mi(`{H}`)),
  "G0": (source:source-stroke(c: palette.ink, thickness: massive, dash: dashed), sink:sink-stroke(c: palette.ink, thickness: massive, dash: dashed), label:mi(`{G_0}`)),
  "G+": (source:source-stroke(c: palette.accent, thickness: massive, dash: dashed), sink:sink-stroke(c: palette.accent, thickness: massive, dash: dashed), label:mi(`{G^+}`)),
  "G-": (source:source-stroke(c: palette.accent, thickness: massive, dash: dashed), sink:sink-stroke(c: palette.accent, thickness: massive, dash: dashed), label:mi(`{G^-}`)),
  "e-": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{e^-}`)) + fermion-flow,
  "e+": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{e^+}`)) + fermion-flow,
  "mu-": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{\mu^-}`)) + fermion-flow,
  "mu+": (source:source-stroke(c: palette.accent, thickness: massless), sink:sink-stroke(c: palette.accent, thickness: massless), label:mi(`{\mu^+}`)) + fermion-flow,
  "ta-": (source:source-stroke(c: palette.accent, thickness: massive), sink:sink-stroke(c: palette.accent, thickness: massive), label:mi(`{\tau^-}`)) + fermion-flow,
  "ta+": (source:source-stroke(c: palette.accent, thickness: massive), sink:sink-stroke(c: palette.accent, thickness: massive), label:mi(`{\tau^+}`)) + fermion-flow,
)

#let map = generated-map

// Merge model-specific particles with neutral aliases and caller overrides.
#let style(map: (:), typst-fields: "plain", ..options) = physics-style(
  map: default-map + generated-map + map,
  typst-fields: typst-fields,
  ..options.named(),
)

#let source-style(edge, typst-fields: "plain", ..options) = {
  let callbacks = style(typst-fields: typst-fields, ..options.named())
  (callbacks.source-style)(edge)
}

#let sink-style(edge, typst-fields: "plain", ..options) = {
  let callbacks = style(typst-fields: typst-fields, ..options.named())
  (callbacks.sink-style)(edge)
}

#let edge-label(edge, typst-fields: "plain", ..options) = {
  let callbacks = style(typst-fields: typst-fields, ..options.named())
  (callbacks.edge-label)(edge)
}
