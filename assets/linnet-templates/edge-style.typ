#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, cetz,edge,hide

#let source_stroke(c:black) = (stroke:(paint:c,thickness:1mm))
#let sink_stroke(c:black) = (stroke:source_stroke(c:c).stroke+(dash: (1mm,2mm)))
#let wave = (decorations:cetz.decorations.wave.with(amplitude: 8pt,segment-length:0.4))
#let double = (extrude:(-1mm, 1mm))
#let arrow = (marks:((inherit:"solid",rev:false,pos:1.1,scale:50%),))
#let coil = (decorations:cetz.decorations.coil.with(amplitude: 8pt,segment-length:0.4))
#let zigzag = (decorations:cetz.decorations.zigzag.with(amplitude: 8pt,segment-length:0.4))


#let map = (
  "e+":(source: sink_stroke()+arrow,sink: source_stroke()),
  "e-":(source: source_stroke()+arrow,sink: sink_stroke()),
  "a":(source: source_stroke()+wave,sink:sink_stroke()+wave),
  "d":(source: sink_stroke()+arrow,sink: source_stroke()),
  "d~":(source: source_stroke()+arrow,sink: sink_stroke())
)
