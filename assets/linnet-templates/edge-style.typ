#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, cetz,edge,hide

#let source_stroke(c:black) = (stroke:(paint:c,thickness:0.5mm))
#let sink_stroke(c:black) = (stroke:source_stroke(c:c).stroke+(dash: (0.5mm,1mm)))
#let wave = (decorations:cetz.decorations.wave.with(amplitude: 4pt,segment-length:0.2))
#let double = (extrude:(-0.5mm, 0.5mm))
#let arrow = (marks:((inherit:"solid",rev:false,pos:1.1,scale:50%),))
#let coil = (decorations:cetz.decorations.coil.with(amplitude: 4pt,segment-length:0.2))
#let zigzag = (decorations:cetz.decorations.zigzag.with(amplitude: 4pt,segment-length:0.2))


#let map = (
  "e+":(source: sink_stroke()+arrow,sink: source_stroke()),
  "e-":(source: source_stroke()+arrow,sink: sink_stroke()),
  "a":(source: source_stroke()+wave,sink:sink_stroke()+wave),
  "g":(source: source_stroke()+coil,sink:sink_stroke()+coil),
  "d":(source: sink_stroke()+arrow,sink: source_stroke()),
  "d~":(source: source_stroke()+arrow,sink: sink_stroke())
)
