#import "../src/lib.typ": draw, graph, layouts, subgraph
#import "@preview/cetz:0.5.1" as cetz

#set page(width: auto, height: auto)
#set text(size: 7pt)

#let diagram-unit = 2.0
#let auto-root-node-label-style = (padding: 0.62)
#let auto-root-node-corner-radius = 0.08


#let display-statement(value, fallback: "") = {
  let value = str(value).replace("\\\"", "\"").trim("\"")
  if value == "" { fallback } else { value }
}

#let basketball-src = ```dot
digraph basketball {
  node [num = "1"]
  edge [particle=scalar_1]
  e [style=invis]

  e -> A:1 [id=5]
  B:0 -> e [id=4]
  A -> B [id=0 lmb_id=0]
  A -> B [id=1 lmb_id=1]
  A -> B [id=2 lmb_id=2]
  A -> B [id=3]
}
```

#let basketball-node-style(node) = (
  radius: 0.14,
  fill: white,
  stroke: rgb("#60656f") + 0.28pt,
)

#let basketball-edge-style(edge) = (
  stroke: rgb("#9aa6b6") + 0.28pt,
  route: "edge-pos",
)

#let basketball = graph.parse(basketball-src.text).first()
#let basketball = layouts.layout(
  basketball,
  layout-algo: "force",
  epochs: 52,
  steps: 60,
  seed: 7,
  viewport-w: 4.0,
  viewport-h: 2.8,
  length-scale: 0.38,
  beta: 80.0,
  k-spring: 12.0,
  gamma-dangling: .1,
  gamma-ee: 0.25,
  gamma-ev: 0.1,
  g-center: 0.01,
  label-steps: 0,
)
#let g = basketball


#let S(g, value, ..rest) = subgraph.label(g, display-statement(value))
#let basketball-highlight-subgraphs(approx, given) = {
  let highlights = ()
  if approx != none and approx != 0 {
    highlights.push((subgraph: approx, edge-style: (stroke: rgb("#3b82f6").transparentize(50%) + 1.45pt)))
  }
  if given != none and given != 0 {
    highlights.push((subgraph: given, edge-style: (stroke: rgb("#e4504f").transparentize(0%) + 1.65pt)))
  }
  highlights
}
#let T(approx, given, expr) = {
  let mini = box(inset: -2mm, scale(42%, draw(
    g,
    unit: 1.0,
    subgraph: basketball-highlight-subgraphs(approx, given),
    node-label: none,
    node-style: basketball-node-style,
    source-style: basketball-edge-style,
    sink-style: basketball-edge-style,
    edge-label: none,
    node-min-radius: 0.12,
    node-label-padding: 0.04,
    padding: 0.12,
    subgraph-edge-style: (stroke: rgb("#e4504f") + 1.45pt),
  )))
  $op("T")_(#mini)(#expr)$
}

#let label-scope = (
  g: g,
  S: S,
  T: T,
  expr: "expr",
)

#let auto-root-node-label(node) = {
  let label = display-statement(node.at("label", default: str(node.at("vid", default: ""))), fallback: "")
  text(size: 0.82em)[#eval(label, mode: "code", scope: label-scope)]
}

#let auto-root-node-style(node) = {
  let cover = display-statement(node.at("cover", default: ""), fallback: "")
  (
    fill: if cover == "0" {
      rgb("#fff2b8")
    } else if cover == "GS" {
      rgb("#dff0ff")
    } else {
      white
    },
    stroke: rgb("#60656f") + 0.28pt,
  )
}

#let draw-auto-root-node(node, box) = {
  let center = box.center
  let half-width = box.width / 2
  let half-height = box.height / 2
  cetz.draw.rect(
    (center.at(0) - half-width, center.at(1) - half-height),
    (center.at(0) + half-width, center.at(1) + half-height),
    radius: auto-root-node-corner-radius,
    name: box.name,
    fill: box.style.fill,
    stroke: box.style.stroke,
  )
  if box.label != none {
    cetz.draw.content(center, box.label, padding: 0, ..box.label-style)
  }
}

#let auto-root-source-edge-style(edge) = (
  stroke: rgb("#7890b4") + 0.45pt,
  route: "direct",
  source-anchor: "south",
)

#let auto-root-sink-edge-style(edge) = (
  stroke: rgb("#7890b4") + 0.45pt,
  route: "direct",
  sink-anchor: "north",
)

#let draw_auto_root_dot(doc) = {
  show raw: it => {
    if it.at("lang") != "dot" {
      return it
    }
    for g in graph.parse(it.text) {
      let g = graph.style(
        g,
        node-label: auto-root-node-label,
        node-style: auto-root-node-style,
        node-label-style: auto-root-node-label-style,
        unit: diagram-unit,
      )
      let g = layouts.layout(
        g,
        layout-algo: "dot",
        tree-dx: 0.42,
        tree-dy: 7.0,
        route-label-width-cap: 0.,
        label-steps: 0,
      )
      [#metadata(graph.nodes(g)) <linnest-auto-root-nodes>]
      draw(
        g,
        unit: diagram-unit,
        draw-node: draw-auto-root-node,
        node-label: auto-root-node-label,
        node-style: auto-root-node-style,
        node-label-style: auto-root-node-label-style,
        source-style: auto-root-source-edge-style,
        sink-style: auto-root-sink-edge-style,
        edge-label: none,
        edge-omega: 0.35,
        node-min-radius: 0.2,
        node-label-padding: 0.1,
        padding: 0.45,
      )
    }
  }
  doc
}

#show: draw_auto_root_dot

```dot
digraph {
  node   [shape=circle,height=0.1,label=""];
  overlap = "scale";
  layout = "neato";
  start=2;
  0   [label="\"expr\"" foata="" cover="0"];
  1   [label="T(S(g,\"3I\",0),0,\"expr\")" foata="3I" cover="3I"];
  2   [label="T(S(g,\"3s\",1),0,\"expr\")" foata="3s" cover="3s"];
  3   [label="T(S(g,\"Ca\",2),0,\"expr\")" foata="Ca" cover="Ca"];
  4   [label="T(S(g,\"Fg\",3),0,\"expr\")" foata="Fg" cover="Fg"];
  5   [label="T(S(g,\"44\",4),0,\"expr\")" foata="44" cover="44"];
  6   [label="T(S(g,\"DA\",5),0,\"expr\")" foata="DA" cover="DA"];
  7   [label="T(S(g,\"GG\",6),0,\"expr\")" foata="GG" cover="GG"];
  8   [label="T(S(g,\"y\",7),0,\"expr\")" foata="y" cover="y"];
  9   [label="T(S(g,\"DM\",8),0,\"expr\")" foata="DM" cover="DM"];
  10  [label="T(S(g,\"FU\",9),0,\"expr\")" foata="FU" cover="FU"];
  11  [label="T(S(g,\"GS\",10),0,\"expr\")" foata="GS" cover="GS"];
  12  [label="T(S(g,\"Fg\",11),S(g,\"3I\"),T(S(g,\"3I\",0),0,\"expr\"))" foata="3I;Fg" cover="Fg"];
  13  [label="T(S(g,\"44\",12),S(g,\"3I\"),T(S(g,\"3I\",0),0,\"expr\"))" foata="3I;44" cover="44"];
  14  [label="T(S(g,\"GS\",13),S(g,\"3I\"),T(S(g,\"3I\",0),0,\"expr\"))" foata="3I;GS" cover="GS"];
  15  [label="T(S(g,\"44\",14),S(g,\"3s\"),T(S(g,\"3s\",1),0,\"expr\"))" foata="3s;44" cover="44"];
  16  [label="T(S(g,\"GG\",15),S(g,\"3s\"),T(S(g,\"3s\",1),0,\"expr\"))" foata="3s;GG" cover="GG"];
  17  [label="T(S(g,\"GS\",16),S(g,\"3s\"),T(S(g,\"3s\",1),0,\"expr\"))" foata="3s;GS" cover="GS"];
  18  [label="T(S(g,\"Fg\",17),S(g,\"Ca\"),T(S(g,\"Ca\",2),0,\"expr\"))" foata="Ca;Fg" cover="Fg"];
  19  [label="T(S(g,\"DM\",18),S(g,\"Ca\"),T(S(g,\"Ca\",2),0,\"expr\"))" foata="Ca;DM" cover="DM"];
  20  [label="T(S(g,\"GS\",19),S(g,\"Ca\"),T(S(g,\"Ca\",2),0,\"expr\"))" foata="Ca;GS" cover="GS"];
  21  [label="T(S(g,\"GS\",20),S(g,\"Fg\"),T(S(g,\"Fg\",3),0,\"expr\"))" foata="Fg;GS" cover="GS"];
  22  [label="T(S(g,\"GS\",21),S(g,\"44\"),T(S(g,\"44\",4),0,\"expr\"))" foata="44;GS" cover="GS"];
  23  [label="T(S(g,\"GG\",22),S(g,\"DA\"),T(S(g,\"DA\",5),0,\"expr\"))" foata="DA;GG" cover="GG"];
  24  [label="T(S(g,\"DM\",23),S(g,\"DA\"),T(S(g,\"DA\",5),0,\"expr\"))" foata="DA;DM" cover="DM"];
  25  [label="T(S(g,\"GS\",24),S(g,\"DA\"),T(S(g,\"DA\",5),0,\"expr\"))" foata="DA;GS" cover="GS"];
  26  [label="T(S(g,\"GS\",25),S(g,\"GG\"),T(S(g,\"GG\",6),0,\"expr\"))" foata="GG;GS" cover="GS"];
  27  [label="T(S(g,\"44\",26),S(g,\"y\"),T(S(g,\"y\",7),0,\"expr\"))" foata="y;44" cover="44"];
  28  [label="T(S(g,\"DM\",27),S(g,\"y\"),T(S(g,\"y\",7),0,\"expr\"))" foata="y;DM" cover="DM"];
  29  [label="T(S(g,\"GS\",28),S(g,\"y\"),T(S(g,\"y\",7),0,\"expr\"))" foata="y;GS" cover="GS"];
  30  [label="T(S(g,\"GS\",29),S(g,\"DM\"),T(S(g,\"DM\",8),0,\"expr\"))" foata="DM;GS" cover="GS"];
  31  [label="T(S(g,\"Fg\",30),S(g,\"FU\"),T(S(g,\"FU\",9),0,\"expr\"))" foata="FU;Fg" cover="Fg"];
  32  [label="T(S(g,\"GG\",31),S(g,\"FU\"),T(S(g,\"FU\",9),0,\"expr\"))" foata="FU;GG" cover="GG"];
  33  [label="T(S(g,\"GS\",32),S(g,\"FU\"),T(S(g,\"FU\",9),0,\"expr\"))" foata="FU;GS" cover="GS"];
  34  [label="T(S(g,\"GS\",20),S(g,\"Fg\"),T(S(g,\"Fg\",11),S(g,\"3I\"),T(S(g,\"3I\",0),0,\"expr\")))" foata="3I;Fg;GS" cover="GS"];
  35  [label="T(S(g,\"GS\",21),S(g,\"44\"),T(S(g,\"44\",12),S(g,\"3I\"),T(S(g,\"3I\",0),0,\"expr\")))" foata="3I;44;GS" cover="GS"];
  36  [label="T(S(g,\"GS\",21),S(g,\"44\"),T(S(g,\"44\",14),S(g,\"3s\"),T(S(g,\"3s\",1),0,\"expr\")))" foata="3s;44;GS" cover="GS"];
  37  [label="T(S(g,\"GS\",25),S(g,\"GG\"),T(S(g,\"GG\",15),S(g,\"3s\"),T(S(g,\"3s\",1),0,\"expr\")))" foata="3s;GG;GS" cover="GS"];
  38  [label="T(S(g,\"GS\",20),S(g,\"Fg\"),T(S(g,\"Fg\",17),S(g,\"Ca\"),T(S(g,\"Ca\",2),0,\"expr\")))" foata="Ca;Fg;GS" cover="GS"];
  39  [label="T(S(g,\"GS\",29),S(g,\"DM\"),T(S(g,\"DM\",18),S(g,\"Ca\"),T(S(g,\"Ca\",2),0,\"expr\")))" foata="Ca;DM;GS" cover="GS"];
  40  [label="T(S(g,\"GS\",25),S(g,\"GG\"),T(S(g,\"GG\",22),S(g,\"DA\"),T(S(g,\"DA\",5),0,\"expr\")))" foata="DA;GG;GS" cover="GS"];
  41  [label="T(S(g,\"GS\",29),S(g,\"DM\"),T(S(g,\"DM\",23),S(g,\"DA\"),T(S(g,\"DA\",5),0,\"expr\")))" foata="DA;DM;GS" cover="GS"];
  42  [label="T(S(g,\"GS\",21),S(g,\"44\"),T(S(g,\"44\",26),S(g,\"y\"),T(S(g,\"y\",7),0,\"expr\")))" foata="y;44;GS" cover="GS"];
  43  [label="T(S(g,\"GS\",29),S(g,\"DM\"),T(S(g,\"DM\",27),S(g,\"y\"),T(S(g,\"y\",7),0,\"expr\")))" foata="y;DM;GS" cover="GS"];
  44  [label="T(S(g,\"GS\",20),S(g,\"Fg\"),T(S(g,\"Fg\",30),S(g,\"FU\"),T(S(g,\"FU\",9),0,\"expr\")))" foata="FU;Fg;GS" cover="GS"];
  45  [label="T(S(g,\"GS\",25),S(g,\"GG\"),T(S(g,\"GG\",31),S(g,\"FU\"),T(S(g,\"FU\",9),0,\"expr\")))" foata="FU;GG;GS" cover="GS"];
  0:0:s  -> 1:1:s  [id=0  color="red:blue;0.5"];
  0:2:s  -> 2:3:s  [id=1  color="red:blue;0.5"];
  0:4:s  -> 3:5:s  [id=2  color="red:blue;0.5"];
  0:6:s  -> 4:7:s  [id=3  color="red:blue;0.5"];
  0:8:s  -> 5:9:s  [id=4  color="red:blue;0.5"];
  0:10:s -> 6:11:s [id=5  color="red:blue;0.5"];
  0:12:s -> 7:13:s [id=6  color="red:blue;0.5"];
  0:14:s -> 8:15:s [id=7  color="red:blue;0.5"];
  0:16:s -> 9:17:s [id=8  color="red:blue;0.5"];
  0:18:s -> 10:19:s [id=9  color="red:blue;0.5"];
  0:20:s -> 11:21:s [id=10  color="red:blue;0.5"];
  1:22:s -> 12:23:s [id=11  color="red:blue;0.5"];
  1:24:s -> 13:25:s [id=12  color="red:blue;0.5"];
  1:26:s -> 14:27:s [id=13  color="red:blue;0.5"];
  2:28:s -> 15:29:s [id=14  color="red:blue;0.5"];
  2:30:s -> 16:31:s [id=15  color="red:blue;0.5"];
  2:32:s -> 17:33:s [id=16  color="red:blue;0.5"];
  3:34:s -> 18:35:s [id=17  color="red:blue;0.5"];
  3:36:s -> 19:37:s [id=18  color="red:blue;0.5"];
  3:38:s -> 20:39:s [id=19  color="red:blue;0.5"];
  4:40:s -> 21:41:s [id=20  color="red:blue;0.5"];
  5:42:s -> 22:43:s [id=21  color="red:blue;0.5"];
  6:44:s -> 23:45:s [id=22  color="red:blue;0.5"];
  6:46:s -> 24:47:s [id=23  color="red:blue;0.5"];
  6:48:s -> 25:49:s [id=24  color="red:blue;0.5"];
  7:50:s -> 26:51:s [id=25  color="red:blue;0.5"];
  8:52:s -> 27:53:s [id=26  color="red:blue;0.5"];
  8:54:s -> 28:55:s [id=27  color="red:blue;0.5"];
  8:56:s -> 29:57:s [id=28  color="red:blue;0.5"];
  9:58:s -> 30:59:s [id=29  color="red:blue;0.5"];
  10:60:s -> 31:61:s [id=30  color="red:blue;0.5"];
  10:62:s -> 32:63:s [id=31  color="red:blue;0.5"];
  10:64:s -> 33:65:s [id=32  color="red:blue;0.5"];
  12:66:s -> 34:67:s [id=33  color="red:blue;0.5"];
  13:68:s -> 35:69:s [id=34  color="red:blue;0.5"];
  15:70:s -> 36:71:s [id=35  color="red:blue;0.5"];
  16:72:s -> 37:73:s [id=36  color="red:blue;0.5"];
  18:74:s -> 38:75:s [id=37  color="red:blue;0.5"];
  19:76:s -> 39:77:s [id=38  color="red:blue;0.5"];
  23:78:s -> 40:79:s [id=39  color="red:blue;0.5"];
  24:80:s -> 41:81:s [id=40  color="red:blue;0.5"];
  27:82:s -> 42:83:s [id=41  color="red:blue;0.5"];
  28:84:s -> 43:85:s [id=42  color="red:blue;0.5"];
  31:86:s -> 44:87:s [id=43  color="red:blue;0.5"];
  32:88:s -> 45:89:s [id=44  color="red:blue;0.5"];
}
```
