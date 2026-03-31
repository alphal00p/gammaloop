#import "@preview/fletcher:0.5.8" as fletcher: cetz, diagram, edge, hide, node
#let mom_arr = (
  stroke: black + 0.3mm,
  marks: ((inherit: "head", rev: false, pos: 1, scale: 40%),),
)
#let p = plugin("./linnest.wasm")

#let eval-dict(data, name, scope) = {
  let dict = data.remove(name, default: "(:)")
  if dict == none {
    dict = "(:)"
  }
  eval(dict, scope: scope, mode: "code")
}
#let graph-info(graph) = cbor(p.graph_info(graph))
#let graph-nodes(graph, subgraph: none) = if subgraph == none {
  cbor(p.graph_nodes(graph))
} else {
  cbor(p.graph_nodes_of(graph, cbor.encode(subgraph)))
}
#let graph-edges(graph, subgraph: none) = if subgraph == none {
  cbor(p.graph_edges(graph))
} else {
  cbor(p.graph_edges_of(graph, cbor.encode(subgraph)))
}
#let layout(
  input,
  split_edge: true,
  scope: (:),
  columns: 1fr,
  unit: 1,
  additional_data: (:),
) = {
  let parsed = p.parse_graph(bytes(input))
  let config = cbor.encode(
    (
      steps: sys.inputs.at("steps", default: "15"),
      seed: sys.inputs.at("seed", default: "14"),
      step: ".81",
      step_shrink: "0.21",
      temp: ".3",
      beta: "46.1",
      k_spring: "11.",
      g_center: "40.0",
      epochs:"30",
      crossing_penalty: "30",
      gamma_dangling: "40",
      gamma_ee: ".1",
      directional_force:"5",
      label_length_scale:".6",
      label_spring:"23",
      label_charge:"3",
      label_steps:"20",
      gamma_ev: ".1",

      z_spring:"0.05",
      z_spring_growth:"1.3",
      length_scale: "0.1",
    ) + additional_data,
  )
  let graphs = cbor(parsed)
  let diags = ()
  for graph in graphs {
    let graph = p.layout_parsed_graph(bytes(graph), config)
    let g = graph-info(graph)
    let nodes = graph-nodes(graph)
    let edges = graph-edges(graph)
    let noed = ()
    let n = (:)


    for (i, v) in nodes.enumerate() {
      n.insert(str(i), v)
      let pos = v.pos
      let x = pos.x
      let y = pos.y
      let ev = v.eval
      if ev == none {
        ev = "(:)"
      }
      noed.push(node(
        pos: (x * unit, y * unit),
        name: label(str(i)),
        ..eval(ev, scope: scope+(vid:i), mode: "code"),
        layer: 2,
      ))
    }


    for (i, e) in edges.enumerate() {




      let start = e.source
      let end = e.sink
      let source = if start == none { none } else { start.statement }
      let sink = if end == none { none } else { end.statement }

      let ext = start == none or end == none;
      let data = e.statements

      let o = e.orientation
      let ev_sink = eval-dict(data, "eval_sink",scope+(orientation:o)+(eid:i)+(ext:ext))
      let ev_label = eval-dict(data,"eval_label",scope+(orientation:o)+(eid:i)+(ext:ext)+(sink:sink)+(source:source)+data)
      // let ev_label =[in]



      // ev_label
      let ev_source = eval-dict(data, "eval_source",scope+(orientation:o)+(eid:i))

      let bend = e.bend
      if bend == none {
        bend = 0.
      }

      let enmlab = label("em" + str(i))
      let (end-node, end-node-pos) = if end != none {
        let nodelab = label(str(end.node))
        if start != none and nodelab == label(str(start.node)){
          bend = bend + 2
        }

          noed.push(edge(
            vertices: ((e.pos.x * unit, e.pos.y * unit), nodelab),
            bend: bend * 0.5rad,
            ..ev_sink,
          ))


        (nodelab, n.at(str(end.node)).pos)
      } else {
        let lab = label("exte" + str(i))
        noed.push(node(
          (e.pos.x * unit, e.pos.y * unit),
          name: lab,
          outset: -5mm,
          radius: 5mm,
          fill: none,
        ))
        (lab, e.pos)
      }

      let snmlab = label("sm" + str(i))
      let (start-node, start-node-pos) = if start != none {
        let nodelab = label(str(start.node))



          noed.push(edge(
            vertices: (nodelab, (e.pos.x * unit, e.pos.y * unit)),
            bend: bend * 0.5rad,
            ..ev_source,
          ))


        (nodelab, n.at(str(start.node)).pos)
      } else {
        let lab = label("exts" + str(i))
        noed.push(node(
          (e.pos.x * unit, e.pos.y * unit),
          name: lab,
          outset: -5mm,
          radius: 5mm,
          fill: none,
        ))
        (lab, e.pos)
      }




      let percentb = 1 + calc.abs(bend / calc.pi)

      let a = (
        calc.sqrt(
          calc.pow(start-node-pos.x - end-node-pos.x, 2)
            + calc.pow(start-node-pos.y - end-node-pos.y, 2),
        )
          * 2.5
          * percentb
          * unit
      )

      noed.push(node(
        pos: (start-node-pos.x * unit, start-node-pos.y * unit),
        name: snmlab,
        outset: a,
      ))
      noed.push(node(
        pos: (end-node-pos.x * unit, end-node-pos.y * unit),
        name: enmlab,
        outset: a,
      ))
      let label-pos = if e.label_pos == none { e.pos } else { e.label_pos }
      noed.push(node(
          (label-pos.x * unit, label-pos.y * unit),
          ev_label,//rotate(data.label_angle*-1rad,ev_label)
          inset:0mm,
          snap:false,
          name: label("e"+str(i)),
          fill: none,
        ))

      let shift = (
        if bend != 0. {
          bend * 0.15
        } else { 1 }
          * 1.5mm
      )


      // noed.push(edge(vertices:(snmlab,enmlab),bend:bend * (percentb -1) * -1rad,shift:shift,..mom_arr,..eval(mev,scope: scope ,mode: "code")))
    }
    diags.push(grid(
      align: center + top,
      gutter: 1em,
      [#g.name],
      diagram(
        debug:0,
        node-shape: circle,
        node-fill: black,
        edge-stroke: 0.1em,
        spacing: 2em,
        ..noed,
      ),
      // {set align(left);raw(parse)},
    ))
  }
  for d in diags{
    d
  }
  // grid(
  //   align: center + top, gutter: 1em,
  //   columns: columns,
  //   ..diags
  // )
}
