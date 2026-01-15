// Default per-graph Typst template bundled with the linnet CLI.
// The CLI passes --input title="..." and --input data="..." for every build.
#{
  import "layout.typ": layout
  import "edge-style.typ": *
  show raw: it => if it.at("lang") == "dot" {
    layout(
      it.at("text"),
      unit: 10mm,
      columns: 1,
      scope: (map: map,arrowmap:arrowmap),
      additional_data: (
                node: (
                  eval: "(stroke:blue,fill :white,
                  radius:6pt,
                  outset: -1pt,label:[#set text(size:3pt);{int_id}])",
                ),
                edge: (
                  eval_source: "map.at(\"{particle}\").source+arrowmap(orientation)+(label:map.at(\"{particle}\").at(\"label\",default:[a])+[ ]+[$e_#eid$])",
                  eval_sink: "map.at(\"{particle}\").sink+(label:map.at(\"{particle}\").at(\"label\",default:[a])+[ ]+[$e_#eid$],label-sep:0mm)",
                ),

        length_scale: "0.2",
      ),
    )
  } else {
    it
  }

  set page(width: auto, height: auto, margin: (x: 1mm, y: 1mm))

  let title = sys.inputs.at("title", default: "A")
  let data_path = sys.inputs.at("data_path", default: none)

  if data_path == none {
    text(fill: gray)[No data path provided.]
  } else {
    let text = read(data_path)
    raw(text, lang: "dot")
  }
}
