// Default per-graph Typst template bundled with the linnet CLI.
// The CLI passes --input title="..." and --input data="..." for every build.
#{
import "layout.typ": layout
import "edge-style.typ": *
show raw: it => if it.at("lang") == "dot"{
  layout(it.at("text"),columns: 3,scope:(map:map),additional_data:(node:(eval:
          "(stroke:blue,fill :black,
          radius:2pt,
          outset: -1pt)"),edge:(
            eval_source:"(label:[{particle}],..(map.at(\"{particle}\").source))",
            eval_sink:"map.at(\"{particle}\").sink"
          // eval_source:"..(map.at(\"{particle}\").source)", //",label:[{particle}],label-pos:100%,label-side:left)"
          // eval_sink:"..(map.at(\"{particle}\").sink)"
        ),length_scale:0.2))
  }else{
    it
  }

set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm))

let title = sys.inputs.at("title", default: "A")
let data_path = sys.inputs.at("data_path", default: none)

if data_path == none {
    text(fill: gray)[No data path provided.]
  } else {
    let text = read(data_path)
    raw(text, lang: "dot")
  }
}
