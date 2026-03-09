// Default per-graph Typst template bundled with the linnet CLI.
// The CLI passes --input title="..." and --input data="..." for every build.
#{
import "layout.typ": layout
show raw: it => if it.at("lang") == "dot"{
    layout(it.at("text"),columns: 1)
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
  // #text(fill: gray)[Placeholder figure for #data_path]
  // TODO: import and visualize the DOT data here.
}
