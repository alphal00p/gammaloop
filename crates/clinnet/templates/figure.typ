// Default per-graph Typst template bundled with the linnet CLI.
// The CLI passes --input title="..." and both spellings of the data-path input.
#import "layout.typ": layout


#set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm))

#let title = sys.inputs.at("title", default: "A")
#let data-path = sys.inputs.at(
  "data-path",
  default: sys.inputs.at("data_path", default: none),
)

#show raw: it => if it.at("lang") == "dot" {
  layout(it.at("text"), columns: 1)
} else {
  it
}

#if data-path == none {
  text(fill: gray)[No data path provided.]
} else {
  let text = read(data-path)
  raw(text, lang: "dot")
}
