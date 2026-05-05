// Default per-graph Typst template bundled with the linnet CLI.
// The CLI passes --input title="..." and --input data-path="..." for every build.
#import "layout.typ": layout


#set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm))

#let title = sys.inputs.at("title", default: "A")
#let data-path = sys.inputs.at("data-path", default: none)
#let typst-fields = sys.inputs.at("typst-fields", default: "plain")

#if data-path == none {
  text(fill: gray)[No data path provided.]
} else {
  let text = read(data-path)
  raw(text, lang: "dot")
}
// #text(fill: gray)[Placeholder figure for #data-path]
// TODO: import and visualize the DOT data here.

#show raw: it => if it.at("lang") == "dot" {
  layout(it.at("text"), columns: 1, typst-fields: typst-fields)
} else {
  it
}


```dot
digraph default {
    ext    [style=invis]
    ext -> v7:1 [particle=a, id=0];
    ext -> v8:2 [particle=a, id=1];
    v9:3 -> ext [particle=a, id=2];
    v10:4 -> ext [particle=a, id=3];
    v11:5 -> ext [particle=a, id=4];
    v12:0 -> ext [particle=a, id=5];
    v7 -> v8 [pdg=6, id=6];
    v8 -> v9 [pdg=6, id=7];
    v9 -> v10 [pdg=6, id=8];
    v10 -> v11 [pdg=6, id=9];
    v11 -> v12 [pdg=6, id=10, lmb_id=0];
    v12 -> v7 [pdg=6, id=11];
}
```
