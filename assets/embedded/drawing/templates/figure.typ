// Default per-graph Typst template bundled with the linnet CLI.
// The CLI passes --input title="..." and --input data-path="..." for every build.
#import "layout.typ": layout


#set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm))

#let title = sys.inputs.at("title", default: "A")
#let data-path = sys.inputs.at("data-path", default: none)
#let typst-fields = sys.inputs.at("typst-fields", default: "plain")
#let bool-input(value) = if type(value) == bool {
  value
} else {
  let value = str(value).trim("\"")
  value in ("true", "True", "TRUE", "on", "On", "ON", "yes", "Yes", "YES", "1")
}
#let mark-input(value) = {
  let value = str(value).trim("\"")
  if value == "auto" {
    auto
  } else if value == "none" {
    none
  } else if value.starts-with("(") {
    eval(value)
  } else {
    value
  }
}
#let edge-style-options = {
  let options = eval(sys.inputs.at("edge-style-options", default: "(:)"))
  let momentum-arrows = sys.inputs.at("momentum-arrows", default: none)
  if momentum-arrows != none {
    options = options + (momentum-arrows: bool-input(momentum-arrows))
  }
  let show-edge-index = sys.inputs.at("show-edge-index", default: none)
  if show-edge-index != none {
    options = options + (show-edge-index: bool-input(show-edge-index))
  }
  let show-half-edge-index = sys.inputs.at("show-half-edge-index", default: none)
  if show-half-edge-index != none {
    options = options + (show-half-edge-index: bool-input(show-half-edge-index))
  }
  let show-particle = sys.inputs.at("show-particle", default: none)
  if show-particle != none {
    options = options + (show-particle: bool-input(show-particle))
  }
  let momentum-arrow-mark = sys.inputs.at("momentum-arrow-mark", default: none)
  if momentum-arrow-mark != none {
    options = options + (momentum-arrow-mark: mark-input(momentum-arrow-mark))
  }
  let momentum-arrow-paint = sys.inputs.at("momentum-arrow-paint", default: none)
  let momentum-arrow-thickness = sys.inputs.at("momentum-arrow-thickness", default: none)
  if momentum-arrow-paint != none or momentum-arrow-thickness != none {
    let stroke = (paint: black, thickness: 0.55pt, cap: "round")
    if momentum-arrow-paint != none {
      stroke = stroke + (paint: eval(str(momentum-arrow-paint).trim("\"")))
    }
    if momentum-arrow-thickness != none {
      stroke = stroke + (thickness: eval(str(momentum-arrow-thickness).trim("\"")))
    }
    options = options + (momentum-arrow-stroke: stroke)
  }
  options
}

#show raw: it => if it.at("lang") == "dot" {
  layout(it.at("text"), columns: 1, typst-fields: typst-fields, edge-style-options: edge-style-options)
} else {
  it
}

#if data-path == none {
  text(fill: gray)[No data path provided.]
} else {
  let text = read(data-path)
  raw(text, lang: "dot")
}
