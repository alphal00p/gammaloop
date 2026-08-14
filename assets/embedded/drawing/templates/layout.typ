#import "crates/linnest/typst/src/lib.typ": draw, graph, layout as apply-layout, physics
#import "edge-style.typ" as edge-style

#let edge-label-style(edge) = {
  let raw-edge = edge.at("edge", default: (:))
  let statements = raw-edge.at("statements", default: (:))
  let anchor = statements.at("label-anchor", default: edge.at("label-anchor", default: none))
  if type(anchor) == str {
    (anchor: anchor)
  } else {
    (:)
  }
}

#let _node-index-label(node) = {
  let index = node.vid
  [$n_#index$]
}

#let _particle-label(edge, index, typst-fields, include-index: true) = {
  let scope = edge.fields + (edge: edge.edge)
  let entry = physics.edge-entry(scope, map: edge-style.map, default: (label: none))
  let label = physics.label-content(
    entry.at("label", default: none),
    scope,
    mode: typst-fields,
    map: edge-style.map,
    scope: scope,
  )
  if label == none {
    return none
  }
  if include-index {
    label + [$(p_#index)$]
  } else {
    label
  }
}

#let _momentum-edge-label(edge, typst-fields, edge-style-options) = {
  let momentum-arrows = edge-style-options.at("momentum-arrows", default: false)
  let show-momentum-index = momentum-arrows and edge-style-options.at(
    "show-edge-index",
    default: true,
  )
  let label-options = if show-momentum-index {
    edge-style-options + (
      show-edge-index: false,
    )
  } else {
    edge-style-options
  }
  let label = edge-style.edge-label(
    edge,
    typst-fields: typst-fields,
    ..label-options,
  )
  if not show-momentum-index {
    return label
  }
  let index = edge.at("momentum-index", default: physics.edge-index(edge))
  if index == none {
    label
  } else if label == none {
    [$p_#index$]
  } else {
    label + [$(p_#index)$]
  }
}

#let _rank(ids, id) = {
  for (rank, value) in ids.enumerate() {
    if value == id {
      return rank
    }
  }
  none
}

#let _field(record, name, default: none) = {
  record
    .at("fields", default: record.at("statements", default: (:)))
    .at(name, default: default)
}

// Match GammaLoop external-edge conventions with outward-facing particle
// labels. A match field pairs the two sides of a cross section without exposing
// its internal sewing tag as the momentum index.
#let autogen-external-edge-fields(
  g,
  typst-fields: "plain",
  match-field: none,
  place: true,
  y-scale: 10,
  include-index: true,
) = {
  let left = ()
  let right = ()
  for edge in graph.edges(g) {
    let id = if match-field == none {
      edge.edge
    } else {
      _field(edge, match-field)
    }
    if id != none and edge.source == none and edge.sink != none {
      left.push(id)
    } else if id != none and edge.source != none and edge.sink == none {
      right.push(id)
    }
  }

  graph.map(g, edge: edge => {
    let side = if edge.source == none and edge.sink != none {
      "left"
    } else if edge.source != none and edge.sink == none {
      "right"
    } else {
      none
    }
    if side == none {
      return none
    }

    let id = if match-field == none {
      edge.edge
    } else {
      _field(edge, match-field)
    }
    if id == none {
      return none
    }
    // Incoming order defines p_i; the sewing tag only maps its partner on the
    // right back to the same rank.
    let ids = if match-field != none {
      left
    } else if side == "left" {
      left
    } else {
      right
    }
    let rank = _rank(ids, id)
    if rank == none {
      return none
    }
    let index = if match-field == none { edge.edge } else { rank }
    let patch = (
      "label-anchor": if side == "left" { "east" } else { "west" },
      "momentum-index": index,
    )
    if place {
      let x = if side == "left" {
        graph.group("left", side: "-")
      } else {
        graph.group("right", side: "+")
      }
      patch.insert(
        "pos",
        graph.pos(x: x, y: graph.start(((ids.len() - 1) / 2 - rank) * y-scale)),
      )
    }
    let explicit-label = edge.fields.at(
      "display-label",
      default: edge.fields.at("label", default: none),
    )
    if explicit-label == none {
      let label = _particle-label(
        edge,
        index,
        typst-fields,
        include-index: include-index,
      )
      if label != none {
        patch.insert("label", label)
      }
    }
    patch
  })
}

#let layout(
  input,
  split-edge: true,
  scope: (:),
  columns: 1fr,
  unit: 1,
  typst-fields: "plain",
  edge-style-options: (:),
  show-node-index: false,
  amplitude-mode: false,
  cross-section-mode: false,
  additional-data: (:),
) = {
  let graphs = graph.parse(input)
  let diags = ()
  for graph-bytes in graphs {
    let momentum-arrows = edge-style-options.at("momentum-arrows", default: false)
    if amplitude-mode {
      graph-bytes = autogen-external-edge-fields(
        graph-bytes,
        typst-fields: typst-fields,
        include-index: not momentum-arrows,
      )
    } else if cross-section-mode {
      graph-bytes = autogen-external-edge-fields(
        graph-bytes,
        typst-fields: typst-fields,
        match-field: "is_cut",
        place: false,
        include-index: not momentum-arrows,
      )
    }
    let layout-options = if amplitude-mode {
      (
        k-spring: 4.5,
        eps: 1e-7,
        step: 0.6,
        gamma-dangling: 2.3,
        label-length-scale: 1.2,
        label-steps: 100,
        directional-force: 4.5,
        label-layout: "dangling-tangent",
      ) + additional-data
    } else if cross-section-mode {
      (label-layout: "dangling-tangent") + additional-data
    } else {
      additional-data
    }
    graph-bytes = apply-layout(graph-bytes, ..layout-options)
    diags.push(draw(
      graph-bytes,
      scope: scope,
      unit: unit,
      title: auto,
      node-label: if show-node-index { _node-index-label } else { auto },
      source-style: edge => edge-style.source-style(edge, typst-fields: typst-fields, ..edge-style-options),
      sink-style: edge => edge-style.sink-style(edge, typst-fields: typst-fields, ..edge-style-options),
      edge-label: edge => _momentum-edge-label(
        edge,
        typst-fields,
        edge-style-options,
      ),
      edge-label-style: edge-label-style,
    ))
  }
  for d in diags {
    d
  }
}
