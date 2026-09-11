// GammaLoop supplies physics styles and external-edge placement. Location-specific
// adapters bind Linnest's shared renderer; it owns measurement, constraints,
// subgraphs, layout passes, generic style composition, and drawing.

#let _field(record, name, default: none) = {
  let data = record.at("data", default: none)
  if type(data) == dictionary and data.keys().contains(name) {
    return data.at(name)
  }
  record
    .at("fields", default: record.at("statements", default: (:)))
    .at(name, default: default)
}

#let _resolved-mode(
  g,
  graph: none,
  amplitude-mode: false,
  cross-section-mode: false,
  auto-mode: false,
) = {
  if not auto-mode {
    return (amplitude: amplitude-mode, cross-section: cross-section-mode)
  }
  let external = graph
    .edges(g)
    .filter(edge => edge.source == none or edge.sink == none)
  let cut = external.any(edge => _field(edge, "is_cut") != none)
  (amplitude: external.len() > 0 and not cut, cross-section: cut)
}

// Match GammaLoop external-edge conventions with outward-facing particle
// labels. Incoming order pairs the two sides of a cross section; the sewing
// tag never becomes a momentum index. Matched legs share fixed Y coordinates
// across the left and right columns. Every dangling endpoint stays at raw
// depth zero, independently of explicit XY placement or an unmatched cut tag.
#let autogen-external-edge-fields(
  g,
  graph: none,
  match-field: none,
  place: true,
  y-scale: 10,
) = {
  let left = ()
  let right = ()
  for edge in graph.edges(g) {
    let id = if match-field == none { edge.edge } else {
      _field(edge, match-field)
    }
    if id != none and edge.source == none and edge.sink != none {
      left.push(id)
    } else if id != none and edge.source != none and edge.sink == none {
      right.push(id)
    }
  }
  graph.map(g, edge: edge => {
    let side = if edge.source == none and edge.sink != none { "left" } else if (
      edge.source != none and edge.sink == none
    ) { "right" } else { return none }
    let generated = (pos: graph.pos(z: graph.pin(0)))
    // Native drawing data retains precedence over mode-generated side labels
    // and XY placements. The external depth pin is deliberately unconditional.
    let data = edge.at("data", default: none)
    if (
      (type(data) != dictionary or not data.keys().contains("label-anchor"))
        and not edge.statements.keys().contains("label-anchor")
    ) {
      // Resolve the anchor after layout: explicit coordinates can put either
      // incoming or outgoing legs on either side of the graph.
      generated.insert("label-anchor", auto)
    }
    let id = if match-field == none { edge.edge } else {
      _field(edge, match-field)
    }
    let ids = if match-field != none or side == "left" { left } else { right }
    let rank = ids.position(value => value == id)
    if rank != none and place {
      let position = (z: graph.pin(0))
      if not edge.at("pos-x-set", default: false) {
        position.x = graph.group(side, side: if side == "left" { "-" } else {
          "+"
        })
        position.dx = if side == "left" { -y-scale } else { y-scale }
      }
      if not edge.at("pos-y-set", default: false) {
        let y = ((ids.len() - 1) / 2 - rank) * y-scale
        position.y = if match-field == none { graph.start(y) } else {
          graph.pin(y)
        }
      }
      generated.pos = graph.pos(..position)
    }
    generated
  })
}

#let layout(
  input,
  graph: none,
  renderer: none,
  physics: none,
  edge-style: none,
  diagram-options: (:),
  scope: (:),
  columns: 1fr,
  unit: 1,
  typst-fields: "plain",
  edge-style-options: (:),
  show-node-index: false,
  debug: false,
  amplitude-mode: false,
  cross-section-mode: false,
  additional-data: (:),
  elements: (:),
  style-options: (:),
  layout-passes: ((:),),
  auto-mode: false,
) = {
  // Explicit template-option map entries override the generated model map,
  // which itself extends the model-neutral aliases.
  let options = (
    (default: edge-style.default-edge)
      + edge-style-options
      + (
        map: physics.default-map
          + edge-style.map
          + edge-style-options.at("map", default: (:)),
        typst-fields: typst-fields,
      )
  )
  let show-momentum = options.at("show-momentum", default: auto)
  if show-momentum == auto {
    show-momentum = options.at("momentum-arrows", default: false)
  }
  let label-length-scale = if show-momentum in (true, "true", "\"true\"") {
    0.45
  } else { 0.30 }
  let styles = (
    (
      scope: scope,
      unit: unit,
      node-label: if debug or show-node-index {
        node => [$n_(#node.vid)$]
      } else {
        auto
      },
      node-label-style: (padding: 0.08),
    )
      + physics.style(..options)
      + style-options
  )
  for g in graph.parse(input) {
    g = renderer.attach-elements(g, elements)
    let mode = _resolved-mode(
      g,
      graph: graph,
      amplitude-mode: amplitude-mode,
      cross-section-mode: cross-section-mode,
      auto-mode: auto-mode,
    )
    if mode.amplitude or mode.cross-section {
      g = autogen-external-edge-fields(
        g,
        graph: graph,
        match-field: if mode.cross-section { "is_cut" } else { none },
      )
    }
    let defaults = if mode.amplitude {
      (
        k-spring: 4.5,
        eps: 1e-7,
        step: 0.6,
        gamma-dangling: 2.3,
        label-length-scale: label-length-scale,
        label-steps: 100,
        directional-force: 4.5,
        label-layout: "dangling-tangent",
      )
    } else if mode.cross-section {
      (
        length-scale: 0.4,
        label-length-scale: label-length-scale,
        label-steps: 100,
        label-layout: "dangling-tangent",
      )
    } else { (:) }
    renderer.layout-graph(
      (
        style: styles,
        draw: (show-half-edge-ids: debug) + diagram-options,
        layouts: layout-passes,
        layout-defaults: defaults + additional-data,
      ),
      g,
    )
  }
}

// V1 template adapter. The caller supplies native Typst dictionaries and
// arrays; this layer never parses strings or evaluates source fragments.
#let render-layout(
  config,
  graph: none,
  renderer: none,
  physics: none,
  edge-style: none,
  diagram-options: (:),
) = {
  let path = config.at("data-path", default: none)
  if path == none { panic("render config requires data-path") }
  let options = config.at("options", default: (:))
  let amplitude = options.at("amplitude-mode", default: false)
  let cross-section = options.at("cross-section-mode", default: false)
  assert(
    not (amplitude and cross-section),
    message: "amplitude-mode and cross-section-mode cannot both be true",
  )
  let mode = options.at("mode", default: if amplitude {
    "amplitude"
  } else if cross-section { "cross-section" } else { "auto" })
  assert(
    mode in ("auto", "generic", "amplitude", "cross-section"),
    message: "config.options.mode must be auto, generic, amplitude, or cross-section",
  )
  let physics-options = options
  for key in (
    "mode",
    "amplitude-mode",
    "cross-section-mode",
    "show-node-index",
    "debug",
    "columns",
    "rows",
  ) {
    if physics-options.keys().contains(key) {
      let _ = physics-options.remove(key)
    }
  }
  if (
    not physics-options.keys().contains("momentum-arrow-stroke")
      and (
        physics-options.keys().contains("momentum-arrow-paint")
          or physics-options.keys().contains("momentum-arrow-thickness")
      )
  ) {
    physics-options.insert("momentum-arrow-stroke", (
      paint: physics-options.at(
        "momentum-arrow-paint",
        default: physics.palette.ink,
      ),
      thickness: physics-options.at("momentum-arrow-thickness", default: 1pt),
      cap: "round",
    ))
  }
  for key in ("momentum-arrow-paint", "momentum-arrow-thickness") {
    if physics-options.keys().contains(key) {
      let _ = physics-options.remove(key)
    }
  }
  let styles = config.at("style", default: (:))
  let drawing = (
    (title: config.at("title", default: auto))
      + diagram-options
      + config.at("draw", default: (:))
  )
  layout(
    read(path),
    graph: graph,
    renderer: renderer,
    physics: physics,
    edge-style: edge-style,
    diagram-options: drawing,
    unit: drawing.at("unit", default: styles.at("unit", default: 1.5)),
    scope: drawing.at("scope", default: styles.at("scope", default: (:))),
    // The V1 boundary accepts native values only; executable field evaluation
    // remains an explicit choice of the standalone physics callback API.
    typst-fields: "plain",
    edge-style-options: physics-options,
    show-node-index: options.at("show-node-index", default: false),
    debug: options.at("debug", default: false),
    amplitude-mode: mode == "amplitude",
    cross-section-mode: mode == "cross-section",
    auto-mode: mode == "auto",
    elements: config.at("elements", default: (:)),
    style-options: styles,
    layout-passes: config.at("layouts", default: ((:),)),
  )
}

// Bind one location's Linnest package and particle-style adapter once. The
// returned function preserves the public GammaLoop layout signature.
#let bind-layout(
  graph: none,
  renderer: none,
  physics: none,
  edge-style: none,
  diagram-options: (:),
) = layout.with(
  graph: graph,
  renderer: renderer,
  physics: physics,
  edge-style: edge-style,
  diagram-options: diagram-options,
)

// Bind the same location-specific dependencies to the V1 render(config)
// contract used by Clinnet-generated entrypoints.
#let bind-render(
  graph: none,
  renderer: none,
  physics: none,
  edge-style: none,
  diagram-options: (:),
) = render-layout.with(
  graph: graph,
  renderer: renderer,
  physics: physics,
  edge-style: edge-style,
  diagram-options: diagram-options,
)
