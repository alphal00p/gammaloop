#import "crates/linnest/typst/src/lib.typ": draw, graph, layout
#import graph: edge, node, pin, pos, sink, source
#import "map-style.typ" as feynman

// Keep the pre-configuration presets explicit so default comparisons do not
// merely compare two callbacks backed by the same new implementation.
#let legacy-stroke = (paint: black, thickness: 0.5pt, cap: "round")
#let legacy-node = (radius: 0.18, fill: white, stroke: legacy-stroke)
#let legacy-fermion = (
  stroke: legacy-stroke,
  mark: (
    end: (
      symbol: ">",
      fill: black,
      stroke: black + 0.2pt,
      anchor: "center",
      shorten-to: auto,
    ),
    scale: .5,
  ),
  mark-position: "center-if-dangling",
  mark-orientation: "edge",
)
#let legacy-photon = (
  stroke: legacy-stroke,
  pattern: "wave",
  pattern-amplitude: 0.10,
  pattern-wavelength: 0.50,
)
#let legacy-gluon = (
  stroke: legacy-stroke,
  pattern: "coil",
  pattern-amplitude: 0.15,
  pattern-wavelength: 0.60,
  pattern-coil-longitudinal-scale: 1.60,
)
#let legacy-scalar = (
  stroke: legacy-stroke + (thickness: 1pt, dash: (0.1em, 0.45em)),
)
#let legacy-particles = (
  a: legacy-photon,
  photon: legacy-photon,
  g: legacy-gluon,
  gluon: legacy-gluon,
  scalar: legacy-scalar,
  ghG: legacy-scalar,
)
#let legacy-geometry = (
  offset: 0.62,
  offset-side: "label",
  label-side: auto,
  label-gap: 0.45,
  label-style: (anchor: auto),
)
#let legacy-arrow = (
  legacy-geometry
    + (
      length: 1.0,
      shift: 0,
      ratio: none,
      resolve-length: "length",
      stroke: (paint: black, thickness: 0.4pt, cap: "round"),
      mark: (
        end: (
          symbol: "straight",
          fill: black,
          stroke: black + 0.4pt,
          anchor: "center",
          shorten-to: auto,
        ),
        scale: 0.50,
      ),
      mark-position: "end",
      mark-orientation: "path",
    )
)
#let legacy-label = (
  legacy-geometry
    + (
      length: none,
      ratio: none,
      resolve-length: "none",
      shift: 0,
      stroke: none,
      mark: none,
      label: [],
      label-shift: 0,
    )
)
#let legacy-graph = (
  unit: 1.35,
  node-label: none,
  node-style: node => if node.fields.at("hidden", default: false) {
    (radius: 0, fill: none, stroke: none)
  } else { legacy-node },
  edge-label: edge => hide([$p_(#edge.eid)$]),
  edge-label-style: (anchor: "center", padding: 0.05),
)

#assert.eq(type(feynman.graph-style), function)
#assert.eq(type(feynman.draw-style), dictionary)
#assert.eq(feynman.draw-style.padding, 1.5)
#assert.eq(feynman.draw-style.edge-dangling-tangent, "horizontal")
#let defaults = feynman.graph-style()
#let empty-edge = (fields: (:), momentum: [])
#let legacy-layers = (
  legacy-fermion + (mark-shift: 0),
  legacy-arrow,
  legacy-label,
)
#assert.eq(defaults.unit, 1.35)
#assert.eq(defaults.node-label, none)
#assert.eq(defaults.edge-label-style, legacy-graph.edge-label-style)
#assert.eq(defaults.scope.feynman, (
  node-style: legacy-node,
  fermion: legacy-fermion,
  particles: legacy-particles,
  momentum-stroke: legacy-arrow.stroke,
  momentum-mark: legacy-arrow.mark,
))
#assert.eq(feynman.node-style((fields: (:))), legacy-node)
#assert.eq(feynman.edge-style(empty-edge), legacy-layers)
#assert.eq(feynman.node-style(defaults.scope + (fields: (:))), legacy-node)
#assert.eq(
  (feynman.draw-style.edge-style)(defaults.scope + empty-edge),
  legacy-layers,
)

// Exercise each width independently as well as all overrides together.
#let wider = feynman.graph-style(line-width: 1.25pt)
#let domain = wider.scope.feynman
#assert.eq(domain.node-style.stroke.thickness, 1.25pt)
#assert.eq(domain.fermion.stroke.thickness, 1.25pt)
#assert.eq(domain.particles.photon.stroke.thickness, 1.25pt)
#assert.eq(domain.particles.gluon.stroke.thickness, 1.25pt)
#assert.eq(domain.particles.scalar.stroke.thickness, 2.5pt)
#assert.eq(domain.fermion.mark.end.stroke.thickness, 0.5pt)
#assert.eq(domain.momentum-stroke.thickness, 1pt)
#assert.eq(domain.momentum-mark.end.stroke.thickness, 1pt)
#for (option, width) in (
  ("node-line-width", 1.75pt),
  ("massive-line-width", 3pt),
  ("fermion-arrow-line-width", 0.625pt),
  ("momentum-line-width", 0.875pt),
) {
  let expected = domain
  if option == "node-line-width" {
    expected.node-style.stroke.thickness = width
  } else if option == "massive-line-width" {
    expected.particles.scalar.stroke.thickness = width
    expected.particles.ghG.stroke.thickness = width
  } else if option == "fermion-arrow-line-width" {
    expected.fermion.mark.end.stroke = black + width
  } else {
    expected.momentum-stroke.thickness = width
    expected.momentum-mark.end.stroke = black + width
  }
  let patch = (:)
  patch.insert(option, width)
  let configured = feynman.graph-style(line-width: 1.25pt, ..patch)
  assert.eq(configured.scope.feynman, expected, message: option)
  patch.insert(option, auto)
  assert.eq(
    feynman.graph-style(line-width: 1.25pt, ..patch).scope.feynman,
    domain,
    message: option + ": explicit auto must derive from line-width",
  )
}
#let overridden = feynman.graph-style(
  unit: 10pt,
  node-radius: 0.3,
  line-width: 1.25pt,
  node-line-width: 1.75pt,
  massive-line-width: 3pt,
  fermion-arrow-line-width: 0.625pt,
  momentum-line-width: 0.875pt,
)
#let variants = (
  defaults: defaults,
  legacy: legacy-graph,
  scaled: feynman.graph-style(unit: 20pt, node-radius: 0.3),
  derived: feynman.graph-style(
    unit: 10pt,
    node-radius: 0.3,
    line-width: 1.25pt,
  ),
  overrides: overridden,
)
#assert.eq(variants.scaled.unit, 20pt)
#assert.eq(
  variants.scaled.scope.feynman,
  defaults.scope.feynman
    + (
      node-style: legacy-node + (radius: 0.3),
    ),
)

// Aliases, hidden nodes, and sparse edge settings must work in every scope,
// including direct calls with no scope and calls after unrelated configurations.
#for scope in (
  (:),
  defaults.scope,
  wider.scope,
  overridden.scope,
  defaults.scope,
) {
  let expected = scope.at("feynman", default: defaults.scope.feynman)
  let particle-only = feynman.edge-style.with(show-momentum: false)
  assert.eq(particle-only(scope + (fields: (:), momentum: [$k$])), (
    expected.fermion + (mark-shift: 0),
  ))
  for (alias, preset) in expected.particles {
    let layers = feynman.edge-style(
      scope + (fields: (particle: alias), momentum: []),
    )
    assert.eq(layers.first(), preset + (mark-shift: 0), message: alias)
    assert.eq(
      feynman.edge-style(
        scope + (fields: (particle: alias), momentum: [$k$]),
        show-momentum: false,
      ),
      (layers.first(),),
    )
  }
  assert.eq(expected.particles.a, expected.particles.photon)
  assert.eq(expected.particles.g, expected.particles.gluon)
  assert.eq(expected.particles.ghG, expected.particles.scalar)
  for particle in ("d", "unknown", "\"unknown\"") {
    assert.eq(
      feynman
        .edge-style(scope + (fields: (particle: particle), momentum: []))
        .first(),
      expected.fermion + (mark-shift: 0),
    )
  }
  for hidden in (true, "true", "\"true\"") {
    assert.eq(feynman.node-style(scope + (fields: (hidden: hidden))), (
      radius: 0,
      fill: none,
      stroke: none,
    ))
  }
  for hidden in (false, "false", "\"false\"") {
    assert.eq(
      feynman.node-style(scope + (fields: (hidden: hidden))),
      expected.node-style,
    )
  }
  let layers = feynman.edge-style(
    scope
      + (
        momentum: [$k$],
        fields: (
          particle: "\"photon\"",
          route: "\"straight-through\"",
          cut: "true",
          cut-gap: "0.25",
          crossing-under: <over>,
          crossing-gap: "0.75",
          fermion-arrow-shift: "-0.5",
        )
          + feynman.momentum(
            side: "right",
            offset: -0.4,
            length: 0.7,
            shift: -0.4,
            label: (gap: 0.2, shift: -0.75, anchor: "south-west"),
          ),
      ),
  )
  assert.eq(
    layers.first(),
    expected.particles.photon
      + (
        route: "straight-through",
        split-gap: 0.25,
        crossing-under: <over>,
        crossing-gap: 0.75,
        mark-shift: -0.5,
      ),
  )
  let geometry = (
    route: "straight-through",
    offset: -0.4,
    offset-side: none,
    label-side: "right",
    label-gap: 0.2,
    label-style: (anchor: "south-west"),
  )
  assert.eq(
    layers.at(1),
    legacy-arrow
      + geometry
      + (
        length: 0.7,
        shift: -0.4,
        stroke: expected.momentum-stroke,
        mark: expected.momentum-mark,
      ),
  )
  assert.eq(
    layers.last(),
    legacy-label + geometry + (label: [$k$], label-shift: -0.75),
  )
  let shifted = feynman.edge-style(
    scope
      + (
        momentum: [],
        fields: feynman.momentum(shift: 0.5),
      ),
  )
  assert.eq(shifted.at(1).shift, 0.5)
  assert.eq(shifted.last().label-shift, 0.5)
}
#assert.eq(feynman.graph-style().scope, defaults.scope)
#assert.eq(feynman.node-style((fields: (:))), legacy-node)
#assert.eq(feynman.edge-style(empty-edge), legacy-layers)

#let base = graph.build(
  {
    node(<a>, pos: pos(x: pin(0), y: pin(0)))
    node(<b>, pos: pos(x: pin(6), y: pin(0)))
    node(<hidden>, hidden: true, pos: pos(x: pin(0), y: pin(3)))
    edge(source(<a>), sink(<b>), pos: pos(x: pin(3), y: pin(0)), momentum: [])
    edge(
      source(<hidden>),
      pos: pos(x: pin(6), y: pin(3)),
      particle: "scalar",
      momentum: [],
    )
    edge(
      source(<a>),
      sink(<b>),
      pos: pos(x: pin(3), y: pin(-2)),
      particle: "a",
      momentum: [],
    )
    edge(
      source(<a>),
      sink(<b>),
      pos: pos(x: pin(3), y: pin(-4)),
      particle: "g",
      momentum: [],
    )
  },
  default-edge-data: feynman.momentum(offset: 0.4, label: (shift: 1)),
)
#let base = graph.map(base, edge: edge => feynman.momentum(label: (gap: 0.2)))

// Each entry compiles independently, without importing the inherited named-map
// fixture. No draw-time scope or unit is supplied: graph.style must carry both.
#let cases = (:)
#for (name, config) in variants {
  cases.insert(name, context {
    let styled = graph.style(base, ..config)
    let expected-node = if name == "legacy" { legacy-node } else {
      config.scope.feynman.node-style
    }
    let unit = if type(config.unit) in (int, float) {
      config.unit * 1em
    } else { config.unit }
    let unit = unit.to-absolute()
    for node in graph.nodes(styled) {
      let diameter = if node.name == <hidden> { 0 } else {
        2 * expected-node.radius
      }
      for dimension in ("layout-width", "layout-height") {
        assert(
          calc.abs(float(node.statements.at(dimension)) - diameter) < 1e-6,
          message: name + ": measured node radius",
        )
      }
    }
    let label-size = measure(text(
      top-edge: "cap-height",
      bottom-edge: "bounds",
      (config.edge-label)((eid: 0)),
    ))
    assert(
      calc.abs(
        float(graph.edges(styled).first().statements.at("label-width"))
          - (label-size.width / unit + 0.1),
      )
        < 1e-6,
      message: name + ": label measurement must use unit",
    )
    let solved = layout(
      styled,
      labels: (steps: 0),
      solver: (steps: 2, seed: 17),
    )
    draw(solved, ..feynman.draw-style, edge-style: edge => {
      assert.eq(edge.fields.momentum-arrow-offset, 0.4)
      assert.eq(edge.fields.momentum-label-shift, 1)
      assert.eq(edge.fields.momentum-label-gap, 0.2)
      if name == "legacy" {
        (
          legacy-particles.at(
            edge.fields.at("particle", default: "d"),
            default: legacy-fermion,
          )
            + (mark-shift: 0),
          legacy-arrow + (offset: 0.4, label-gap: 0.2),
          legacy-label + (offset: 0.4, label-gap: 0.2, label-shift: 1),
        )
      } else {
        assert.eq(
          edge.feynman,
          config.scope.feynman,
          message: name + ": inherited drawing scope",
        )
        let layers = (feynman.draw-style.edge-style)(edge)
        assert.eq(layers.at(1).offset, 0.4)
        assert.eq(layers.last().label-gap, 0.2)
        assert.eq(layers.last().label-shift, 1)
        layers
      }
    })
  })
}
