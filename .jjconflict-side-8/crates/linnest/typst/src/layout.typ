#import "graph.typ" as graph-module
#import "subgraph.typ" as subgraph-module

#let _plugin = plugin("../linnest.wasm")

#let _option-rules = (
  spring: (strength: "number", length: "number"),
  repulsion: (
    strength: "number",
    centering: "number",
    "edge-node": "number",
    "edge-edge": "number",
    dangling: "number",
    "dangling-centroid": "number",
    softening: "number",
  ),
  constraints: (
    "side-strength": "number",
    "node-movement": ("fixed", "layout"),
    direction: ("down", "right", "left-to-right", "left-right", "lr"),
    "rank-alignment": (
      "center",
      "start",
      "end",
      "top",
      "north",
      "left",
      "west",
      "bottom",
      "south",
      "right",
      "east",
    ),
    roots: "indices",
    "same-rank": "rank-groups",
  ),
  labels: (
    distance: "number",
    spring: "number",
    repulsion: "number",
    steps: "integer",
    model: ("normal", "dangling-tangent", "fixed-length"),
    step: "number",
    tolerance: "number",
    "max-movement": "number",
  ),
  solver: (
    algorithm: ("force", "anneal", "tree", "dot", "stable-layered", "railroad"),
    steps: "integer",
    epochs: "integer",
    seed: "integer",
    step: "number",
    "step-shrink": "number",
    cooling: "number",
    "acceptance-floor": "number",
    tolerance: "number",
    temperature: "number",
    "max-movement": "number",
    "incremental-energy": "boolean",
    "crossing-penalty": "number",
    "depth-scale": "non-negative-finite",
    "flattening-end": "unit-interval",
  ),
)

#let _option-error(group, field, expected, value) = panic(
  "layout "
    + group
    + "."
    + field
    + " must be "
    + expected
    + ", got "
    + repr(value),
)

#let _check-option-value(group, field, value, rule) = {
  if (
    rule in ("number", "non-negative-finite", "unit-interval")
      and type(value) not in (int, float)
  ) {
    _option-error(group, field, "a number", value)
  } else if (
    rule in ("non-negative-finite", "unit-interval")
      and (
        value != value
          or value < 0
          or value >= calc.inf
          or (rule == "unit-interval" and value > 1)
      )
  ) {
    _option-error(
      group,
      field,
      if rule == "unit-interval" { "a finite fraction between 0 and 1" } else {
        "a non-negative finite number"
      },
      value,
    )
  } else if rule == "integer" and (type(value) != int or value < 0) {
    _option-error(group, field, "a non-negative integer", value)
  } else if rule == "boolean" and type(value) != bool {
    _option-error(group, field, "a boolean", value)
  } else if (
    rule == "indices"
      and (
        type(value) != array
          or not value.all(index => type(index) == int and index >= 0)
      )
  ) {
    _option-error(group, field, "an array of non-negative node indices", value)
  } else if (
    rule == "rank-groups"
      and (
        type(value) != array
          or not value.all(group => (
            type(group) == function
              or (
                type(group) == dictionary
                  and group.at("linnest-kind", default: none)
                    == "linnest-subgraph"
              )
              or (
                type(group) == array
                  and group.all(index => type(index) == int and index >= 0)
              )
          ))
      )
  ) {
    _option-error(
      group,
      field,
      "an array of subgraphs, node-index arrays, or module functions",
      value,
    )
  } else if type(rule) == array and (type(value) != str or value not in rule) {
    _option-error(group, field, "one of " + rule.map(repr).join(", "), value)
  }
}

#let _checked-group(name, value) = {
  if name not in _option-rules {
    panic("layout options do not have a " + repr(name) + " group")
  }
  if value == none {
    return (:)
  }
  if type(value) != dictionary {
    panic("layout " + name + " options must be a dictionary")
  }
  let rules = _option-rules.at(name)
  for key in value.keys() {
    if key not in rules {
      panic(
        "layout "
          + name
          + " options do not have a "
          + repr(key)
          + " field; expected one of "
          + rules.keys().map(repr).join(", "),
      )
    }
    let _ = _check-option-value(name, key, value.at(key), rules.at(key))
  }
  value
}

#let _rank-subgraph(graph, value) = {
  let value = if type(value) == function { value(graph) } else { value }
  if type(value) == dictionary {
    return subgraph-module._impl.validate(graph, value)
  }
  if type(value) != array or not value.all(item => type(item) == int) {
    panic(
      "layout rank-same entries must be subgraph objects, node-index arrays, or module functions",
    )
  }
  subgraph-module.select(graph, nodes: value)
}

/// Construct reusable semantic layout options or sparsely update existing ones.
///
/// Each option group is a plain dictionary. Supplying `base` preserves every
/// field not mentioned by a patch, including fields inside a patched group.
/// Pass the result to @layout with argument spreading.
///
/// ```example
/// #let common = options(
///   spring: (strength: 8, length: 0.5),
///   repulsion: (edge-node: 0.05, dangling: 2),
///   solver: (algorithm: "force", steps: 50),
/// )
/// #let spacious = options(
///   base: common,
///   spring: (length: 0.7),
///   repulsion: (dangling-centroid: 3),
/// )
/// #let g = layout(g, ..spacious)
/// ```
/// -> dictionary
#let options(
  /// Semantic option-group patches. -> arguments
  ..patches,
  /// Existing option dictionary to update. -> dictionary
  base: (:),
) = {
  if type(base) != dictionary {
    panic("layout options base must be a dictionary")
  }
  let result = base
  if patches.pos().len() > 0 {
    panic("layout options only accepts named group patches")
  }
  for (name, patch) in patches.named() {
    let patch = _checked-group(name, patch)
    let previous = result.at(name, default: (:))
    if type(previous) != dictionary {
      panic("layout options base " + name + " field must be a dictionary")
    }
    result.insert(name, previous + patch)
  }
  result
}

/// Apply the linnest layout pass to a graph object.
///
/// This is intentionally a second step: construct or parse a graph first, then
/// call `layout`. Set `layout-algo` to `"force"` for deterministic force
/// integration, `"anneal"` for simulated annealing, `"tree"` for a traversal
/// tree placement, `"dot"` for a Graphviz-like layered placement, or
/// `"stable-layered"` for a stable railroad-inspired layered placement.
///
/// Force layout uses one `steps` × `epochs` iteration budget and one cooling
/// schedule. Effective auxiliary depth is `scale * raw-z`: `scale` follows a
/// smoothstep from `depth-scale` to zero over the initial `flattening-end`
/// fraction of that budget. Remaining iterations are planar, without restarting
/// `step` or cooling. `early-tol` can stop the run only after scale reaches zero.
/// Unpinned raw depths remain bounded; hard raw-depth pins are preserved even
/// when their effective depth is zero. Output positions remain 2D, and auxiliary
/// depth does not change draw order; label layout runs after the single graph pass.
///
/// ```example
/// #let g = graph.parse("digraph partial { a -> b;a -> b;a:s -> b:s; b:s -> c:s; c:s -> d:s; d:s -> a:s }").at(0)
/// #let south = subgraph.compass(g,"s")
/// #let gf = layout(layout(g, layout-algo: "force"), layout-algo: "tree", layout-nodes: "fixed",subgraph:south)
/// #let ga = layout(g, layout-algo: "force")
/// #let gt = layout(layout(g, layout-algo: "tree",layout-roots: (2)), layout-algo: "anneal", layout-nodes: "fixed",gamma-ee:0.1,gamma-ev:.75,beta:5,length-scale:0.4)
/// #grid(columns: 3, gutter: 2cm, draw(gf), draw(ga), draw(gt))
///
/// ```
///
/// ```example
/// #let g = graph.parse("digraph partial { a -> b; b -> c; c -> d; d -> a }").at(0)
/// #let tree = graph.forests(g).at(0)
/// #let g = layout(g, layout-algo: "tree", subgraph: tree)
/// #graph.edges(g).map(edge => edge.pos)
/// ```
/// -> dictionary
#let layout(
  /// Graph object returned by `graph.build` or `graph.parse`.
  /// -> dictionary
  graph,
  /// Semantic incidence-spring options. Supported fields are `strength` and
  /// `length`. These override the corresponding flat parameters below.
  /// -> none | dictionary
  spring: none,
  /// Semantic repulsion options. Supported fields are `strength`, `centering`,
  /// `edge-node`, `edge-edge`, `dangling`, `dangling-centroid`, and
  /// `softening`. These override the corresponding flat parameters below.
  /// -> none | dictionary
  repulsion: none,
  /// Semantic constraint options. Supported fields are `side-strength`,
  /// `node-movement`, `direction`, `rank-alignment`, `roots`, and `same-rank`.
  /// These override the corresponding flat parameters below.
  /// -> none | dictionary
  constraints: none,
  /// Semantic label-layout options. Supported fields are `distance`, `spring`,
  /// `repulsion`, `steps`, `model`, `step`, `tolerance`, and `max-movement`.
  /// These override the corresponding flat parameters below.
  /// -> none | dictionary
  labels: none,
  /// Semantic solver options. Supported fields are `algorithm`, `steps`,
  /// `epochs`, `seed`, `step`, `step-shrink`, `cooling`, `acceptance-floor`,
  /// `tolerance`, `temperature`, `max-movement`, `incremental-energy`,
  /// `crossing-penalty`, `depth-scale`, and `flattening-end`. These override the
  /// corresponding flat parameters below.
  /// -> none | dictionary
  solver: none,
  /// Optional subgraph object to lay out. With `"tree"`, other edges are drawn
  /// from the resulting node positions. With `"dot"` and `"stable-layered"`,
  /// the subgraph determines rank constraints, while all paired edges between
  /// included nodes get dummy routing vertices and edge positions. With
  /// `"force"` and `"anneal"`, nodes and edges outside the subgraph are fixed
  /// boundary points during optimization.
  /// The selection must have compatible topology. -> none | dictionary
  subgraph: none,
  /// Width of the layout viewport used to derive the natural spring length.
  /// Applies to both `"force"` and `"anneal"`. -> float
  viewport-w: 10.0,
  /// Height of the layout viewport used to derive the natural spring length.
  /// Applies to both `"force"` and `"anneal"`. -> float
  viewport-h: 10.0,
  /// Horizontal spacing multiplier for traversal-tree and layered placement.
  /// For `"force"` and `"anneal"`, this scales the initial placement. -> float
  tree-dx: 0.9,
  /// Vertical spacing multiplier for traversal-tree and layered placement.
  /// For `"force"` and `"anneal"`, this scales the initial placement. -> float
  tree-dy: 1.2,
  /// Iterations per epoch. In `"force"` mode this is the number of force
  /// integration steps within the single depth-flattening schedule;
  /// in `"anneal"` mode this is the number of proposals per temperature epoch.
  /// -> int
  steps: 30,
  /// Seed for deterministic initialization, force-mode jitter, and annealing
  /// proposals. Applies to both modes. -> int
  seed: 2,
  /// Initial movement scale. `"force"` multiplies computed forces by this
  /// value, with no restart when depth reaches zero; `"anneal"` uses it as a
  /// proposal step size in natural spring-length units. -> float
  step: 0.81,
  /// Anneal-only step shrink factor, applied when an epoch's acceptance ratio
  /// falls below `accept-floor`. -> float
  step-shrink: 0.21,
  /// Cooling factor applied once per epoch. `"anneal"` cools `temp`; `"force"`
  /// shrinks `step`. -> float
  cool: 0.85,
  /// Anneal-only acceptance-ratio threshold below which `step` is shrunk by
  /// `step-shrink`. -> float
  accept-floor: 0.15,
  /// Force-only early stop threshold for maximum movement in one step, as a
  /// multiple of the natural spring length, eligible only after effective depth
  /// scale reaches zero. The annealing schedule stores this value but does not
  /// currently use it for stopping. -> float
  early-tol: 1e-6,
  /// Anneal-only initial temperature used in the Metropolis acceptance test,
  /// scaled by natural spring length squared. -> float
  temp: 0.3,
  /// Force-mode maximum movement clamp per point and per step, as a multiple
  /// of the natural spring length. -> float
  delta: 0.4,
  /// Base repulsion strength for vertex-vertex interactions. Also scales
  /// `gamma-ev`, `gamma-ee`, `gamma-dangling`,
  /// `gamma-dangling-centroid`, and `g-center`. Applies to both modes through
  /// the shared spring energy.
  /// -> float
  beta: 50.0,
  /// Spring stiffness for node-to-edge incidence lengths. Applies to both
  /// modes. -> float
  k-spring: 11.0,
  /// Centering strength relative to `beta`. Applies to both modes.
  /// -> float
  g-center: 0.002,
  /// Number of epochs. Both modes run up to `steps` iterations inside
  /// each epoch. Force mode shares this budget and its cooling schedule across
  /// depth flattening and the remaining planar relaxation. -> int
  epochs: 30,
  /// Anneal-only fixed energy penalty per detected edge crossing. The direct
  /// force integrator does not currently add a crossing force. -> float
  crossing-penalty: 30.0,
  /// Repulsion for dangling half edges, relative to `beta`. Applies to
  /// both modes through the shared spring energy. -> float
  gamma-dangling: 5.0,
  /// Repulsion of every dangling endpoint from the current node centroid,
  /// relative to `beta`. The equal-and-opposite reaction is shared over the
  /// nodes, avoiding translational drift. Applies to both modes. -> float
  gamma-dangling-centroid: 0.0,
  /// Local edge-edge repulsion, relative to `beta`. Applies to both
  /// modes. -> float
  gamma-ee: 0.1,
  /// Bias that pushes points in directions implied by pin/port constraints.
  /// Force mode treats this as a force and scales it by natural spring length;
  /// anneal mode treats it as a dimensionless multiplier on proposal steps.
  /// Applies to both modes. -> float
  directional-force: 5.0,
  /// Edge-label target offset as a multiple of the graph spring length. Label
  /// layout runs after both graph layout modes. -> float
  label-length-scale: 0.6,
  /// Spring strength pulling each label toward its target offset in the
  /// spring-based label layouts. -> float
  label-spring: 23.0,
  /// Repulsion strength between labels and graph points, scaled by spring
  /// length squared. Label layout runs after both modes. -> float
  label-charge: 3.0,
  /// Maximum number of post-layout label relaxation steps. Set to `0` to
  /// skip label placement. Applies after both modes. -> int
  label-steps: 20,
  /// Edge-label relaxation model. `"normal"` uses a perpendicular offset,
  /// `"dangling-tangent"` uses the edge direction for dangling half-edge labels
  /// and a perpendicular offset for paired edges, and `"fixed-length"` keeps
  /// each label at a fixed distance from its edge point and only lets that
  /// segment rotate. -> string
  label-layout: "normal",
  /// Label relaxation step size. Applies after both modes. -> float
  label-step: 0.15,
  /// Label relaxation early stop threshold. Applies after both modes. -> float
  label-early-tol: 1e-3,
  /// Label movement clamp as a multiple of spring length. Applies after both
  /// modes. -> float
  label-max-delta-scale: 0.5,
  /// Edge-vertex repulsion, relative to `beta`. Applies to both modes.
  /// -> float
  gamma-ev: 0.01,
  /// Softening epsilon used in inverse-square force/energy terms. Applies to
  /// both modes. -> float
  eps: 1e-4,
  /// Whether annealing updates cached energy by local deltas. This is
  /// anneal-only; force mode computes direct forces instead of energies.
  /// -> bool
  incremental-energy: true,
  /// Layout algorithm. Use `"force"` for direct force integration, `"anneal"`
  /// for simulated annealing against the spring energy, `"tree"` for a
  /// traversal-tree placement, `"dot"` for a Graphviz-like layered placement,
  /// or `"stable-layered"` for a stable railroad-inspired layered placement.
  /// -> string
  layout-algo: "force",
  /// Node movement policy. `"layout"` lets the layout algorithm move nodes.
  /// `"fixed"` keeps every node at its current position for this layout pass
  /// and only moves edge control points. With `subgraph`, only edges in the
  /// subgraph are moved; other edge control points stay at their current
  /// positions. -> string
  layout-nodes: "layout",
  /// Direction for traversal-tree and layered rank placement. `"down"` places
  /// increasing ranks downward; `"right"` swaps the layout axes so increasing
  /// ranks go left-to-right and measured node/label widths reserve rank-axis
  /// space. -> string
  layout-direction: "down",
  /// Alignment of real nodes inside a rank along the rank axis. `"center"` keeps
  /// node centers aligned, while `"start"` / `"left"` and `"end"` / `"right"`
  /// align the corresponding measured node-box side. -> string
  rank-align: "center",
  /// Ordered node indices used as preferred roots for `"tree"`, `"dot"`, and
  /// `"stable-layered"`. Roots outside the selected node set are ignored.
  /// Remaining components are laid out afterward in graph order. -> array
  layout-roots: (),
  /// Subgraph objects whose incident nodes should share a dot/stable-layered
  /// rank; node-index arrays and callbacks receiving the graph are also accepted.
  /// These are layout hints supplied by Typst rather than parsed graph
  /// structure. -> array
  rank-same: (),
  /// Dot/stable-layered also honors a node statement `layout-rank` as an exact
  /// non-negative integer rank. Nodes with the same `layout-rank` are placed on
  /// the same horizontal layer, and larger ranks are placed lower.
  /// Relative layout weight for paired edges outside the dot/stable-layered
  /// rank subgraph. Lower values make these edges guide routing without
  /// dominating the rank tree. -> float
  route-edge-weight: 0.15,
  /// Extra horizontal straightening weight for the first or last segment of an
  /// edge with `source-route-exit` or `sink-route-exit` set to a vertical side.
  /// -> float
  route-exit-weight: 4.0,
  /// Multiplier for measured edge-label width when sizing non-rank dummy
  /// routing vertices in dot/stable-layered layout. -> float
  route-label-width-scale: 1.0,
  /// Maximum non-rank dummy label width as a multiple of `tree-dx`. Set to
  /// `0` or a negative value to disable the cap. -> float
  route-label-width-cap: 2.0,
  /// Force-only initial multiplier of raw auxiliary depth, finite and >= 0.
  /// Effective z is this smoothly decreasing scale times raw z, not a spring
  /// pulling raw coordinates toward the plane. A value of 0 starts in 2D;
  /// flattening preserves hard raw-depth pins while removing their invisible
  /// separation from the visible 2D force calculation. -> float
  depth-scale: 1.0,
  /// Force-only fraction of the total `steps` × `epochs` iteration budget over
  /// which the depth scale reaches zero by smoothstep, finite and in 0..=1.
  /// Zero starts in 2D; 1 still evaluates the final iteration on the exact
  /// plane. Neither value adds a phase or restarts the cooling schedule. -> float
  flattening-end: 0.5,
  /// Natural spring-length multiplier. This scales the graph's preferred edge
  /// length and dimensional force/energy terms so changing only this value
  /// mostly zooms the result instead of retuning the force ratios. Applies to
  /// both modes. -> float
  length-scale: 0.35,
) = {
  spring = _checked-group("spring", spring)
  repulsion = _checked-group("repulsion", repulsion)
  constraints = _checked-group("constraints", constraints)
  labels = _checked-group("labels", labels)
  solver = _checked-group(
    "solver",
    (
      depth-scale: depth-scale,
      flattening-end: flattening-end,
    )
      + _checked-group("solver", solver),
  )
  let rank-groups = constraints.at("same-rank", default: rank-same)
  if type(rank-groups) != array {
    panic("layout rank-same/constraints.same-rank must be an array")
  }
  let settings = (
    viewport-w: str(viewport-w),
    viewport-h: str(viewport-h),
    tree-dx: str(tree-dx),
    tree-dy: str(tree-dy),
    steps: str(solver.at("steps", default: steps)),
    seed: str(solver.at("seed", default: seed)),
    step: str(solver.at("step", default: step)),
    step-shrink: str(solver.at("step-shrink", default: step-shrink)),
    cool: str(solver.at("cooling", default: cool)),
    accept-floor: str(solver.at("acceptance-floor", default: accept-floor)),
    early-tol: str(solver.at("tolerance", default: early-tol)),
    temp: str(solver.at("temperature", default: temp)),
    delta: str(solver.at("max-movement", default: delta)),
    beta: str(repulsion.at("strength", default: beta)),
    k-spring: str(spring.at("strength", default: k-spring)),
    g-center: str(repulsion.at("centering", default: g-center)),
    epochs: str(solver.at("epochs", default: epochs)),
    crossing-penalty: str(solver.at(
      "crossing-penalty",
      default: crossing-penalty,
    )),
    gamma-dangling: str(repulsion.at("dangling", default: gamma-dangling)),
    gamma-dangling-centroid: str(
      repulsion.at("dangling-centroid", default: gamma-dangling-centroid),
    ),
    gamma-ee: str(repulsion.at("edge-edge", default: gamma-ee)),
    directional-force: str(constraints.at(
      "side-strength",
      default: directional-force,
    )),
    label-length-scale: str(labels.at("distance", default: label-length-scale)),
    label-spring: str(labels.at("spring", default: label-spring)),
    label-charge: str(labels.at("repulsion", default: label-charge)),
    label-steps: str(labels.at("steps", default: label-steps)),
    label-layout: labels.at("model", default: label-layout),
    label-step: str(labels.at("step", default: label-step)),
    label-early-tol: str(labels.at("tolerance", default: label-early-tol)),
    label-max-delta-scale: str(
      labels.at("max-movement", default: label-max-delta-scale),
    ),
    gamma-ev: str(repulsion.at("edge-node", default: gamma-ev)),
    eps: str(repulsion.at("softening", default: eps)),
    incremental-energy: solver.at(
      "incremental-energy",
      default: incremental-energy,
    ),
    layout-algo: solver.at("algorithm", default: layout-algo),
    layout-nodes: constraints.at("node-movement", default: layout-nodes),
    layout-direction: constraints.at("direction", default: layout-direction),
    rank-align: constraints.at("rank-alignment", default: rank-align),
    layout-roots: constraints.at("roots", default: layout-roots),
    rank-same: rank-groups.map(group => subgraph-module.to-label(_rank-subgraph(
      graph,
      group,
    ))),
    route-edge-weight: str(route-edge-weight),
    route-exit-weight: str(route-exit-weight),
    route-label-width-scale: str(route-label-width-scale),
    route-label-width-cap: str(route-label-width-cap),
    depth-scale: str(solver.depth-scale),
    flattening-end: str(solver.flattening-end),
    length-scale: str(spring.at("length", default: length-scale)),
  )
  if subgraph != none {
    settings.insert("subgraph", subgraph-module.to-label(
      subgraph-module._impl.validate(graph, subgraph),
    ))
  }
  let graph-bytes = _plugin.layout_parsed_graph(
    graph-module.graph-bytes(graph),
    cbor.encode(settings),
  )
  graph-module.with-bytes(graph, graph-bytes)
}

/// Apply multiple layout passes in order.
///
/// Each pass is a dictionary of named arguments accepted by `layout`, excluding
/// the graph itself.
/// -> dictionary
#let sequence(
  /// Graph object returned by `graph.build` or `graph.parse`. -> dictionary
  graph,
  /// Array of layout option dictionaries. -> array
  passes,
) = {
  let result = graph
  for pass in passes {
    result = layout(result, ..pass)
  }
  result
}
