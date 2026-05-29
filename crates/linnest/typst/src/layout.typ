#import "graph.typ" as graph-module
#import "subgraph.typ" as subgraph-module

#let _plugin = plugin("../linnest.wasm")

/// Apply the linnest layout pass to a graph object.
///
/// This is intentionally a second step: construct or parse a graph first, then
/// call `layout`. Set `layout-algo` to `"force"` for deterministic force
/// integration, `"anneal"` for simulated annealing, `"tree"` for a traversal
/// tree placement, `"dot"` for a Graphviz-like layered placement, or
/// `"stable-layered"` for a stable railroad-inspired layered placement.
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
  /// Optional subgraph object to lay out. With `"tree"`, other edges are drawn
  /// from the resulting node positions. With `"dot"` and `"stable-layered"`,
  /// the subgraph determines rank constraints, while all paired edges between
  /// included nodes get dummy routing vertices and edge positions. With
  /// `"force"` and `"anneal"`, nodes and edges outside the subgraph are fixed
  /// boundary points during optimization.
  /// -> none | bytes
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
  /// integration steps; in `"anneal"` mode this is the number of proposals per
  /// temperature epoch. -> int
  steps: int(sys.inputs.at("steps", default: "30")),
  /// Seed for deterministic initialization, force-mode jitter, and annealing
  /// proposals. Applies to both modes. -> int
  seed: int(sys.inputs.at("seed", default: "2")),
  /// Initial movement scale. `"force"` multiplies computed forces by this
  /// value; `"anneal"` uses it as the proposal step size. -> float
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
  /// Force-only early stop threshold for maximum movement in one step. The
  /// annealing schedule stores this value but does not currently use it for
  /// stopping. -> float
  early-tol: 1e-6,
  /// Anneal-only initial temperature used in the Metropolis acceptance test.
  /// -> float
  temp: 0.3,
  /// Force-mode maximum movement clamp per point and per step. -> float
  delta: 0.4,
  /// Base repulsion strength for vertex-vertex interactions. Also scales
  /// `gamma-ev`, `gamma-ee`, `gamma-dangling`, and `g-center`. Applies to both
  /// modes through the shared spring energy.
  /// -> float
  beta: 50.0,
  /// Spring stiffness for node-to-edge incidence lengths. Applies to both
  /// modes. -> float
  k-spring: 11.0,
  /// Centering strength relative to `beta`. Applies to both modes.
  /// -> float
  g-center: 0.002,
  /// Number of epochs. Both modes run up to `steps` iterations inside
  /// each epoch. -> int
  epochs: 30,
  /// Anneal-only fixed energy penalty per detected edge crossing. The direct
  /// force integrator does not currently add a crossing force. -> float
  crossing-penalty: 30.0,
  /// Repulsion for dangling half edges, relative to `beta`. Applies to
  /// both modes through the shared spring energy. -> float
  gamma-dangling: 5.0,
  /// Local edge-edge repulsion, relative to `beta`. Applies to both
  /// modes. -> float
  gamma-ee: 0.1,
  /// Bias that pushes points in directions implied by pin/port constraints.
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
  /// Ordered node indices used as preferred roots for `"tree"`, `"dot"`, and
  /// `"stable-layered"`. Roots outside the selected node set are ignored.
  /// Remaining components are laid out afterward in graph order. -> array
  layout-roots: (),
  /// Subgraphs whose incident nodes should share a dot/stable-layered rank.
  /// These are layout hints supplied by Typst rather than parsed graph
  /// structure. -> array
  rank-same: (),
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
  /// Force-only spring pulling temporary z coordinates back toward the layout
  /// plane. Use with `z-spring-growth` to help separate overlapping
  /// points during integration. -> float
  z-spring: 0.05,
  /// Force-only per-epoch multiplier for `z-spring`. -> float
  z-spring-growth: 1.3,
  /// Natural spring-length multiplier. This scales the graph's preferred edge
  /// length and the repulsive coefficients derived from it. Applies to both
  /// modes. -> float
  length-scale: 0.35,
) = {
  let settings = (
    viewport-w: str(viewport-w),
    viewport-h: str(viewport-h),
    tree-dx: str(tree-dx),
    tree-dy: str(tree-dy),
    steps: str(steps),
    seed: str(seed),
    step: str(step),
    step-shrink: str(step-shrink),
    cool: str(cool),
    accept-floor: str(accept-floor),
    early-tol: str(early-tol),
    temp: str(temp),
    delta: str(delta),
    beta: str(beta),
    k-spring: str(k-spring),
    g-center: str(g-center),
    epochs: str(epochs),
    crossing-penalty: str(crossing-penalty),
    gamma-dangling: str(gamma-dangling),
    gamma-ee: str(gamma-ee),
    directional-force: str(directional-force),
    label-length-scale: str(label-length-scale),
    label-spring: str(label-spring),
    label-charge: str(label-charge),
    label-steps: str(label-steps),
    label-layout: label-layout,
    label-step: str(label-step),
    label-early-tol: str(label-early-tol),
    label-max-delta-scale: str(label-max-delta-scale),
    gamma-ev: str(gamma-ev),
    eps: str(eps),
    incremental-energy: incremental-energy,
    layout-algo: layout-algo,
    layout-nodes: layout-nodes,
    layout-roots: layout-roots,
    rank-same: rank-same.map(subgraph-module.to-label),
    route-edge-weight: str(route-edge-weight),
    route-exit-weight: str(route-exit-weight),
    route-label-width-scale: str(route-label-width-scale),
    route-label-width-cap: str(route-label-width-cap),
    z-spring: str(z-spring),
    z-spring-growth: str(z-spring-growth),
    length-scale: str(length-scale),
  )
  if subgraph != none {
    settings.insert("subgraph", subgraph-module.to-label(subgraph))
  }
  let graph-bytes = _plugin.layout_parsed_graph(graph-module.graph-bytes(graph), cbor.encode(settings))
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
