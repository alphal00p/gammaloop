#import "graph.typ" as graph-module
#import "subgraph.typ" as subgraph-module
#import "curve.typ" as curve-module
#import "draw.typ": draw

#let _plugin = plugin("./linnest.wasm")

/// Apply the linnest layout pass to a graph object.
///
/// This is intentionally a second step: construct or parse a graph first, then
/// call `layout`. Set `layout-algo` to `"force"` for deterministic force
/// integration or `"anneal"` for simulated annealing.
///
/// ```example
/// #let g = graph.build(
///   nodes: ((name: "a"), (name: "b")),
///   edges: ((source: (node: 0), sink: (node: 1)),),
/// )
/// #let g = layout(g)
/// #graph.edges(g).map(edge => edge.pos)
/// ```
/// -> bytes
#let layout(
  /// Graph object returned by `graph.build`, `graph.finish`, or
  /// `graph.parse`.
  /// -> bytes
  graph,

  /// Width of the layout viewport used to derive the natural spring length.
  /// Applies to both `"force"` and `"anneal"`. -> float
  viewport-w: 10.0,

  /// Height of the layout viewport used to derive the natural spring length.
  /// Applies to both `"force"` and `"anneal"`. -> float
  viewport-h: 10.0,

  /// Horizontal spacing multiplier for the initial traversal-tree placement.
  /// Applies before both `"force"` and `"anneal"`. -> float
  tree-dx: 0.9,

  /// Vertical spacing multiplier for the initial traversal-tree placement.
  /// Applies before both `"force"` and `"anneal"`. -> float
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
  g-center: 0.005,

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

  /// Spring strength pulling each label toward its target offset. Label layout
  /// runs after both modes. -> float
  label-spring: 23.0,

  /// Repulsion strength between labels and graph points, scaled by spring
  /// length squared. Label layout runs after both modes. -> float
  label-charge: 3.0,

  /// Maximum number of post-layout label relaxation steps. Set to `0` to
  /// skip label placement. Applies after both modes. -> int
  label-steps: 20,

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

  /// Layout algorithm. Use `"force"` for direct force integration or `"anneal"`
  /// for simulated annealing against the spring energy. -> string
  layout-algo: "force",

  /// Force-only spring pulling temporary z coordinates back toward the layout
  /// plane. Use with `z-spring-growth` to help separate overlapping
  /// points during integration. -> float
  z-spring: 0.05,

  /// Force-only per-epoch multiplier for `z-spring`. -> float
  z-spring-growth: 1.3,

  /// Natural spring-length multiplier. This scales the graph's preferred edge
  /// length and the repulsive coefficients derived from it. Applies to both
  /// modes. -> float
  length-scale: 0.5,
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
    label-step: str(label-step),
    label-early-tol: str(label-early-tol),
    label-max-delta-scale: str(label-max-delta-scale),
    gamma-ev: str(gamma-ev),
    eps: str(eps),
    incremental-energy: incremental-energy,
    layout-algo: layout-algo,
    z-spring: str(z-spring),
    z-spring-growth: str(z-spring-growth),
    length-scale: str(length-scale),
  )
  _plugin.layout_parsed_graph(bytes(graph), cbor.encode(settings))
}

/// Graph namespace.
///
/// Graph objects are opaque zero-copy values. Keep them opaque and pass them back to
/// this module, @subgraph, or @layout.
/// -> module
#let graph = graph-module

/// Curve namespace.
///
/// Helpers for splitting and drawing Bezier edge geometry.
/// -> module
#let curve = curve-module

/// Subgraph namespace.
///
/// Subgraph objects are opaque zero-copy values. They can be passed to
/// `graph.nodes` and `graph.edges`.
/// -> module
#let subgraph = subgraph-module
