#import "graph.typ" as graph-module
#import "subgraph.typ" as subgraph-module
#import "curve.typ" as curve-module
#import "draw.typ": draw

#let _plugin = plugin("./linnest.wasm")

/// Apply the linnest layout pass to a graph object.
///
/// This is intentionally a second step: construct or parse a graph first, then
/// call `layout`. Set `layout_algo` to `"force"` for deterministic force
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
  viewport_w: 10.0,

  /// Height of the layout viewport used to derive the natural spring length.
  /// Applies to both `"force"` and `"anneal"`. -> float
  viewport_h: 10.0,

  /// Horizontal spacing multiplier for the initial traversal-tree placement.
  /// Applies before both `"force"` and `"anneal"`. -> float
  tree_dx: 0.9,

  /// Vertical spacing multiplier for the initial traversal-tree placement.
  /// Applies before both `"force"` and `"anneal"`. -> float
  tree_dy: 1.2,

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
  /// falls below `accept_floor`. -> float
  step_shrink: 0.21,

  /// Cooling factor applied once per epoch. `"anneal"` cools `temp`; `"force"`
  /// shrinks `step`. -> float
  cool: 0.85,

  /// Anneal-only acceptance-ratio threshold below which `step` is shrunk by
  /// `step_shrink`. -> float
  accept_floor: 0.15,

  /// Force-only early stop threshold for maximum movement in one step. The
  /// annealing schedule stores this value but does not currently use it for
  /// stopping. -> float
  early_tol: 1e-6,

  /// Anneal-only initial temperature used in the Metropolis acceptance test.
  /// -> float
  temp: 0.3,

  /// Force-mode maximum movement clamp per point and per step. -> float
  delta: 0.4,

  /// Base repulsion strength for vertex-vertex interactions. Also scales
  /// `gamma_ev`, `gamma_ee`, `gamma_dangling`, and `g_center`. Applies to both
  /// modes through the shared spring energy.
  /// -> float
  beta: 50.0,

  /// Spring stiffness for node-to-edge incidence lengths. Applies to both
  /// modes. -> float
  k_spring: 11.0,

  /// Centering strength relative to `beta`. Applies to both modes.
  /// -> float
  g_center: 0.005,

  /// Number of epochs. Both modes run up to `steps` iterations inside
  /// each epoch. -> int
  epochs: 30,

  /// Anneal-only fixed energy penalty per detected edge crossing. The direct
  /// force integrator does not currently add a crossing force. -> float
  crossing_penalty: 30.0,

  /// Repulsion for dangling half edges, relative to `beta`. Applies to
  /// both modes through the shared spring energy. -> float
  gamma_dangling: 5.0,

  /// Local edge-edge repulsion, relative to `beta`. Applies to both
  /// modes. -> float
  gamma_ee: 0.1,

  /// Bias that pushes points in directions implied by pin/port constraints.
  /// Applies to both modes. -> float
  directional_force: 5.0,

  /// Edge-label target offset as a multiple of the graph spring length. Label
  /// layout runs after both graph layout modes. -> float
  label_length_scale: 0.6,

  /// Spring strength pulling each label toward its target offset. Label layout
  /// runs after both modes. -> float
  label_spring: 23.0,

  /// Repulsion strength between labels and graph points, scaled by spring
  /// length squared. Label layout runs after both modes. -> float
  label_charge: 3.0,

  /// Maximum number of post-layout label relaxation steps. Set to `0` to
  /// skip label placement. Applies after both modes. -> int
  label_steps: 20,

  /// Label relaxation step size. Applies after both modes. -> float
  label_step: 0.15,

  /// Label relaxation early stop threshold. Applies after both modes. -> float
  label_early_tol: 1e-3,

  /// Label movement clamp as a multiple of spring length. Applies after both
  /// modes. -> float
  label_max_delta_scale: 0.5,

  /// Edge-vertex repulsion, relative to `beta`. Applies to both modes.
  /// -> float
  gamma_ev: 0.01,

  /// Softening epsilon used in inverse-square force/energy terms. Applies to
  /// both modes. -> float
  eps: 1e-4,

  /// Whether annealing updates cached energy by local deltas. This is
  /// anneal-only; force mode computes direct forces instead of energies.
  /// -> bool
  incremental_energy: true,

  /// Layout algorithm. Use `"force"` for direct force integration or `"anneal"`
  /// for simulated annealing against the spring energy. -> string
  layout_algo: "force",

  /// Force-only spring pulling temporary z coordinates back toward the layout
  /// plane. Use with `z_spring_growth` to help separate overlapping
  /// points during integration. -> float
  z_spring: 0.05,

  /// Force-only per-epoch multiplier for `z_spring`. -> float
  z_spring_growth: 1.3,

  /// Natural spring-length multiplier. This scales the graph's preferred edge
  /// length and the repulsive coefficients derived from it. Applies to both
  /// modes. -> float
  length_scale: 0.5,
) = {
  let settings = (
    viewport_w: str(viewport_w),
    viewport_h: str(viewport_h),
    tree_dx: str(tree_dx),
    tree_dy: str(tree_dy),
    steps: str(steps),
    seed: str(seed),
    step: str(step),
    step_shrink: str(step_shrink),
    cool: str(cool),
    accept_floor: str(accept_floor),
    early_tol: str(early_tol),
    temp: str(temp),
    delta: str(delta),
    beta: str(beta),
    k_spring: str(k_spring),
    g_center: str(g_center),
    epochs: str(epochs),
    crossing_penalty: str(crossing_penalty),
    gamma_dangling: str(gamma_dangling),
    gamma_ee: str(gamma_ee),
    directional_force: str(directional_force),
    label_length_scale: str(label_length_scale),
    label_spring: str(label_spring),
    label_charge: str(label_charge),
    label_steps: str(label_steps),
    label_step: str(label_step),
    label_early_tol: str(label_early_tol),
    label_max_delta_scale: str(label_max_delta_scale),
    gamma_ev: str(gamma_ev),
    eps: str(eps),
    incremental_energy: incremental_energy,
    layout_algo: layout_algo,
    z_spring: str(z_spring),
    z_spring_growth: str(z_spring_growth),
    length_scale: str(length_scale),
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
