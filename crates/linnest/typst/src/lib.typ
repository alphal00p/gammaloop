#import "graph.typ" as graph-module
#import "subgraph.typ" as subgraph-module

/// Default layout configuration used by @layout.
///
/// These are the defaults for @layout's named settings. Numeric values are
/// native Typst integers/floats here and are converted to strings only when the
/// wrapper sends settings to the plugin.
/// -> dictionary
#let config = (
  viewport_w: 10.0,
  viewport_h: 10.0,
  tree_dx: 0.9,
  tree_dy: 1.2,
  steps: int(sys.inputs.at("steps", default: "15")),
  seed: int(sys.inputs.at("seed", default: "14")),
  step: 0.81,
  step_shrink: 0.21,
  cool: 0.85,
  accept_floor: 0.15,
  early_tol: 1e-6,
  temp: 0.3,
  delta: 0.4,
  beta: 46.1,
  k_spring: 11.0,
  g_center: 40.0,
  epochs: 30,
  crossing_penalty: 30.0,
  gamma_dangling: 40.0,
  gamma_ee: 0.1,
  directional_force: 5.0,
  label_length_scale: 0.6,
  label_spring: 23.0,
  label_charge: 3.0,
  label_steps: 20,
  label_step: 0.15,
  label_early_tol: 1e-3,
  label_max_delta_scale: 0.5,
  gamma_ev: 0.1,
  eps: 1e-4,
  incremental_energy: true,
  layout_algo: "force",
  z_spring: 0.05,
  z_spring_growth: 1.3,
  length_scale: 0.1,
)

#let _plugin = plugin("../linnest.wasm")

/// Apply the linnest layout pass to an archived graph.
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
/// #let g = layout(g, seed: "2", steps: "5")
/// #graph.edges(g).map(edge => edge.pos)
/// ```
/// -> bytes
#let layout(
  /// Archived graph bytes returned by `graph.build`, `graph.finish`, or
  /// `graph.parse`.
  /// -> bytes
  graph,

  /// Width of the layout viewport used to derive the natural spring length.
  /// Applies to both `"force"` and `"anneal"`. -> float
  viewport_w: config.viewport_w,

  /// Height of the layout viewport used to derive the natural spring length.
  /// Applies to both `"force"` and `"anneal"`. -> float
  viewport_h: config.viewport_h,

  /// Horizontal spacing multiplier for the initial traversal-tree placement.
  /// Applies before both `"force"` and `"anneal"`. -> float
  tree_dx: config.tree_dx,

  /// Vertical spacing multiplier for the initial traversal-tree placement.
  /// Applies before both `"force"` and `"anneal"`. -> float
  tree_dy: config.tree_dy,

  /// Iterations per epoch. In `"force"` mode this is the number of force
  /// integration steps; in `"anneal"` mode this is the number of proposals per
  /// temperature epoch. -> int
  steps: config.steps,

  /// Seed for deterministic initialization, force-mode jitter, and annealing
  /// proposals. Applies to both modes. -> int
  seed: config.seed,

  /// Initial movement scale. `"force"` multiplies computed forces by this
  /// value; `"anneal"` uses it as the proposal step size. -> float
  step: config.step,

  /// Anneal-only step shrink factor, applied when an epoch's acceptance ratio
  /// falls below `accept_floor`. -> float
  step_shrink: config.step_shrink,

  /// Cooling factor applied once per epoch. `"anneal"` cools `temp`; `"force"`
  /// shrinks `step`. -> float
  cool: config.cool,

  /// Anneal-only acceptance-ratio threshold below which `step` is shrunk by
  /// `step_shrink`. -> float
  accept_floor: config.accept_floor,

  /// Force-only early stop threshold for maximum movement in one step. The
  /// annealing schedule stores this value but does not currently use it for
  /// stopping. -> float
  early_tol: config.early_tol,

  /// Anneal-only initial temperature used in the Metropolis acceptance test.
  /// -> float
  temp: config.temp,

  /// Force-mode maximum movement clamp per point and per step. -> float
  delta: config.delta,

  /// Base repulsion strength for vertex-vertex interactions. Also scales
  /// `gamma_ev`, `gamma_ee`, `gamma_dangling`, and `g_center`. Applies to both
  /// modes through the shared spring energy.
  /// -> float
  beta: config.beta,

  /// Spring stiffness for node-to-edge incidence lengths. Applies to both
  /// modes. -> float
  k_spring: config.k_spring,

  /// Centering strength relative to `beta`. Applies to both modes.
  /// -> float
  g_center: config.g_center,

  /// Number of epochs. Both modes run up to `steps` iterations inside
  /// each epoch. -> int
  epochs: config.epochs,

  /// Anneal-only fixed energy penalty per detected edge crossing. The direct
  /// force integrator does not currently add a crossing force. -> float
  crossing_penalty: config.crossing_penalty,

  /// Repulsion for dangling half edges, relative to `beta`. Applies to
  /// both modes through the shared spring energy. -> float
  gamma_dangling: config.gamma_dangling,

  /// Local edge-edge repulsion, relative to `beta`. Applies to both
  /// modes. -> float
  gamma_ee: config.gamma_ee,

  /// Bias that pushes points in directions implied by pin/port constraints.
  /// Applies to both modes. -> float
  directional_force: config.directional_force,

  /// Edge-label target offset as a multiple of the graph spring length. Label
  /// layout runs after both graph layout modes. -> float
  label_length_scale: config.label_length_scale,

  /// Spring strength pulling each label toward its target offset. Label layout
  /// runs after both modes. -> float
  label_spring: config.label_spring,

  /// Repulsion strength between labels and graph points, scaled by spring
  /// length squared. Label layout runs after both modes. -> float
  label_charge: config.label_charge,

  /// Maximum number of post-layout label relaxation steps. Set to `"0"` to
  /// skip label placement. Applies after both modes. -> int
  label_steps: config.label_steps,

  /// Label relaxation step size. Applies after both modes. -> float
  label_step: config.label_step,

  /// Label relaxation early stop threshold. Applies after both modes. -> float
  label_early_tol: config.label_early_tol,

  /// Label movement clamp as a multiple of spring length. Applies after both
  /// modes. -> float
  label_max_delta_scale: config.label_max_delta_scale,

  /// Edge-vertex repulsion, relative to `beta`. Applies to both modes.
  /// -> float
  gamma_ev: config.gamma_ev,

  /// Softening epsilon used in inverse-square force/energy terms. Applies to
  /// both modes. -> float
  eps: config.eps,

  /// Whether annealing updates cached energy by local deltas. This is
  /// anneal-only; force mode computes direct forces instead of energies.
  /// -> bool
  incremental_energy: config.incremental_energy,

  /// Layout algorithm. Use `"force"` for direct force integration or `"anneal"`
  /// for simulated annealing against the spring energy. -> string
  layout_algo: config.layout_algo,

  /// Force-only spring pulling temporary z coordinates back toward the layout
  /// plane. Use with `z_spring_growth` to help separate overlapping
  /// points during integration. -> float
  z_spring: config.z_spring,

  /// Force-only per-epoch multiplier for `z_spring`. -> float
  z_spring_growth: config.z_spring_growth,

  /// Natural spring-length multiplier. This scales the graph's preferred edge
  /// length and the repulsive coefficients derived from it. Applies to both
  /// modes. -> float
  length_scale: config.length_scale,
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
/// Graph values are archived byte arrays. Keep them opaque and pass them back to
/// this module, @subgraph, or @layout.
/// -> module
#let graph = graph-module

/// Subgraph namespace.
///
/// Subgraph values are archived byte arrays. They can be passed to
/// `graph.nodes` and `graph.edges`.
/// -> module
#let subgraph = subgraph-module
