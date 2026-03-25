# Sampling Settings Structure

## Summary

The runtime settings tree is rooted at `RuntimeSettings`, and `sampling` is one branch under it alongside `general`, `kinematics`, `integrator`, `observables`, and the other runtime domains. The REPL and settings completion UX are generated directly from the serialized default `RuntimeSettings` value and schema, so the Rust type layout is also the user-facing settings tree layout. There is no separate hand-written tree model.

For `sampling`, the current structure is:

```text
sampling:
  type = default | multi_channeling | discrete_graph_sampling

  if type=default:
    mode, mapping, b

  if type=multi_channeling:
    alpha
    parameterization_settings:
      mode, mapping, b

  if type=discrete_graph_sampling:
    sample_orientations
    sampling_type:
      subtype = default | multi_channeling | discrete_multi_channeling | tropical

      if subtype=default:
        mode, mapping, b

      if subtype=multi_channeling or discrete_multi_channeling:
        alpha
        parameterization_settings:
          mode, mapping, b

      if subtype=tropical:
        upcast_on_failure
        matrix_stability_test
```

Run-card examples in `tests/resources/run_cards/test_qqx_aaa_pentagon_generate_dario.toml` reflect this intended shape.

## Why The Structure Exists

Some of the split is real and not just accidental configuration complexity.

- The outer split is meaningful. `Default` and `MultiChanneling` have no discrete graph axis, while `DiscreteGraphs` does. That changes sample topology, validation, and grid construction.
- The inner split is also partly meaningful. `tropical` is not just a parameter tweak. It changes how x-space points are mapped into loop momenta, and it also changes how the continuous dimension is computed.
- `discrete_multi_channeling` is structurally different from `multi_channeling`: it introduces another discrete axis for channel choice, which changes the discrete depth and the nested grid shape.
- The sample data structures mirror this branching, so the tree is reflected all the way into runtime sample representation rather than being just a parsing artifact.

## Assessment

The current design is partly structurally needed and partly more nested than it has to be.

- A dedicated `sampling` branch is needed.
- The outer distinction between continuous-only modes and discrete-graph modes is mostly needed.
- The inner `sampling_type` nesting under `discrete_graph_sampling` is not fully needed as a user-facing structure.

The current inner nesting mixes three separate concerns:

1. graph-group selection strategy,
2. orientation and channel discrete axes,
3. continuous coordinate backend (`standard` versus `tropical`).

That makes the settings shape harder to understand than the runtime problem actually is.

## Simpler Flow

A simpler public configuration could make these concerns explicit and orthogonal, then normalize them into one internal runtime plan:

```text
sampling:
  graphs: summed | monte_carlo
  orientations: summed | monte_carlo
  channels:
    mode: summed | weighted_sum | monte_carlo
    alpha: ...
  coordinates:
    kind: standard | tropical
    if standard: mode, mapping, b
    if tropical: upcast_on_failure, matrix_stability_test
```

The internal normalized plan would then answer:

- how many discrete axes exist,
- what each axis means,
- how continuous coordinates are parameterized,
- whether tropical rules apply.

## Practical Recommendation

If the goal is low-risk cleanup, keep the current external schema but add an internal normalized `SamplingPlan` or `SamplingLayout`, and make the runtime branch on that instead of repeatedly matching the nested config enums in multiple places.

If the goal is also to simplify the public settings model, split coordinate parameterization from discrete sampling axes. That is the more natural conceptual boundary.

## Bottom Line

The current tree is not arbitrary, but it is not the simplest possible decomposition either. The runtime really does need:

- a distinction between continuous-only and discrete-graph sampling,
- a special tropical path,
- and explicit handling of extra discrete axes such as orientation or channel sampling.

What is not structurally required is expressing all of that as a nested sum-of-sums in the user-facing settings tree.
