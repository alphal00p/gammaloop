= Factorized numerator restoration validation

Measured on 2026-09-21 at `20974534251153e805c00ed017dfd874a29626d5`, using native
binary SHA256 `415c94c32ac8297ac22fca1169fc5d5e5a19df7c8b48d236dd27d85eb300d874`.

== Representation

Numerators remain shared functions of affine energy arguments, including signs,
integer sampling nodes, and external-energy coefficients. Complete energy maps
define branches; arguments may span the whole graph. This retains the structure
needed for future LTD maps without implementing LTD.

`NumeratorAtomExt::series_preserving_factors` retains expansion-independent
factors and shared numerator coefficient families. Expansion-dependent energy
arguments remain visible to differentiation; Symbolica supplies the Laurent
algebra. Both local UV routes use this machinery.

Evaluator preparation contracts closed parametric numerator contexts before
component enumeration. Closed sums are prepared from the inside out and remain
opaque in their parent products. Shared bodies retain factorized numerator sums,
separate formals and private indices per occurrence, and existing chain/trace
scopes. Scalar weights, inverse denominators, selectors, and independent closed
products keep their boundaries. Open interfaces use ordinary component
preparation; only reachable functions and scalar aliases are retained.

Integrated UV preparation reduces color algebra before lowering chain/trace
shorthands for Vakint. This removes reducible traces before repeated dummy-index
materialization, without globally expanding the numerator.

== Comparison protocol

Live main was checked as `312cf6aaf4cd1414f6742ad87ce39922c497a0b6` immediately
before the runs. Its reused native binary was built from
`226b7bea45398023f767fb4e53a1b5c6c3a43a2d`, through unchanged child checkout
`6fbfe7cd89839ca6b9f10bdd1e5bbc9eba1152f1`. Rust sources, assets, Cargo manifests,
lockfile, and toolchain match main312; the eight intervening changed files concern
CI and contributor guidance. Timing measurements were fresh.

Both sides use Rust 1.98.1, Symbolica 3.0.0, native `dev-optim`, one generation
core, and no evaluator compilation. Main uses `sparse_atom_aware`; the candidate
uses `intermediate_cost`. Binary, card, graph, and effective-setting identities
were checked. Runs were sequential with summary logging; this task ran no concurrent builds
or tests. Each row is one observation on a shared host.

== GL1 generation

The direct-3D card keeps local UV subtraction enabled and toggles integrated UV.
Durations are seconds. Peak memory is native process RSS sampled every 100 ms,
not the separate watchdog's process-tree measurement. Symbolica time means
evaluator construction.

#table(
  columns: (1.8fr, 1fr, 1fr, 1fr, 1fr, 1fr),
  table.header([Case], [Expression], [Spenso], [Symbolica], [Total], [Peak GiB]),
  [Main, integrated off], [7.387], [11.389], [26.741], [45.517], [2.653],
  [Candidate, integrated off], [2.732], [4.464], [0.831], [8.027], [0.191],
  [Main, integrated on], [24.536], [15.599], [31.382], [71.517], [2.849],
  [Candidate, integrated on], [12.322], [10.946], [2.745], [26.013], [0.215],
)

Native generation is 5.67 times faster with integrated UV off and 2.75 times
faster with it on. These derivative-UV timings are not equal-physics numerical
comparisons with main, which lacks the required higher numerator energy powers.

== Numerical controls

Eligible comparisons with main use GL1 without UV and the logarithmic-only GL0
top-quark photon box. The latter has UV degree zero, no proper loop subgraphs,
and no Taylor derivatives. Graphs, runtime settings, and model parameters match;
each side produced two identical, finite, nonzero inspections at the fixed point.

#table(
  columns: (2fr, 1fr, 1fr, 1fr, 1fr),
  table.header([Control], [Main native s], [Candidate native s], [Main wall s], [Candidate wall s]),
  [GL1, no UV], [0.281], [0.086], [1.291], [1.285],
  [GL0, logarithmic UV], [4.166], [0.929], [5.165], [2.119],
)

Raw complex results, written as `(real, imaginary)`, are:

- GL1: main `(0, 1.818545415743298e-6)`, candidate `(0, 1.818545415743152e-6)`.
- GL0: main `(-2.4834308873547831e-29, -8.71033584620754e-10)`, candidate
  `(-1.6648052328444753e-29, 8.71033584620754e-10)`.

The predeclared phase `(-1)^L` is +1 for GL1 and -1 for GL0. Raw relative
differences are `8.03e-14` and `2.0`; phase-adjusted differences are `8.03e-14`
and `4.76e-20`, both below `1e-12`. GL0's opposite raw sign is retained explicitly.
Native peak RSS rises slightly in these controls: 0.0396/0.0433 GiB for
main/candidate GL1 and 0.0801/0.0973 GiB for GL0.

Same-physics derivative GL1 inspections against raised-energy revision
`11fe63d82a13abc0fc16c2fc6296c880a52c0572` agree to `3.4933e-16` with integrated
UV off and on, using phase +1. Both repeat inspections agree, inputs match,
and saved states remain unchanged. This is a single-point sanity check.

== Projected local-4D and validation status

The established projected-4D card, with integrated UV off and explicit orientation
summation, takes 12.337 s native generation, 13.658 s command wall, and 1.290 GiB
native peak RSS. Its expression/Spenso/Symbolica phases are 6.600/4.336/1.401 s.
Main lacks this route and explicit-sum option; its 45.517 s direct-3D reference
above is a generation-cost comparison, not equal-route numerical validation.

Local formatting and Clippy passed, together with 1,005 selected Idenso/GammaLoop
unit tests; 84 tests were skipped. Final review additionally requires a successful
`just ci-checks-and-upload` run before publishing CI-enabled work; the PR records
that result. These measurements establish no numerical evaluator
throughput claim or isolated speedup attributable to any single preparation step.
