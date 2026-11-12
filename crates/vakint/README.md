# Vakint

Vakint matches, canonicalizes, tensor-reduces, and evaluates single-scale vacuum integrals for
high-energy physics. It combines Symbolica expression handling with the FORM-independent RustRed
reduction, analytic FORM-backed methods and numerical pySecDec evaluation.

By default Vakint evaluates with RustRed, which needs neither FORM nor pySecDec. When FORM
(4.2.1 or newer) is available, the default evaluation order keeps AlphaLoop, MATAD and FMFT after
RustRed as fallbacks, followed by pySecDec when it can be imported. Explicitly requested
FORM-based backends report a missing FORM installation. Massless propagator integrals
(p-integrals, one to four loops; `docs/p_integrals.md`): RustRed is the default for 1L-4L
p-integrals whenever its library is available (built in for one to three loops; for four loops
the library of vakint-external-data's `p_integral` branch, unreleased until M4, otherwise four
loops stay on Forcer); Forcer is the fallback. The vendored Forcer program
(`EvaluationOrder::forcer_only()`) needs FORM 5 and is distributed with its authors' permission
(see `form_src/forcer/NOTICE`).

Two settings control how the ordered reduction computes. `n_cores` is the single source of every
worker count (RustRed's ordered scheduler and modular replay, the five-loop frontend, the retained
post-processing and the five-loop route verification). It defaults to the cores available to the
process, at most ten, and to one without a Symbolica license. With `n_cores = 1` every computation
runs on the calling thread. `ordered_weight_arithmetic` selects `modular_polynomial` (the default
and the fastest), `modular_sampled` or `exact` frontier weights. Advanced options use RustRed's
syntax, for example `modular_sampled:batch_points=64`. All choices return identical exact results:

```python
vakint = Vakint(n_cores=8, ordered_weight_arithmetic="modular_sampled")
```

The [ordered-reduction notes](docs/ordered_reduction.typ) describe the modes and their options.

Use the canonical Typst documentation for maintained workflows:

- [overview and backend boundaries](https://alphal00p.github.io/gammaloop/products/vakint/latest/)
- [matching-only tutorial](https://alphal00p.github.io/gammaloop/products/vakint/latest/tutorial/)
- [staged evaluation, numerical substitution, result interpretation, and citations](https://alphal00p.github.io/gammaloop/products/vakint/latest/guides/evaluation/)
- [Rust and Python interfaces](https://alphal00p.github.io/gammaloop/products/vakint/latest/reference/interfaces/)

The crate also contains source examples under [`examples/`](examples/); check them against the
versioned API reference before adapting them. Vakint is exposed to Python through
`symbolica.community.hepkit.vakint` when that community module is included in the installed Symbolica
assembly, not as an independent Python distribution.

The optional `feynkit-ingress` feature exposes the native
`VakintExpression::from_diagram` constructor and `DiagramIntegralOptions`.
It consumes FeynKit's graph, integral family and existing momentum routing,
with native Symbolica scalar products and simultaneous parameter substitutions.
The community-module feature enables this automatically and exports the same
operation as `integral_from_diagram`; no Python-side graph or algebra adapter is
involved. Neither entry point rematches the graph or changes evaluation defaults.

RustRed provides FORM-independent scalar reduction. Its 1L–3L programs and master
coefficients are available offline. Additional data live in
[vakint-external-data](https://github.com/alphal00p/vakint-external-data).
Explicitly loading a maximum loop count installs every lower loop count as well:

```python
engine.load_official_external_data(4, local_artifacts_path="./data")
```

A local external-data checkout can be passed as the same path; it is validated
and loaded without downloading or changing its contents. Installed releases are
reused offline until an explicit refresh. See [external-data loading and master
capabilities](external_data.typ) for the Rust API, release selection, and errors.
[Master definitions and capabilities](masters.typ) records the normalization,
finite coefficient ranges, and distinction between terminals and computed masters.

RustRed's `master_mode` selects `retain`, `symbolic`, or `numerical` (default).
Symbolic expansion retains labeled coefficients when only numerical data exists.
Numerical evaluation rejects insufficient coefficient precision or Laurent depth,
and lists unresolved integrals explicitly. Candidate campaigns can return partial
reductions with warnings; their terminal inventories are not complete-basis proofs.

FORM reductions and their master formulas remain available. HPL numerical
substitution uses Symbolica's arbitrary-precision sixth-root algorithms in both
FORM and RustRed paths. Other finite FORM tables retain their existing precision
warnings. Working precision alone does not certify the accuracy of a final result.

Vakint reduces with RustRed's ordered scheduler only: all contributions of a request share one
exact frontier per compatible coefficient channel, with exact or modular weight arithmetic. See
[ordered reduction](docs/ordered_reduction.typ). Depth-first reduction and its completed-result
caches were removed from Vakint; prepared depth-first caches shipped by older external-data
releases are ignored. Instead, an Ordered step cache records each expanded integral's checked
reduction step once and replays it exactly (never weights or results), in memory by default or
in a directory that later processes reuse (`VakintSettings::ordered_cache`; Python
`ordered_cache*` keywords and `load_`/`save_`/`warm_`/`clear_ordered_cache`).

Licensed Symbolica execution uses `SYMBOLICA_LICENSE` in the environment. Neither
building nor running Vakint requires Mathematica.

## Testing

Two recipes run the Vakint tests from the workspace root. Both use CI's build
(`ci-optim`, feature `feynkit-ingress`) and pass extra arguments to nextest.

- `just test_vakint` is the fast local gate (nextest profile `test_vakint`,
  under ten minutes with 16 test threads). It covers the vendored 1L–3L
  evaluation with the MATAD/FMFT/alphaLoop/RustRed comparisons, the ordered
  reduction rungs through three loops, the Ordered step and shipped caches,
  external-data loading and validation, and the monitor. The Python stub test
  needs `--features vakint/python_stubgen,vakint/symbolica_community_module`
  (a separate build; CI checks the stubs with the docs Python exporter).
- `just test_vakint_full` runs exactly what CI runs for Vakint (nextest profile
  `ci_gammaloop`, ignored tests excluded). It adds the four-loop rungs and the
  five-loop closure and acceptance tests, and like CI it runs the Vakint tests
  one at a time (about 20 minutes).

CI sets `VAKINT_TEST_EXTERNAL_DATA` to the vakint-external-data commit pinned as
`vakintExternalTestData` in `nix/rust-workspace.nix`; set it to a checkout of
that commit to match CI. Without it, tests that need 4L/5L data install the
release named in `data/manifest.json`. pySecDec comparisons run only with
`RUN_PYSECDEC_TESTS=1`. Ignored tests are slow fixtures, timings and offline
FORM oracles; each `#[ignore]` reason says what it needs.
