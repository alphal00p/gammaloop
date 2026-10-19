# Vakint

Vakint matches, canonicalizes, tensor-reduces, and evaluates single-scale vacuum integrals for
high-energy physics. It combines Symbolica expression handling with analytic FORM-backed methods
and numerical pySecDec evaluation.

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
coefficients are available offline. Additional data and prepared caches live in
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

Prepared caches are optional and read-only during evaluation. See [cache
preparation](cache_preparation.typ), [the cache example](examples/native_cache.rs),
and the [historical program refresh record](rustred_package_refresh.typ).

Licensed Symbolica execution uses `SYMBOLICA_LICENSE` in the environment. Neither
building nor running Vakint requires Mathematica.
