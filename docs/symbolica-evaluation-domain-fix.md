# Symbolica evaluation-domain bug and fix

After replacing exponential workarounds with Symbolica's native hyperbolic
functions, GammaLoop's interpreted evaluator could panic during thermal CFF
evaluation, including `thermal_bugblatter`:

```text
External function 'symbolica_coth' does not have an implementation for gammalooprs::utils::F<f64>
```

Symbolica resolves these functions through numeric callbacks registered for an
exact Rust type. A callback for `f64` therefore does not automatically apply to
GammaLoop's `F<f64>`. The wrappers implemented constant conversion through
`try_from_complex_float`, but lacked callback forwarding through
`EvaluationDomain::resolve_function`. Implementing numeric methods such as
`tanh` on a wrapper was insufficient for this dispatch path. The panic came from
missing callback forwarding in GammaLoop and its bundled Spenso wrappers.

Commit `951972fd3` adds callback resolution to the existing domain implementations:

| Wrapper | Callback delegation |
| --- | --- |
| `F<T>` | Unwrap arguments, resolve through `T`, and wrap the result. |
| `QuadFloat` | Delegate through Symbolica's `DoubleFloat` domain. |
| `VarFloat<N>` | Convert through Symbolica's arbitrary-precision `Float`. |
| Spenso `Complex<T>` | Prefer the corresponding Symbolica complex callback; otherwise evaluate through `Complex<Float>`. |

Callbacks registered directly for a wrapper retain priority. Symbolic tags and
all numeric arguments are forwarded, and higher-precision values avoid an
intermediate conversion through `f64`. Complex callbacks evaluate the complete
complex value. The Spenso implementation adds an `Into<Float>` bound, satisfied
by GammaLoop's existing numeric wrappers. This fix uses the pinned Symbolica API.

Eight new regression tests cover native `tanh`, `coth`, `sech`, and `csch`
evaluation, complex values, precision beyond `f64`, callback precedence, tags,
multiple arguments, and unresolved functions. Together with nine thermal CFF
tests, all **17 tests passed with retries disabled**. Formatting, `cargo check`,
and Clippy also passed.

Testing exposed an intermittent `thermal_ring` relative error of approximately
`2.02e-11`, just above its previous `2e-11` limit. With user approval, that
tolerance was increased to `3e-11`. Cancellation and varying evaluation order
are plausible contributors; the underlying numerical sensitivity remains to be
investigated.

The separate Symbolica C++ export/header issue was outside this runtime fix.
Switching compilation to symjit does not supply callbacks for GammaLoop's
interpreted numeric domains.

The targeted test run was:

```sh
cargo nextest run -p gammalooprs -p spenso --lib \
  --cargo-profile dev-optim --locked --offline -P test_gammaloop \
  --retries 0 \
  -E 'test(evaluation_domain) or (test(cff::generation::tests_cff::thermal) and not test(failing))'
```
