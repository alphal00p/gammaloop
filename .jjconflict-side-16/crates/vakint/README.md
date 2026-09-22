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
`symbolica.community.vakint` when that community module is included in the installed Symbolica
assembly, not as an independent Python distribution.

Symbolica has its own license terms. FORM, MATAD, FMFT, and pySecDec are required only by the
selected evaluation path; matching can use an empty evaluation order. Contributor policy is in
[`CONTRIBUTING.typ`](../../CONTRIBUTING.typ).
