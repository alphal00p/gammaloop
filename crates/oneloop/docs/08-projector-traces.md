# Chirality-projector traces — closed fermion loops now collapse

This documents the change on this branch that makes **closed-fermion-loop
numerators fully evaluate** in the app, instead of leaving an inert Dirac trace
`Tr(γ … ℙ_p … γ)` in the reduced numerator. It affects any process with a
closed fermion loop carrying a chirality projector or a `γ5` — `gg→h`, `H→γγ`,
and everything electroweak.

The change is **entirely in `idenso`** (the Dirac/tensor-algebra crate). It flows
to the app through the existing `simplify_gamma()` call in
`gammalooprs`; **`spenso` is untouched** and no public API changed.

## The symptom

Before this branch, reducing `gg→h` in the app produced a `reduced_num` whose
coefficients still contained an *unevaluated* closed-loop trace, e.g.

```
… Tr(γ^α ℙ_p^(in out) γ^α1 γ^(e3) γ^(e4)) p₀ k₀ ϵ₀ ϵ₁ …
```

The top-quark loop's Dirac trace never collapsed to a Lorentz structure, so the
numerator wasn't a clean `Σ coeff × master` — it was blocked at the
gamma-algebra step. (This was the "closed loops show an unevaluated `Tr(...)`
gamma-algebra frontier" noted in the app coverage.)

## Root cause (`crates/idenso/src/dirac/simplify.rs`)

`simplify_trace_node` evaluates traces of pure gammas, and has dedicated
handling for `γ5` and `γ0`, **but a trace containing a chirality projector
`ℙ±` is never decomposed**. When it reaches
`gamma_mink_index_sequence_for`, the projector is not a gamma, that helper
returns `None`, and the *entire* trace is left inert. A closed fermion loop that
`collect_gamma_chains` had correctly closed into `trace(cyclic(ℙ_p, γ, γ, γ, γ))`
therefore stayed symbolic forever.

(A false lead worth recording: a `parse_lit!("gamma", default_namespace="spenso")`
reproduction created a gamma symbol that is **not** `AGS.gamma` — which is a
`tensor_symbol!` — so even a pure-`γ5` loop "failed to evaluate." That was a test
artifact, not the bug. Tests must build gammas with the `gamma!` macro and
projectors with `FunctionBuilder::new(AGS.projp)`.)

## The fix — three surgical additions to `dirac/simplify.rs`

1. **`expand_chirality_projector`** (chain nodes) and
   **`expand_chirality_projector_trace`** (trace nodes) — Rust-side rules
   **guarded on four-dimensional bispinor endpoints**, mirroring every other
   dimension-sensitive rule in this file:

   ```
   ℙ₊(a,b) → ½ ( δ(a,b) + γ5(a,b) )
   ℙ₋(a,b) → ½ ( δ(a,b) − γ5(a,b) )
   ```

   where `δ = id_atom` is the bispinor identity linking the endpoints. This is
   the textbook definition of the chirality projectors; expressing them in the
   gamma basis lets the existing pure-gamma and `γ5` trace machinery finish the
   job. The projector-free (`δ`) part and the `γ5` part are both things the
   evaluator already knows how to reduce.

   The **four-dimensional guard is load-bearing**: every downstream `γ5` rule
   (anticommutation, adjacent-pair contraction, trace recursion) is itself
   4D-only, so expanding a projector on `bis(d, …)` endpoints would strand an
   irreducible `γ5` and leave the expression *worse* than the inert `projp` it
   replaced. In `d` dimensions the projector is therefore left untouched, and
   output is byte-identical to before this change.

2. **`chain(i,i) → trace`** in `simplify_chain_node` — a chain whose two
   endpoints coincide *is* a trace over that bispinor index, so it is
   re-expressed as `trace!(rep; factors)` and the trace evaluator takes over.
   This fires only for a genuine equal-endpoint bispinor chain, so open chains
   are untouched, and it is gated on `settings.evaluate_traces` so
   `without_trace_evaluation()` keeps its historical output shape.

   Note: now that the expansion in (1) splices the projector out *in place*
   rather than injecting a free `δ`, chain endpoints no longer change, and
   `collect_gamma_chains` already closes genuine loops itself. This rule is
   consequently **unexercised by the current test suite** — it is retained as
   defensive, and is a candidate for removal.

3. **`bispinor_rep_of_slot`** — small helper recovering the representation
   `bis(dim)` from a slot `bis(dim, index)`, needed to build the trace in (2).

Scope of the behaviour change — stated honestly, because this is a shared
crate. In `d` dimensions and on projector-free input, every gamma/trace path is
byte-identical. **In four dimensions it is not**: a `projp`/`projm` on a 4D
fermion line now leaves `simplify_gamma` as a two-term `γ5` polynomial instead
of a compact projector factor. That is mathematically equal but a different
representation, and it reaches every caller of the public `simplify_gamma()` —
i.e. every UFO `ProjM`/`ProjP`, which is every chiral electroweak and Yukawa
vertex. The results below are unchanged; downstream consumers that pattern-match
on `projp` are not.

## Why it is correct

`Tr(ℙ₊ γ^μ γ^ν γ^ρ γ^σ)` reduces to

```
2 ( g^{μν} g^{ρσ} − g^{μρ} g^{νσ} + g^{μσ} g^{νρ} )  −  2 ε^{μνρσ}
```

i.e. exactly `½ Tr(γ^μγ^νγ^ργ^σ) + ½ Tr(γ5 γ^μγ^νγ^ργ^σ)` (idenso's `epsilon4`
convention folds the `i`). For the parity-even `gg→h` scalar coupling the `γ5`
(ε) piece cancels, leaving a clean parity-even numerator — which is what we see.

## Validation

| check | result |
|-------|--------|
| idenso unit test `projector_closed_loop_reduces_through_trace` | `Tr(ℙ₊ γμγνγργσ)` → exact metrics + ε |
| 21 FORM / FeynCalc trace/chain/γ5/γ0 reference tests | pass (no regression) |
| `gammalooprs` + full `gammaloop` binary | compile clean |
| **real `gg→h`, end-to-end** (`generate amp g g > h [ {{1}} ]` → `save dot --reduce`) | reduces to `C0(…; MT²,MT²,MT²) + 3·B0`, **zero residual `Tr`/`ℙ`/`chain`**; the `½ (…)` projector-expansion factors are visible in the coefficients, and there is no ε term (parity-even) |

> Note on running the idenso suite: the full `cargo test -p idenso --lib`
> SIGABRTs after ~13 tests — a **pre-existing** Symbolica one-instance-per-process
> limit (clean `main` aborts at the same test, unrelated to this change). Verify
> individual tests with `cargo test <exact> -- --exact`.

## Files touched

- `crates/idenso/src/dirac/simplify.rs` — the three additions above.
- `crates/idenso/src/dirac/test/mod.rs` — the unit test + `bis!` import.
- `crates/idenso/src/dirac/test/chirality_projector.rs` — regression tests for
  the 4D guard, the `evaluate_traces` gate, idempotency, the lone-projector
  (Yukawa `ū ℙ± v`) case, and projector-position epsilon sign.

See [the app integration](05-app.md) for how the reduced numerator reaches the
browser, and the app-coverage note for what this un-blocks.
