# oneloop benchmarks & validation harnesses

A practical how-to-run guide for the validation and benchmarking harnesses that
sit alongside the `oneloop` crate. The crate itself is the *reduction* step: it
turns a one-loop integral with an arbitrary polynomial numerator into a linear
combination of the four scalar masters A0/B0/C0/D0 with coefficients rational in
`d = 4 - 2ε` (see [the overview](../docs/01-overview.md) and
[the reduction algorithm](../docs/02-reduction.md)). These harnesses put numbers
on those reductions and check them against independent engines.

`benchmarks/` is committed for reference and reproducibility, but it is **not
part of the shipped library** — only `src/` merges upstream. The Rust harnesses
live in `benchmarks/rust/` and are wired as Cargo examples via `[[example]]`
entries in `Cargo.toml`; the Python scripts are validation and oracle drivers
that run directly against external one-loop libraries. Together they are the
evidence trail behind the claims in
[the benchmarks writeup](../docs/06-benchmarks.md) and
[the frontier writeup](../docs/04-frontier.md).

---

## Layout

```
benchmarks/
├── README.md                 # this file
├── madloop_reference.md      # MadLoop/MG5 reproduction record (validated numbers)
├── MONDAY_AGENDA.md          # status + speed tables condensed for a review call
├── rust/                     # Cargo example harnesses (cargo run -p oneloop --example NAME)
│   ├── bench.rs
│   ├── golden_master.rs
│   ├── emit_reductions.rs
│   ├── time_reduce.rs
│   ├── frontier_map.rs
│   ├── massless_legs.rs
│   ├── symbolic_delta.rs
│   └── reduce_regularized_draft.rs
└── python/                   # cross-engine + oracle validation (python3 SCRIPT)
    ├── crosscheck.py
    ├── verify_pentagon_reduction.py
    ├── verify_dotted_pentagon.py
    ├── verify_pentagon_rsp.py
    ├── verify_heptagon_degenerate.py
    ├── _tensor_oracle.py
    ├── _probe_conventions.py
    └── _probe_fjt.py
```

---

## Rust harnesses (`rust/`)

Each is a Cargo example over the public reducer API (`reduce()`,
`IntegralFamily`, `MasterIntegral`, …). Run any of them with:

```bash
cargo run --release --example NAME -p oneloop
```

Use `--release` for the timing harnesses so the numbers are meaningful. All of
them activate the Symbolica license on startup (via `gammalooprs`), and note the
one-instance-per-process constraint: a harness that deliberately triggers a
`reduce()` panic must run each degenerate case in its own process, because a
caught panic poisons global Symbolica state (see `massless_legs.rs` below).

| Example | What it does |
|---------|--------------|
| `bench` | Stress + timing sweep of the full reducer over every topology (tadpole…heptagon) at ranks 0–4, scalar and dotted, on physical 4D-offset kinematics (deliberately degenerate-Gram for N≥6). For each reduction it catches panics, asserts every coefficient is finite, checks the masters are valid, and times it. |
| `golden_master` | Characterization test. Reduces a battery of families and prints every `(master, coefficient)` evaluated at one fixed numeric point, so two algebraically different but equal reductions produce the *same* dump. Capture on known-good code, then require a byte-identical dump after a behaviour-preserving refactor (`… > /tmp/golden.txt` then `… \| diff /tmp/golden.txt -`). |
| `emit_reductions` | Emits machine-readable reductions (`coeff(d) * master(numeric args)`) for a battery of scalar + dotted families over two fixed spacelike massive geometries, for consumption by `python/crosscheck.py`. Redirect to a file: `… > /tmp/oneloop_reductions.txt`. |
| `time_reduce` | Per-topology absolute timing of `reduce()` for triangle/box/pentagon, scalar and tensor. Source of the sub-millisecond reduce-speed numbers reported in the reference records. |
| `frontier_map` | Maps the on-shell-massless-leg frontier: for a massive-internal triangle/box, varies the number of on-shell (zero-invariant) legs and the numerator rank, catching panics so one run prints the full OK/PANIC map. Backs the "triangles break at ≥2 on-shell legs, rank≥2; boxes are robust" finding in [the frontier writeup](../docs/04-frontier.md). |
| `massless_legs` | Deep dive on the off-shell-δ regularization fix as the number of massless legs grows (and for massless *internal* lines). Each config runs in its **own** process (config index passed as an argument) so one panic can't poison Symbolica state for the rest: `for i in $(seq 0 N); do cargo run --release --example massless_legs -p oneloop -- $i; done`. Prints a parseable `RESULT`/`TERM` reduction, or `PANIC`. |
| `symbolic_delta` | Blast-radius probe for the fix: checks that `reduce()` already accepts a *symbolic* invariant `delta` (on-shell leg → `delta`), so the coefficients come out rational in `delta` and the limit can be taken downstream — i.e. the fix needs no change to the core reducer. |
| `reduce_regularized_draft` | Prototype of the regularization wrapper, implemented entirely *outside* the core (it calls the unchanged public `reduce()`): replace the degeneracy-causing zero invariants with a symbol `delta`, reduce, then substitute `delta → 0` in the coefficients and master arguments. Demonstrates the on-shell-massless path before any `src/` integration. |

---

## Python harnesses (`python/`)

These are validation and oracle scripts, not part of the crate. They compare the
reducer's output (and the underlying identities) against **independent** one-loop
engines. Run them directly:

```bash
python3 SCRIPT.py [args]
```

### Prerequisites

Set up a Python virtualenv with:

- **OneLOopBridge** — Python bindings to OneLOop / `avh_olo` (van Hameren),
  imported as `oneloop_bridge`, used for numeric A0/B0/C0/D0 master values
  (Ellis–Zanderighi normalization; this is the same scalar-master library MadLoop
  links).
- **feynalg** — Denner-style analytic closed forms
  (`analytic_A0/B0/C0/D0`), used as a second, independent master engine. Scripts
  degrade gracefully if it is not importable.
- **scipy** / **numpy** — direct Feynman-parameter (Dirichlet / Monte-Carlo)
  integration of the *original* integral, and general numerics.
- **sympy** — symbolic manipulation in `crosscheck.py` / `_tensor_oracle.py`.

| Script | What it checks |
|--------|----------------|
| `crosscheck.py` | The cross-engine driver. Reads the reductions from `emit_reductions` and, for each family, forms the full Laurent series `V(ε) = Σ_i c_i(d=4−2ε)·M_i(ε)` with masters from OneLOop, then requires (a) the 1/ε and 1/ε² poles cancel and (b) the finite part matches a direct scipy integration. A running OneLOop-vs-feynalg tally on the masters is a third independent check. Run: `python3 crosscheck.py /tmp/oneloop_reductions.txt`. |
| `verify_pentagon_reduction.py` | Independent numeric check of the van Neerven–Vermaseren / Melrose / FJT reduction of scalar N-point integrals (N≥5) to (N−1)-point ones (pentagon→boxes, hexagon→pentagons, …) using the bordered modified-Cayley coefficients, at d=4 via Feynman-parameter integration. |
| `verify_dotted_pentagon.py` | Verifies the FJT/Tarasov index-lowering recurrence generalizes verbatim from the box (N=4) to the pentagon (N=5): checks LHS vs RHS at generic `d` by direct Feynman-parameter integration of every (dotted, and pinched) term. |
| `verify_pentagon_rsp.py` | Derives and verifies the reducible-scalar-product substitution rules for a pentagon numerator (`l.l`, `l.q_i`) as exact algebraic identities in the loop momentum, checked to machine precision at random Euclidean kinematics; confirms the pentagon has no ISP at the top topology. |
| `verify_heptagon_degenerate.py` | Verifies the degenerate-Cayley reduction used for N≥7, where `det(Y)=0` identically in 4D: builds a null vector of the bordered Cayley and checks the resulting partial-fractioning integrand identity to machine precision — validating the exact coefficients the Rust `degenerate_coeffs` computes. |
| `_tensor_oracle.py` | An independent, convention-free tensor one-loop integrator for numerators `(k²)^p`: analytic loop-momentum integration via symmetric-moment rules, leaving a numeric Feynman-parameter integral. A standalone oracle for the dotted/tensor cross-checks. |
| `_probe_conventions.py` | Empirical convention probe: computes A0/B0/C0/D0 at shared kinematic points with three engines (OneLOop, feynalg, scipy MC) side by side, to pin down the exact normalization and sign relations before trusting any cross-check. |
| `_probe_fjt.py` | Tests the FJT index-lowering identity for the dotted triangle `[2,1,1]` independently of the crate: evaluates both sides of the derived IBP relation with exact Feynman-parameter quadrature, isolating "is the formula right" from "is the Rust execution right". |

---

## Reference records

Two Markdown records capture the validated numbers behind the docs. Prose in
[the benchmarks writeup](../docs/06-benchmarks.md) and
[the frontier writeup](../docs/04-frontier.md) should trace back to these.

### `madloop_reference.md`

The MadLoop / MG5_aMC reproduction record. Anchored on Valentin's benchmark
`e+ e- > a > d d~ [virt=QCD]` (MG5_aMC v3.7.2), reproduced to ~14 digits. Contains,
among others:

- The reference numbers at the default PS point (μ_R = M_Z = 91.188 GeV,
  s = 10⁶ GeV²): Born = `3.4754514769164148e-03`; virtual normalized by
  `Born·α_s/(2π)`: Finite = `-8.9363792407373648e+00`,
  Single pole = `8.7724371707012754e+00`,
  Double pole = `-2.6666666666666670e+00` (= −8/3 = −2·C_F, matching exactly).
- The scope boundary: oneloop does the *reduction* of the massless on-shell
  triangle to C0(0,0,s;0,0,0) + bubbles + tadpoles; it does not build the spinor
  trace, evaluate the IR-divergent masters, or do the Born interference.
- Target-A results at `det(Y)=0`: scalar → C0(0,0,100;0,0,0); k² → B0(100;0,0);
  `dot(k,q1)` → ½B0(0;0,0) − ½B0(100;0,0) (hand-verified against OneLOop), with
  the tensor cases that hit the degenerate-Gram wall.
- The on-shell-massless frontier map (triangle tolerates 1 on-shell leg; boxes
  robust; the ≥2-on-shell-leg triangle at rank≥2 breaks) and the off-shell-δ
  **fix**: gg→h rank-2 `(k.q1)²` → (1/20)[B0(0;1,1) − B0(2/5;1,1)] = `-3.4748e-3`,
  matching an independent Feynman-parameter value to ~3e-10.
- The honest speed framing: oneloop's per-topology symbolic reduction is
  sub-millisecond (triangle rank-2 ≈ 0.27 ms), while MadLoop's warmed-up
  per-point *numerical* eval is 13.6 µs — a different operation. The value
  proposition is analytic closed-form reductions, not raw numerical speed.

### `MONDAY_AGENDA.md`

A condensed status record for a review call: the 20 reproduced
MadGraph/MadLoop processes (triangle, box, pentagon, 4-gluon, loop-induced),
the frontier + fix headline, and the per-operation speed table (symbolic
reduction vs per-point re-reduce vs the projected reduce-once + eval-masters mode).

---

## See also

- [`../docs/01-overview.md`](../docs/01-overview.md) — what the crate is and does.
- [`../docs/02-reduction.md`](../docs/02-reduction.md) — the reduction algorithm.
- [`../docs/03-numerators.md`](../docs/03-numerators.md) — numerator handling.
- [`../docs/04-frontier.md`](../docs/04-frontier.md) — the on-shell-massless frontier and the δ fix.
- [`../docs/05-app.md`](../docs/05-app.md) — the app / graph-bridge path.
- [`../docs/06-benchmarks.md`](../docs/06-benchmarks.md) — the validation story in prose.
- [`../docs/CHANGELOG.md`](../docs/CHANGELOG.md) — change history.
