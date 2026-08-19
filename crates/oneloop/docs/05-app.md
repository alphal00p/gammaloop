# The FeynmanEngine app integration

This document describes the **end-to-end app path**: how a Feynman diagram drawn
in the browser becomes a one-loop numerator reduced to the four scalar master
integrals (A0/B0/C0/D0), rendered as typeset math. It is the "application"
counterpart to the reducer internals — see [the overview](01-overview.md),
[the reduction algorithm](02-reduction.md), [numerators](03-numerators.md), and
the [known frontier](04-frontier.md) for what happens inside the box.

The integration is live at **<https://gammaloop.hirschi.lu/feynmangraph/>**.

---

## The pipeline at a glance

```
Browser (React/ReactFlow canvas)
  │  draw diagram, click "Reduce to masters"
  ▼
FastAPI  POST /api/reduce            (feyngraph/api/reduce.py)
  │  diagram → gammaloop .dot, run gammaloop CLI
  ▼
gammaloop  save dot --output-full-numerator --reduce
  │  build contracted numerator, then:
  ▼
reduce_graph_numerator  (crates/gammalooprs/src/reduce_bridge.rs)
  │  Dirac trace → scalar  ─►  Q→loop basis (K/P)  ─►  IntegralFamily
  ▼
oneloop::reduce::reduce  →  Σ coeff × master
  │  emitted as dot attributes reduced_num / reduce_status
  ▼
FastAPI parses the .dot, returns { raw, warnings, reason }
  ▼
Frontend sanitizes for Typst, renders  Σ (coeff) B_0(...) + …
```

Each stage is described below.

---

## 1. The glue: `reduce_graph_numerator`

The bridge from a gammaloop `Graph` to the reducer lives in
`crates/gammalooprs/src/reduce_bridge.rs`. Its entry point,
`reduce_graph_numerator(graph, num)`, returns a **typed outcome** rather than a
string-or-error, so callers can distinguish "reduced" from the several expected
"can't reduce" cases:

```rust
pub(crate) enum ReduceOutcome {
    Reduced(String),   // Σ coeff × master, formatted for the dot `reduced_num`
    NotOneLoop(usize), // tree / multi-loop; carries the loop count
    ZeroNumerator,     // one-loop, but the contracted numerator is identically 0
    Unsupported,       // one-loop, non-zero numerator, but reducer produced no terms
}
```

The function performs the full symbolic path:

1. **Loop-count gate.** It reads `graph.loop_momentum_basis.loop_edges.len()`.
   If it is not exactly `1`, it short-circuits to `NotOneLoop(n)` — the reducer
   only handles a single loop.

2. **Complete the Dirac traces to a scalar.** The incoming numerator still
   carries open γ-chains. It is collapsed with
   `num.collect_gamma_chains().simplify_gamma().expand().simplify_metrics()`,
   producing a scalar expression in dot products and the loop momentum (see
   [numerators](03-numerators.md) for the algebra the reducer expects).

3. **Rewrite edge momenta into the loop basis.** Physical edge momenta
   `Q(edge, …)` are replaced by the loop-momentum basis `K(0, …) + P(j, …)` via
   `graph.integrand_replacement(...)`, then re-expanded and metric-simplified.
   This puts every propagator and every numerator dot product in terms of the
   single loop momentum `K` and external momenta `P`.

4. **Collect loop propagators with masses.** For each edge it forms the
   loop-basis representative (`loop_atom` + `ext_atom`) and the squared mass
   `mass * mass`, skipping edges whose loop part is zero (external / tree edges).
   The result is a `Vec<GammaloopEdge>` of `{ lmb_rep, mass_sq }`.

5. **Bridge to an `IntegralFamily` and reduce.** With a `GammaloopHeads`
   descriptor (which symbols mean loop momentum, external momentum, Lorentz
   index, metric), it calls `oneloop::bridge::family_from_gammaloop(...)` to
   build the family, then `oneloop::reduce::reduce(&family)`.

6. **Classify the result.** If `reduction.terms` is empty, the numerator was
   either identically zero (`ZeroNumerator`) or non-zero but not reducible by
   the current engine (`Unsupported`). Otherwise each `(coeff, master)` pair is
   formatted `(coeff) MASTER(args…)` and joined with `+`.

### Master formatting

Masters are printed by `format_master` in the canonical A0/B0/C0/D0 shapes the
frontend and downstream tooling expect:

| Master   | Rendered form                                             |
|----------|-----------------------------------------------------------|
| Tadpole  | `A0(m²)`                                                  |
| Bubble   | `B0(p²; m1², m2²)`                                         |
| Triangle | `C0(p1², p2², p12²; m1², m2², m3²)`                        |
| Box      | `D0(p1², p2², p3², p4², s, t; m1², m2², m3², m4²)`         |

Coefficients and invariants are printed through Symbolica's Typst printer
(`SpensoPrintSettings::typst().typst_symbolica()`), and the whole string has its
inner quotes escaped so it survives embedding in a `.dot` attribute.

---

## 2. The CLI surface: `save dot --output-full-numerator --reduce`

The reducer is exposed through gammaloop's dot exporter. The relevant setting is
`DotExportSettings::reduce` (in `crates/gammalooprs/src/processes/mod.rs`),
documented as: *emit a `reduced_num` attribute — the full numerator reduced to
master integrals via the `oneloop` IBP reducer.* Turning it on **forces gamma
and color algebra on** (the reducer needs a contracted numerator), independent
of the `--do-gamma-algebra` / `--do-color-algebra` flags.

During serialization (`graph/parse/serialization.rs`) the export matches on the
`ReduceOutcome` and writes graph-level dot attributes:

| Outcome              | Dot attributes emitted                          |
|----------------------|-------------------------------------------------|
| `Reduced(s)`         | `reduced_num = "…"`                             |
| `NotOneLoop(n)`      | `reduce_status = "not_one_loop"`, `reduce_loops = n` |
| `ZeroNumerator`      | `reduce_status = "zero_numerator"`             |
| `Unsupported`        | `reduce_status = "unsupported"`                |

The unreduced contracted numerator is always emitted alongside as `full_num`, so
a diagram that fails to reduce still shows its (Typst-printed) numerator.

---

## 3. The FastAPI endpoint: `POST /api/reduce`

`feyngraph/api/reduce.py` turns a canvas diagram into a reduced expression by
driving the CLI. In outline:

1. Locate the gammaloop binary; `503 GAMMALOOP_NOT_FOUND` if missing.
2. Load the model and convert the posted `GraphSpec` to a gammaloop `.dot`
   (`to_gammaloop_dot`), returning `422` for unassigned edges or missing
   external legs.
3. Write a temporary run-card TOML whose command block is:
   ```
   import model <model>.json
   import graphs graph.dot -p amp
   save dot --output-full-numerator --reduce
   ```
   and run `gammaloop <toml> run g` in a temp directory.
4. Read the first emitted `.dot` from
   `state/processes/amplitudes/amp/default/` and scrape it with two regexes:
   `reduced_num = "…"` and `reduce_status = "…"`.

The response model is:

```python
class ReduceResponse(BaseModel):
    raw: str
    format: str = "typst-symbolica"
    warnings: list[str] = []
    reason: str | None = None
```

Crucially, a diagram that **does not** reduce is not treated as an error. If
there is no `reduced_num` but there is a known `reduce_status`, the endpoint
returns `200` with `raw=""`, `reason=<status>`, and a friendly warning string.
The status→message map mirrors the reducer's typed outcomes:

| `reduce_status`   | User-facing message                                                                          |
|-------------------|-----------------------------------------------------------------------------------------------|
| `not_one_loop`    | "Reduce to masters only works for one-loop diagrams."                                          |
| `zero_numerator`  | "This diagram vanishes identically — its numerator is zero, so there is nothing to reduce."    |
| `unsupported`     | "This one-loop diagram isn't supported by the reducer yet."                                    |

Only a genuinely unexpected shape — a diagram emitted with no `reduced_num` and
no recognized `reduce_status` — is surfaced as an error (`422 NO_REDUCTION`).
If gammaloop emits more than one diagram, the endpoint reduces the first and
adds a warning saying so.

---

## 4. The frontend: "Reduce to masters"

The reduce action is guarded and post-processed on the client, in
`frontend/src/panels/reduceGuard.ts`.

### Loop-count gate → yellow warning

The client first fetches the diagram's loop count from the backend
(`POST /api/validate-graph`) and, if it is not exactly one, warns without
attempting a reduction. `reduceLoopGuard(loopCount)` returns `null` for a
one-loop diagram (proceed), and otherwise a message like:

> Reduce to masters only works for one-loop diagrams — this diagram has *N* loops.

This is shown as a non-blocking yellow warning rather than letting the user
trigger a backend error.

### Friendly backend reasons

If the request does come back with a `reason` (zero numerator, unsupported,
etc.), `reduceReasonMessage(status)` maps it to the same friendly text the
backend uses, keeping the two in sync. An unknown status falls through to the
normal error path.

### Typst sanitization

The reducer prints math with conventions that Typst's math mode misreads as
undefined variables or stray function calls: `B0(...)`, `dot(a,b)`, flat momenta
like `q1`, and bare identifiers like `Tr` and `ZERO`. `sanitizeReducedTypst(raw)`
rewrites them into valid Typst math:

| Reducer output | Typst math      | Rule                                             |
|----------------|-----------------|--------------------------------------------------|
| `dot(a, b)`    | `(a dot b)`     | function-call dot → infix `dot` operator         |
| `B0(`          | `B_0 (`         | master name → subscripted symbol (A/B/C/D)       |
| `q1`           | `q_(1)`         | flat momentum → subscript                        |
| `Tr`, `ZERO`   | `"Tr"`, `"ZERO"`| quote residual multi-letter identifiers          |

The final rule quotes any remaining bare multi-letter identifier that isn't
already quoted and isn't the `dot` operator; it uses lookarounds (not `\b`) so a
trailing Unicode superscript such as `ZERO²` still matches. Per the source, this
sanitizer was validated with the Typst compiler across **43 diagrams**
(bubble/triangle/box, open and closed fermion lines).

---

## 5. Deployment

The app is deployed to Valentin's EC2 server and reverse-proxied by nginx under
`location /feynmangraph/` → `127.0.0.1:8000` (the trailing slash strips the
prefix, so the app runs at root internally; the subpath is injected only at
build time via `FEYNGRAPH_FRONTEND_BASE=/feynmangraph/`).

- **gammaloop binary.** The server builds gammaloop from the `oneloop` branch
  (the branch that carries `reduce_bridge.rs` and the `--reduce` CLI surface).
  A wrapper picks the `target/release` binary, and the FastAPI service points at
  it via `FEYNGRAPH_GAMMALOOP_BIN`.
- **Redeploy.** `scripts/deploy-server.sh` (run on the server) does the app-side
  redeploy: `git pull --ff-only` → rebuild the frontend with the subpath base →
  swap it into `feyngraph/data/frontend` → restart the systemd service → health
  check. Because the Python side is an editable install, `git pull` updates it
  automatically; only the frontend is rebuilt. The deploy script does **not**
  run `npm install`, so a new frontend dependency must be installed on the
  server manually first.
- **Health check.** `GET /feynmangraph/api/health` returns `{"status":"ok"}`.

---

## What works, and where the edges are

- **Works end-to-end.** Draw a one-loop diagram, hit "Reduce to masters", and
  get `Σ (coeff) × master(args)` typeset in the browser — this is the
  demonstrated path for scalar and dotted bubble/triangle/box numerators.
- **Honest failure modes.** Multi-loop and tree diagrams are gated (yellow
  warning); vanishing numerators and reducer-unsupported numerators return a
  clear reason, not a stack trace.
- **Frontier.** Numerators the current engine cannot yet reduce surface as
  `unsupported`; the reach of the reducer itself — including the tensor /
  degenerate-Gram cases at the IR boundary — is documented in
  [the frontier](04-frontier.md). For the validated reference numbers behind the
  "works" claims, see [benchmarks](06-benchmarks.md) and
  [../benchmarks/README.md](../benchmarks/README.md).
```