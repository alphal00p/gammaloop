# oneloop — summary & document map

*One-page orientation — start here. Last updated 2026-08-21.*

## What it is (two sentences)
`oneloop` is a symbolic **one-loop IBP reducer**: it turns any one-loop Feynman
integral with an arbitrary polynomial (tensor) numerator into a linear combination
of the four scalar master integrals **A0/B0/C0/D0**, with coefficients **rational
in `d = 4 − 2ε`**. It is the **reduction** step only — *evaluating* the masters
(putting numbers on A0–D0) is delegated to external libraries (OneLOop/avh_olo, the
same one MadLoop links, and feynalg).

## Where it's at (one line)
**Complete and cross-validated for its scope** — any N-point, any tensor rank,
deployed live in the FeynmanEngine app; the remaining work is *different in kind*
(full-amplitude assembly, independent Mathematica checks, 2-loop), not more one-loop
validation.

## Results at a glance
| What | Result |
|---|---|
| Reducer completeness | any N-point, any tensor rank, raised powers — no guarded panic left in the tensor path |
| **Cross-engine correctness** | **132 / 132** integral families vs OneLOop + feynalg + scipy (+ 2 tensor oracles); 944/944 masters |
| High-rank + singular-Gram | rank-3–6 mixed tensors and the heptagon (N=7) validated (0.0–1.6σ) |
| Deployed app | reduces 5+ point diagrams live; **848** physical diagrams, **0** degenerate-Cayley walls |
| MadLoop reproduction | **~110** MG5_aMC processes run (96 fresh), self-consistent to ~14-digit accuracy; anchors match published values |
| Speed | sub-millisecond symbolic reduction (triangle → pentagon, scalar → rank-2) |
| Unit tests | 38 pass, 1 ignored, 0 fail |

## Two things to be precise about (for Valentin)
1. **What "validated" means for the reducer.** The *direct* check is the
   **132/132 cross-engine** result: `oneloop`'s reductions agree with OneLOop
   (MadLoop's own scalar-master library), feynalg, and a direct scipy
   Feynman-parameter integration to ~13-digit / machine precision. That is the
   rigorous validation of the tool's output.
2. **What the MadLoop batch is.** `oneloop` does the *reduction*, **not** the full
   amplitude, so the ~110 MG5_aMC runs are a **reproduction record** — MadLoop runs
   cleanly on our setup, the published anchors (`e⁺e⁻→dd̄ = −8.936`, `gg→gg = −66.63`,
   the loop-induced values) match to 14 digits, and the IR / colour-factor pole
   structure comes out textbook — mapping the *process space* the reductions cover.
   A per-process `oneloop`-vs-MadLoop *full-amplitude* diff would need the amplitude
   assembly (spinor trace + colour + Born interference) = "target B", a separate
   project.

## Document map — where to find what
| For… | Read |
|---|---|
| what it is, scope, status | [01-overview](01-overview.md) |
| **the reduction math** — masters, per-topology IBP recursions, N>4 bordered-Cayley | [02-reduction](02-reduction.md) |
| **tensor / dotted numerators** — RSP rules, PV transverse average, the gammaloop→family bridge | [03-numerators](03-numerators.md) |
| **the on-shell-massless frontier** — why massless legs are hard, the off-shell-δ fix, literature | [04-frontier](04-frontier.md) |
| the live app — the "Reduce to masters" button, end-to-end | [05-app](05-app.md) |
| **chirality-projector traces** — how closed fermion loops (`gg→h`, electroweak) now collapse instead of leaving an inert `Tr(…ℙ…)` (idenso fix) | [08-projector-traces](08-projector-traces.md) |
| the validation story in prose | [06-benchmarks](06-benchmarks.md) |
| **the results + numbers** — 132/132, the three MadLoop batches, tables | [07-benchmark-report](07-benchmark-report.md) |
| the MadLoop reference numbers (~110 processes, 3 batches) | [../benchmarks/madloop_reference.md](../benchmarks/madloop_reference.md) |
| how to re-run every benchmark | [../benchmarks/README](../benchmarks/README.md) |

## What's next / what's not done (honest)
- **Full amplitudes ("target B").** Assemble the reduced integrand into a full
  amplitude number (spinor trace + colour + Born interference) to compare `oneloop`
  end-to-end against MadLoop per process. A real project, best scoped with Valentin.
- **Independent Mathematica cross-check.** An expression-level check of the
  reductions — the natural next validation, done separately.
- **The threshold frontier.** Genuinely singular / threshold kinematics (different
  invariants vanishing at *different rates*) need a systematic Denner–Dittmaier
  expansion; the current single-δ regularization covers the tested on-shell-massless
  cases but not those.
- **2-loop.** Out of scope — a separate, long-horizon effort (no Laporta engine here).

---
*Code = `src/` (the mergeable library). Everything else is documentation plus
reference / gitignored benchmark harnesses.*
