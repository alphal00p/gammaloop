# oneloop reducer — Monday agenda (Valentin, 16h00 CEST / 10h00 ET)

## 1. Where it stands
One-loop IBP reducer: any N, tensor numerators, raised powers → A0/B0/C0/D0, symbolic (rational in d).
- **Validated generically** vs OneLOop + feynalg + scipy + 2 tensor oracles: scalar N=3–7, tensor to
  rank-6, UV-div, massless lines, timelike/threshold/near-degenerate. 27 unit tests, clippy clean.
- **20 MadGraph/MadLoop processes reproduced** (below). PR #86 open.

## 2. MadLoop cross-check — 20 processes reproduced (matches to ~14 digits where checked)
| topology | processes |
|---|---|
| triangle | e+e−→dd̄ (−8.936), Drell-Yan uū→e+e− (=−8.936, crossing), uū→Z, ud̄→W+, tt̄→g, e+e−→bb̄ |
| box | uū→dd̄, dd̄→dd̄, uū→gg=dd̄→gg=gg→dd̄ (−54.03, crossing), gu→gu, uū→tt̄ (massive top) |
| pentagon | e+e−→dd̄g |
| 4-gluon | gg→gg (123 loop diagrams, −66.63) |
| loop-induced | gg→h, gg→hh, **H→γγ (28 W+top loops, 6.64e-2)**, γγ→γγ (186 loops, 1.05e-3) |

## 3. The frontier + THE FIX (headline)
Reducer is a **massive** method (mass-derivative/Tarasov, divides by Gram/Cayley dets). It's complete
off-shell; **on-shell massless legs** (universal: gluons, photons, light quarks) make pinched
sub-Grams/Cayleys vanish → it panics. Frontier is narrow: triangle tolerates 1 on-shell leg; **boxes
are robust**; scalar+rank-1 fine on-shell; only the ≥2-on-shell-leg **triangle at rank≥2** breaks.

**FIX (Eli's idea, validated to 1e-10):** regularize with off-shellness δ, reduce, take δ→0. The
reducer's exact rational arithmetic **cancels the 1/δ inverse-Gram poles automatically** → coeffs
finite at δ=0. gg→h rank-2 → (1/20)[B0(0;1,1)−B0(2/5;1,1)] = −3.47483e-3 vs independent Feynman-param
−3.47483e-3 (diff 3e-10). **Blast radius on the reducer = ZERO**: reduce() already accepts symbolic
δ; the fix is a ~50-line wrapper (working prototype: `reduce_regularized_draft.rs`). Connects to
Valentin's dimension-shift / vacant work — ask his read on the analytic-δ limit + multi-massless case.

## 4. Speed — "why slower than MadLoop?" (the honest answer)
We measured a **different operation**. oneloop does SYMBOLIC exact-rational algebra in d (makes a
*formula*); MadLoop does optimized FLOAT numerical eval (makes a *number* at one point). Symbolic is
inherently ~10–100× slower — that's the price of analyticity, not a defect.

| operation (triangle, rank-2) | oneloop | MadLoop | ratio |
|---|---|---|---|
| **symbolic reduction** (once → analytic coeffs) | 0.27 ms | — (float codegen at build) | different task |
| **per-point, re-reducing each point** (current use) | 0.27 ms | 0.0136 ms | **~20× slower (+1900%)** |
| **per-point, reduce-once + eval masters** (future mode) | ~0.003 ms (proj.) | 0.0136 ms | **~4× faster (proj.)** |

- Absolute reduce speed (generic massive): tri scalar 0.016 / rank-1 0.097 / rank-2 0.27 ms; box
  0.029 / 0.14 / 0.54; pentagon 0.50 / 0.73. **Sub-millisecond across the board.**
- Per-point floor is the scalar masters (avh_olo: C0 ~0.5 µs), which **MadLoop links too** — so a
  "reduce-once-symbolically-in-the-invariants, evaluate-per-point" mode would be master-dominated
  (~µs) and competitive. Report absolute ("0.02–0.7 ms/reduction, analytic output"); let Valentin
  place it in the 1-loop landscape and reason about 2-loop scaling.

## 5. Asks / open questions for Valentin
1. Read on the off-shell-δ fix: analytic limit (symbolic δ → expand → δ→0) vs numerical; is the
   single-shared-δ valid for multiple massless legs? (his dimension-shift expertise).
2. Path to real amplitudes: the app/graph bridge (spenso-contracted numerator → oneloop dot form).
3. Where does this sit / timeline / most useful next direction?
4. (later) 2-loop vacuum via parametric IBP he floated — is that the natural growth path?
