# MadLoop reference — e+ e- > a > d d~ [virt=QCD]  (Valentin's benchmark)

Reproduced with MG5_aMC v3.7.2 (python3.11, gfortran 15.2), 2026-07-30. Matches Valentin's
numbers to ~14 digits.  Install: ~/Documents/GitHub/MG5_aMC_v3_7_2 (NOT in git).
Rerun: `python3.11 bin/mg5_aMC /tmp/mg5_bench.txt` with
  import model loop_sm-default
  generate e+ e- > a > d d~ [virt=QCD]
  output standalone BENCHMARK_TEST
  launch -f

## Reference numbers (default PS point, mu_R = M_Z = 91.188 GeV, s = 10^6 GeV^2)
Born        =  3.4754514769164148e-03   (GeV^0)
virtual normalized by Born * alpha_s/(2*pi):
  Finite      = -8.9363792407373648e+00
  Single pole =  8.7724371707012754e+00
  Double pole = -2.6666666666666670e+00

## Physics content
- Diagrams: 1 Born, 1 loop, 1 R2, 0 UV.  The single loop = the QCD vertex correction to the
  photon-quark-quark vertex: a gluon exchanged between the d and d~ lines => a TRIANGLE with 3
  massless loop propagators (2 quark + 1 gluon).
- Kinematics: quarks on-shell massless (p3^2 = p4^2 = 0), photon off-shell q^2 = s = (p3+p4)^2
  = 10^6 GeV^2 (two 500 GeV beams head-on).  => the massless on-shell triangle C0(0,0,s;0,0,0),
  IR-divergent (1/eps^2, 1/eps, finite).

## Analytic cross-check of the pole structure (V/B in units of alpha_s/2pi, C_F = 4/3)
V/B ~ C_F (mu^2/s)^eps [ -2/eps^2 - 3/eps - 8 + pi^2 ] + ...
- Double pole: -2 C_F = -8/3 = -2.66667                          => matches MadLoop EXACTLY.
- Single pole: C_F[-3 - 2 ln(mu^2/s)] = (4/3)[-3 - 2 ln(0.0083)] ~ 8.773  => matches ~8.7724.
- ln(mu^2/s) = ln(8315.25/1e6) = -4.7899.

## What oneloop's part IS (scope)
oneloop does the REDUCTION only: given the vertex-correction loop integrand (tensor numerator from
the fermion trace + gluon), reduce it to C0(0,0,s;0,0,0) + bubbles + tadpoles with rational-in-d
coefficients.  It does NOT build the spinor-trace numerator (spenso/MadGraph), evaluate the
IR-divergent masters (OneLOop / analytic), or do the Born interference + color/spin sum (MadGraph).
=> Reducer-level validation (target A): reduce the massless on-shell triangle & check vs OneLOop's
IR-divergent masters (the one regime the OneLOop harness hadn't covered).  Full -8.936 end-to-end
(target B) needs the whole amplitude assembly = the app/graph-bridge integration (plan w/ Valentin).

## Subtlety to test (massless on-shell IR)
Pinched sub-topologies of the massless on-shell triangle are SCALELESS massless bubbles
B0(0;0,0) = 0 in dim reg (1/eps_UV - 1/eps_IR), and tadpoles A0(0) = 0.  The sharp test is whether
oneloop's numerator reduction + master evaluation handles these vanishing scaleless pieces correctly.

## Target-A results (2026-07-30): oneloop on the massless on-shell triangle (det(Y)=0)
Cayley Y = [[0,s01,s02],[s01,0,s12],[s02,s12,0]] with masses 0; onsh lex [s01,s02,s12]=[0,0,100]
=> Y has a zero row (two on-shell legs) => det(Y)=0 (DEGENERATE, the classic inverse-Gram problem).

WORKS (reduces correctly at det(Y)=0, verified by hand + OneLOop, mu^2=1):
- scalar   -> C0(0,0,100;0,0,0).  OneLOop eps^-2 = +0.01 = 1/s  (textbook massless triangle). OK
- k^2      -> B0(100;0,0)  exactly (the k^2=D0 identity, drop line 0).  No spurious scaleless terms. OK
- dot(k,q1)-> 1/2 B0(0;0,0) - 1/2 B0(100;0,0).  Hand-check: r1^2=0 => D1=k^2+2k.q1 => 2 dot(k,q1)=
  D1-D0 => 1/2[B0(s02;0,0)-B0(s12;0,0)] = 1/2 B0(0;0,0) - 1/2 B0(100;0,0).  MATCHES oneloop exactly.
  B0(0;0,0)=scaleless=0 (OneLOop returns 0) => physical finite = +1.3026 - 1.5708j = 1/2(2.605 - pi i). OK
  => the scaleless-bubble handling is CORRECT (right coeff, vanishes on evaluation).

WALL (panics "degenerate Cayley reduction: exceptional kinematics", reduce.rs:319):
- dot(k,q2) and raised-power [2,1,1].  dot(k,q2) is actually SCALELESS=0 (both bubble pieces have
  s=0: 2 dot(k,q2)=D2-D1 => B0{0,1}-B0{0,2}, both invariants 0), but the reducer routes the C0 piece
  through degenerate_coeffs() whose bordered-Cayley nullspace has no v[n]!=0 => panic instead of
  returning the (zero) value.  This is the known EXCEPTIONAL-KINEMATICS / inverse-Gram-determinant
  problem: the integral is well-defined (standard tools compute it via special IR reduction, cf.
  Denner-Dittmaier), oneloop's current algorithm doesn't cover it and cleanly refuses.

HONEST SUMMARY for Valentin: the reducer reproduces MadLoop's benchmark SETUP, understands the pole
structure analytically (double pole -2 C_F = -8/3, single ~8.77), and correctly reduces the massless
on-shell triangle's scalar + reducible/linear numerators at det(Y)=0.  The genuinely-tensor numerators
at this exceptional (on-shell massless) point hit the degenerate-Gram wall -> a clear next work item
(special IR/degenerate-kinematics reduction).  Full -8.936 end-to-end still needs the amplitude
assembly (spinor trace + color + Born interference) = the app/graph-bridge project.

## Speed comparison (2026-07-30, rough) + gg>hh
- oneloop reduce() on a rank-1 triangle numerator: **0.067 ms/reduction** (3 terms), 200-run avg.
  This is the SYMBOLIC reduction, done ONCE per topology+numerator structure -> coefficients rational
  in d. Per phase-space point afterward = evaluate those rationals (ns) + scalar masters (OneLOop, us).
- MadLoop ./check: **20.4 ms/invocation** — but this is STARTUP-DOMINATED (loads param cards, inits
  ninja/collier/cuttools, reads filters) + 1 loop-ME eval. NOT a fair per-point number; MadLoop's
  warmed-up per-point loop eval is a fraction of this. A rigorous per-point comparison (loop check_sa.f
  over N points after warmup) is a Monday item -- do NOT overclaim a speedup from the 20ms.
- HONEST framing: the real structural advantage is symbolic-reduce-ONCE + trivial-per-point vs
  MadLoop's numerical-reduce-EVERY-point; amortizes over a phase-space scan. Raw reduce() at 0.067ms
  shows the reduction itself is cheap.

## gg > h h / b b~ [virt=QCD]  (loop-induced, Valentin's 2nd example) -- reproduced
0 Born, 8 loops, 2 R2 (top-quark boxes+triangles; loop-induced Higgs pair production).
Loop |M|^2: Finite = 3.4123262140874510e-05, poles = 0 (finite, as expected). accuracy 5.2e-13.
=> MG5 handles the harder box-level process cleanly; a good future reducer target (4-point boxes).

## RIGOROUS per-point timing (2026-07-30, corrected) -- IMPORTANT, do not overclaim
Edited check_sa.f: warmed up (1 eval), then timed 2000 loop-ME evals with cpu_time (numbers still
match -8.9363792407). Result: MadLoop per-point loop-ME eval = **0.0136 ms (13.6 us)** per point.
- MadLoop 13.6 us/point (numerical) is FASTER than oneloop reduce() 67 us (symbolic, one-time).
- So the earlier "symbolic amortizes -> oneloop faster per point" hypothesis is WRONG. MadLoop is
  mature optimized compiled code; its per-point numerical eval is very fast (it did the expensive
  process-specific codegen once at `output` time, ~minutes).
- oneloop's per-point eval would ALSO call OneLOop for the scalar masters (the SAME avh_olo library
  MadLoop links) => oneloop per-point ~ same order as MadLoop, not dramatically faster.
- HONEST value proposition: oneloop is NOT about beating MadLoop on raw numerical speed. Its
  differentiator is producing ANALYTIC closed-form reductions (rational-in-d coefficients x masters):
  formulas not numbers, reusable across kinematics, human-readable, composable/symbolic. That is the
  story for Valentin -- analyticity, not speed.

## UNIFIED FINDING (2026-07-30): on-shell massless legs are the reducer's frontier
Benchmarked 5 processes with MadLoop (all reproduced): e+e->ddx (triangle), gg>hh (boxes,
loop-induced), gg>h (massive top triangle, loop-induced, finite 9.37e-3), uux>ddx (boxes,
finite -42.34), ddx>ddx (boxes, finite -24.38).

Reducer behaviour, pinned down definitively via the gg>h massive top triangle (m^2=1) at
on-shell vs off-shell external legs:
- OFF-SHELL legs (0.05,0.05,0.4): reduces EVERYTHING - scalar, rank-1, rank-2 (q1q1, q1q2),
  rank-3 (llq1), raised power [2,1,1]. All clean, 0 panics. (Also the whole generic-kinematics
  benchmark suite, 100+ families.) The triangle itself is non-degenerate (det Y = -0.32).
- ON-SHELL massless legs (p^2=0, i.e. s01=s02=0 -- gluons/photons/light quarks): scalar + rank-1
  numerators still reduce, but
    * rank-2 numerators PANIC "singular Gram matrix" (a pinched sub-bubble's external leg is the
      null on-shell leg => 1x1 sub-Gram = q^2 = 0), and
    * raised power [2,1,1] PANICS "degenerate Cayley" (pinched bubble B0(0;m,m) has det(2x2 Cayley)=0).
  Both guards fire CLEANLY (no garbage).

ROOT CAUSE (unified): on-shell massless external legs make PINCHED sub-topologies degenerate
(zero sub-Gram / sub-Cayley). This is the classic inverse-Gram-determinant / exceptional-kinematics
problem. Fix = special IR reduction / dimensional recurrence (Tarasov d->d+2, Denner-Dittmaier).

WHY IT MATTERS (the Monday headline): on-shell massless legs are UNIVERSAL in real collider
processes (external gluons, photons, light quarks). So despite full validation on generic
kinematics, THIS is the concrete gap between the reducer and real-process amplitudes -- precisely
what Valentin's MadLoop suggestion was meant to surface. Scalar + rank-1 already work on-shell;
rank>=2 and raised powers need the degenerate-kinematics reduction. Well-scoped next work item.

## FRONTIER MAP refined (2026-07-31) — it's more nuanced/encouraging than "on-shell breaks it"
frontier_map.rs (catch_unwind per case), massive internal m^2=1, rank-2 (k.q1)^2 unless noted:
TRIANGLE vs # on-shell (zero-invariant) legs:  0 -> OK(4); 1 -> OK(3); 2 -> PANIC(singular Gram);
  3 -> PANIC. => a triangle tolerates ONE on-shell leg; it takes TWO to break rank-2.
TRIANGLE at 2 on-shell legs, by rank:  scalar OK; rank-1 dot(k,q1) OK; rank-1 k^2 OK; rank-2 PANIC
  (singular Gram); raised power [2,1,1] PANIC (degenerate Cayley). => boundary is rank-1 -> rank-2.
BOX (massive internal): all off-shell -> OK(8); 2 on-shell legs -> **OK(7)**. => BOXES are ROBUST
  to on-shell legs at rank-2 (richer topology: pinches to triangles, not directly to null-leg
  bubbles). The degenerate frontier is WORST for the triangle (fewest legs).
NB "OK" = reduces WITHOUT PANIC (n terms); the triangle rank-1/k^2/scalar cases were also HAND-
  VERIFIED (see Target-A). The box on-shell numbers aren't numerically cross-checked yet (on-shell
  box masters are IR-divergent; needs OneLOop IR + assembly).

6th process reproduced: e+ e- > d d~ g [virt=QCD] (pentagon, 5 ext): 4 Born, 30(+4) loops, 18 R2,
52 UV; finite = 2.4369. Full process set now spans triangle, box, pentagon, loop-induced.

REFINED TAKEAWAY for Monday: the reducer is complete off-shell; on-shell massless legs degrade it
only in a SPECIFIC, characterized way -- triangles with >=2 on-shell legs at rank>=2 (or raised
power). One on-shell leg is fine; boxes are fine; scalar+rank-1 are fine on-shell. So the gap to
real processes is narrower than "on-shell breaks it" -- it's the multi-on-shell-leg TRIANGLE at
rank>=2, exactly the IR/exceptional-kinematics corner needing special reduction.

## THE FIX (2026-07-31): off-shell regularization of on-shell tensor integrals -- WORKS
Eli's insight: the reducer is a MASSIVE-configuration method (mass-derivative/Tarasov recursions
that divide by Gram/Cayley dets). On-shell massless legs make those dets vanish -> breaks. FIX =
regularize: reduce at small off-shellness p^2=delta (where the massive reducer works perfectly, and
its coefficients are EXACT rationals so the 1/delta inverse-Gram poles are explicit), then take
delta->0. The 1/delta poles CANCEL in sum_i c_i M_i.

Demonstrated on gg>h massive top triangle (m^2=1, legs (delta,delta,0.4)), rank-2 (k.q1)^2:
delta-ladder eps^0 (finite): 1e-1:+4.2e-4, 1e-2:-2.89e-3, 1e-3:-3.415e-3, 1e-4:-3.469e-3,
1e-5:-3.4742e-3. Fit A + B*delta*ln(delta) + C*delta => on-shell limit A = -3.4752e-3.
The UV pole eps^-1 = 0.25*delta -> 0 (physically correct: (k.q1)^2 UV div ~ q1^2). No numerical
breakdown to delta=1e-5 (exact-rational coeffs keep the cancellation clean).

INDEPENDENT VALIDATION (convention-free): on-shell the trace piece ~ q1^2 = 0, so the integral =
C0 * <(P.q1)^2>_F (F-weighted Feynman-param simplex average), C0 from OneLOop. = -3.474834e-3.
vs off-shell->0 extrapolation -3.475174e-3.  DIFFERENCE = 3.4e-7.  MATCH.

=> The reducer's LAST issue (on-shell massless legs, universal in real processes) is FIXABLE with NO
new degenerate-reduction algorithm, just off-shell regularization + limit. Production version:
ANALYTIC (keep the leg virtuality symbolic, reduce, expand coeffs in delta, take delta->0 -- avoids
even the numerical residual) or NUMERICAL (delta~1e-4 gives ~4 digits). This IS Eli's massive-method
design working as intended. Great Valentin discussion: connects to his dimension-shift/vacant work.

## FIX BLAST RADIUS (2026-07-31): ZERO core changes -- symbolic-delta already works
Probe (symbolic_delta.rs): reduce() ACCEPTS a SYMBOLIC invariant delta (no panic). With the on-shell
legs kept = delta, the coefficients come out POLYNOMIAL in delta (finite at delta=0): the reducer's
EXACT rational arithmetic cancels the 1/delta inverse-Gram poles automatically. gg>h rank-2 (k.q1)^2
-> 1/4 delta^2 * C0(delta,2/5,delta) + (..delta..) B0(delta) + (..) B0(2/5). At delta=0 this collapses
to (1/20)[B0(0;1,1) - B0(2/5;1,1)]; evaluated: eps^-1=0 (UV pole cancels), eps^0 = -3.4748337e-3, vs
validated -3.474834e-3, DIFF 3.2e-10. => the fix needs NO change to reduce()/reduce_cayley/gram_solve/
degenerate_coeffs. It's a thin WRAPPER: (1) replace degeneracy-causing zero invariants/masses with a
symbol delta, (2) call reduce() unchanged, (3) substitute delta->0 in coeffs + master args, (4) eval.
Also covers massless INTERNAL lines (regularize m^2=delta). Multiple massless legs: one shared delta.
Design Qs to settle: how to DETECT which zeros to regularize (pre-check dets vs catch_unwind), robust
delta->0 (substitute + assert no residual 1/delta = genuine divergence), same-delta-for-all validity.
14 processes batch-reproduced (u/dx>gg, gg>ddx, Drell-Yan uux>epem = -8.936 by crossing, epem>bbx
massive b, udx>wp, uux>z, gu>gu, uux>ttx massive top, gg>gg 4-gluon 123 loops -66.63, ...). ~18 total.

## SPEED TABLE (2026-07-31): oneloop reduce() per topology (absolute, generic massive kinematics)
  triangle scalar   0.016 ms | rank-1  0.097 ms | rank-2  0.272 ms
  box      scalar   0.029 ms | rank-1  0.142 ms | rank-2  0.543 ms
  pentagon scalar   0.504 ms | rank-1  0.727 ms
=> sub-millisecond symbolic reduction across triangle/box/pentagon, scalar->rank-2. This is the
ABSOLUTE reducer speed (produces analytic rational-in-d coefficients, ONCE per topology+numerator).
Context for Valentin: MadLoop's warmed-up numerical per-point eval was 13.6 us (faster, but numerical
& process-specific-codegen). Different value props: oneloop = analytic formulas, sub-ms; MadLoop =
optimized numbers. Report as absolute ("~0.02-0.7 ms/reduction"), let Valentin contextualize vs the
one-loop landscape and how it may scale to 2-loop.

## Loop-induced (2026-07-31) + reduce_regularized prototype WORKS + speed explanation
- h > a a [virt=QED] = H->gamma gamma MILESTONE: 0 Born, 28 loops (W+top), finite 6.6360e-2.
- a a > a a [virt=QED] light-by-light: 0 Born, 186(+30) loops, finite 1.0538e-3. (loop-induced needs
  the [virt=QED] bracket; bare "generate h > a a" -> NoDiagram). => 20 processes total.
- reduce_regularized_draft.rs (thin wrapper, ZERO core changes): on the gg>h on-shell rank-2 case
  (bare reduce PANICS "singular Gram") it returns 1/20 B0(0;1,1) - 1/20 B0(2/5;1,1) = -3.4748e-3.
  WORKS. Gotcha: a caught reduce() panic POISONS Symbolica global state, so don't catch_unwind a
  reduce then reduce again in the same process; detect degeneracy with a pre-check (is_zero), not
  catch. Uses s.is_zero() for detection, expr.replace(delta).with(Zero) for the limit (reduce.rs:530).
- SPEED explanation (why "slower" than MadLoop): we measured SYMBOLIC reduction (0.27ms tri rank-2,
  makes a formula in d); MadLoop's 13.6us is NUMERICAL float eval (makes a number). Symbolic ~10-100x
  slower inherently. Per-point re-reduce = ~20x slower (+1900%); but reduce-once + eval masters (proj
  ~3us, C0=0.48us x ~6) = ~4x FASTER, master-dominated (avh_olo, which MadLoop links too). Table in
  benchmarks/MONDAY_AGENDA.md.

## DEEP-DIVE: on-shell-massless fix across 2/3/all massless legs (2026-08-03, massless_legs.rs)
Tested 13 configs (each in its own process to avoid the catch_unwind poison). ZERO PANICS across all
— the fix (reduce() regularizes zero invariants -> delta -> 0) handles every case. Results (all
physically sensible, verified by identity/limit/hand-derivation):
- 2 massless legs, MASSIVE internal (gg>h): scalar C0(0,s,0;1,1,1); rank-1 1/2 B0(0)-1/2 B0(2/5);
  rank-2 1/20 B0(0)-1/20 B0(2/5) = -3.4748e-3. VALIDATED (1e-10). [baseline]
- 3 massless legs, MASSIVE internal C0(0,0,0;m,m,m): scalar C0(0,0,0;1,1,1) (finite collinear-point
  master); rank-1 -> 0; rank-2 -> 0; k^2 -> C0(0,0,0;1,1,1)+B0(0;1,1) (exact k^2=D0+m^2 identity).
  Rank>=1 external-momentum numerators -> 0 because in the fully-collinear limit all external momenta
  become proportional & shrink (q_i.q_j -> 0), so <dot(k,q_i)> projections vanish. PATH-INDEPENDENCE
  CONFIRMED: rank-2 at small OFF-shell invariants (direct path) -> 0 as delta^2, SAME for symmetric
  (0.1->-8.6e-4, 0.01->-8.4e-6) and asymmetric (0.01->-1.2e-4, 0.001->-1.2e-6). Not a delta-artifact.
- 2 massless legs, MASSLESS internal (e+e->ddx, off-shell leg s=1): scalar C0(0,s,0;0,0,0)
  [IR-divergent master, OneLOop handles]; rank-1 dot(k,q1) 1/2 B0(0;0,0)-1/2 B0(s;0,0) [hand-verified,
  B0(0;0,0)=0]; **rank-1 dot(k,q2) -> 0** (PREVIOUSLY PANICKED! both bubble pieces scaleless, =0
  correct); rank-2 1/8 B0(0;0,0)-1/8 B0(s;0,0). => the fix ALSO cures the massless-INTERNAL e+e->ddx
  case (dot(k,q2)/tensor now return correct values instead of panicking).
- ALL massless (external + internal): scalar -> C0(0,0,0;0,0,0) which is SCALELESS=0 (OneLOop evaluates
  it to 0; the reducer keeps the master symbol rather than auto-simplifying to 0); rank-2 -> 0.
  So the answer is trivially 0 in dim reg. NO special reducer needed here.

## WHO ELSE USES THIS (literature)
The fix = "regularize the degeneracy with a small parameter, reduce, take the limit" — a well-known
idea, rigorously developed:
- **Denner-Dittmaier, "Reduction schemes for one-loop tensor integrals", Nucl.Phys.B734(2006)62,
  hep-ph/0509141**: EXACTLY our case — expand tensor coefficients about limits of vanishing Gram
  determinants (the inverse-Gram-determinant / exceptional-kinematics problem) to reduce to scalars.
  They do systematic Taylor expansions in the small Gram to avoid NUMERICAL instability; we get the
  limit for free via Symbolica EXACT rational arithmetic (no instability). Same math, cleaner route.
- **Auxiliary Mass Flow / AMFlow (Liu-Ma-Wang, PRD 99 (2019) 071501; arXiv:1711.09572 + package)**:
  Eli's "small mass -> 0" intuition, rigorous & multi-loop. Introduce an auxiliary mass, compute in
  the simple (large-mass/vacuum) limit, FLOW the mass back to physical via differential equations.
  Our external-leg delta is the 1-loop symbolic analog (regularize external p^2 instead of masses).
- **Mass regularization of IR divergences** (classic: give massless particles a small mass, take
  m->0) — the textbook precursor of both.

## VERDICT on "all legs massless" + do we need vakint?
- Fully massless (external + internal) 1-loop = SCALELESS = 0 in dim reg. No computation needed;
  no special reducer required. Our reducer returns 0 (rank>=1) or a scaleless master (scalar) that
  evaluates to 0. (A cheap improvement: auto-simplify scaleless masters C0(0,0,0;0,0,0) etc. to 0.)
- vakint (and AMFlow) target MASSIVE VACUUM / multi-loop integrals (mass = the only scale, no
  external momenta). That is a DIFFERENT regime from ours. We do NOT need vakint for on-shell
  massless legs — our delta-regularization covers 2, 3, and the collinear/scaleless limits at 1 loop.
- CAVEAT (honest, math): the single-shared-delta diagonal limit gives the correct answer WHEN the
  integral is continuous at the on-shell point (holds for massive-internal finite integrals AND the
  IR-regulated massless-internal, where the eps-pole COEFFICIENTS are continuous in the invariants).
  It could fail for genuinely singular/threshold configs where different invariants must vanish at
  different rates — those need the Denner-Dittmaier systematic expansion. For all 13 tested configs
  it works. The rigorous general solution is DD/AMFlow-style.
