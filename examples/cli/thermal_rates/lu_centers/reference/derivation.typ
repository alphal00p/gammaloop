= Independent thermal spectral reference

#link("thermal_reference.py")[`thermal_reference.py`] (driver),
#link("matsubara.py")[`matsubara.py`] (exact Matsubara sums) and
#link("engine.py")[`engine.py`] (finite-$epsilon$ integrand, radial quadrature) compute the
retarded spectral kernel of the massive scalar double triangle (DT) and
triangle-box-triangle (TBT). Nothing from GammaLoop is used: no CFF or 3D expression, no LU
root or residue, no exported evaluator. The only repository inputs are the conventions of
#link("../../massive_double_triangle/benchmark.typ")[`benchmark.typ`] and
#link("../../massive_double_triangle/reference_derivation.typ")[`reference_derivation.typ`],
and the topologies of #link("../../massive_double_triangle/double_triangle.dot")[`double_triangle.dot`]
and #link("../triangle_box_triangle.dot")[`triangle_box_triangle.dot`]. Every output JSON
records `code_sha256`, the SHA-256 of the three scripts. Results are in
#link("results.json")[`results.json`], and #link("../reference.json")[`reference.json`] holds
the campaign references as `[value, total uncertainty]`.

== Model, observable, normalization

The interaction is $-lambda\/2 (X H^2 + L H^2 + X L^2)$ with $lambda=1$, $m_H=0.6$,
$m_L=0.25$, zero chemical potential and Bose statistics. Each cubic vertex gives $-lambda$
(the $1\/2$ cancels against the two identical legs), so a self-energy graph with an even
number $V$ of vertices carries $(-lambda)^V=1$. The positive Euclidean master sum-integral of
a graph with $L$ loops and lines $k$ is

$ I(K) = T^L sum_({omega_l}) integral product_l (dif^3 bold(l)_l)/(2 pi)^3
  product_k 1/(omega_k^2 + E_k^2), quad E_k = sqrt(abs(bold(k)_k)^2 + m_k^2). $

With $Pi_R (P) = -kappa I_R (P)$ and $kappa$ the `overall_factor` of the dot graph ($1\/2$ for
both graphs, from the $b <-> c$ reflection of the DT and the $v_8 <-> v_12$,
$v_9 <-> v_11$ reflection of the TBT), the observable is

$ Gamma_rho (P) = -(op("Im") Pi_R (P))/E = kappa (op("Im") I_R (E + i 0, bold(p)))/E,
  quad P = (E, bold(p)). $

Both graphs share the Feynman rules, so the DT validation also fixes the TBT normalization.
The vacuum gates check the sign convention ($i omega_K -> E + i 0$, $op("Im") I_R > 0$ for
the vacuum decay).

=== Graphs and routing

DT loops $(P, Q)$ with $K$ entering at $a$ and leaving at $d$; the momenta of `p3` and `p4`
run against the dot edge direction, which is immaterial for scalar propagators:

#table(
  columns: 4,
  table.header([Line], [Dot edge], [Momentum], [Mass]),
  [`p1`], [`0`: $a -> b$], [$P$], [$H$],
  [`p2`], [`3`: $b -> d$], [$Q$], [$L$],
  [`p3`], [`1`: $a -> c$], [$P-K$], [$H$],
  [`p4`], [`4`: $c -> d$], [$Q-K$], [$L$],
  [`p5`], [`2`: $b -> c$], [$P-Q$], [$H$],
)

TBT loops $P = e_0$, $Q = e_3$, $R = e_6$ with $K$ entering at $v_7$ and leaving at $v_10$.
The routing, masses and vertex types match `triangle_box_triangle.dot` edge by edge, and its
topology is that of
#link("../../../../../tests/resources/graphs/triangle_box_triangle_phys.dot")[`triangle_box_triangle_phys.dot`]:

#table(
  columns: 4,
  table.header([Edge], [Vertices], [Momentum], [Mass]),
  [$e_0$], [$v_7 -> v_8$], [$P$], [$H$],
  [$e_1$], [$v_7 -> v_12$], [$K-P$], [$H$],
  [$e_2$], [$v_8 -> v_12$], [$P-Q$], [$L$],
  [$e_3$], [$v_8 -> v_9$], [$Q$], [$H$],
  [$e_4$], [$v_12 -> v_11$], [$K-Q$], [$H$],
  [$e_5$], [$v_11 -> v_9$], [$R-Q$], [$H$],
  [$e_6$], [$v_9 -> v_10$], [$R$], [$L$],
  [$e_7$], [$v_11 -> v_10$], [$K-R$], [$L$],
)

== Exact Matsubara sums

At fixed spatial loop momenta write $z_l = i omega_l$, $z_K = i omega_K$ and
$1\/(omega_k^2 + E_k^2) = 1\/(E_k^2 - Omega_k^2)$ with
$Omega_k = sum_l c_(k l) z_l + c_(k K) z_K$. Each loop sum uses

$ T sum_n f(i omega_n) = - sum_(z^* in "poles"(f)) op("Res")_(z=z^*) [f(z) n_B (z)],
  quad n_B (x) = 1/(e^(x\/T) - 1), $

valid because every summand decays at least like $abs(z)^(-2)$. At the pole
$Omega_k = sigma E_k$ ($sigma = plus.minus 1$) of a propagator containing $z_l$ with
coefficient $c = plus.minus 1$,

$ -op("Res") = (sigma c)/(2 E_k) n_B (z_l^*) times (text("other factors at") z_l = z_l^*),
  quad z_l^* = c (sigma E_k - [Omega_k - c z_l]). $

$z_l^*$ is a signed sum of on-shell energies plus a combination of the other loop variables
and $z_K$, all bosonic Matsubara frequencies, so $n_B (z_l^*) = n_B (op("Re") z_l^*)$. *The
replacement $n_B (x + i omega_"ext") = n_B (x)$ is made here, before any continuation.* The
substitution is propagated into the remaining propagators, line $k$ is marked cut, and the
next loop is summed. After all $L$ sums,

$ R(z_K) = sum_t s_t product_(j in "cut") 1/(2 E_j) product_a n_B (bold(a)_t dot bold(E))
  product_(k in "real") 1/(E_k^2 - (bold(b)_(t k) dot bold(E))^2)
  product_(k in "z-dep.") 1/(E_k^2 - (c_(t k) z_K + bold(b)_(t k) dot bold(E))^2), $

with integer vectors $bold(a), bold(b)$ and $c = plus.minus 1$. This rational function of
$z_K$ equals $T^L sum product G$ at every bosonic $omega_K$. The DT has 40 terms, the TBT
368 (summation order $P, Q, R$).
#link("check_matsubara.py")[`check_matsubara.py`] compares $R(i omega_K)$ with direct
truncated sums at random momenta with $bold(p) = 0$ and $bold(p) != 0$, $T = 0.5, 0.3$ and
$omega_K = 2 pi T, 4 pi T$, removing the $n_"max"^(-5)$ truncation tail by Richardson
extrapolation. The largest relative deviations are $4.4 times 10^(-15)$ (DT, roundoff) and
$3.4 times 10^(-10)$ (TBT at $n_"max" = 48$; $1.6 times 10^(-11)$ at $n_"max" = 80$). An
independent review with its own direct sums at $omega_K = 2 pi T, 6 pi T$ found
$1.7 times 10^(-15)$ and $5.5 times 10^(-11)$; the reversed summation order agrees to
$2.1 times 10^(-14)$.

=== Retarded continuation

$R(z)$ decays for $abs(z) -> infinity$, is analytic off the real axis and coincides with the
Euclidean sum at all non-zero bosonic frequencies, so (Baym–Mermin) the retarded function is
$R(E + i 0)$. Its poles lie at $z = bold(h) dot bold(E)$, where $bold(h)$ runs over signed sums
of the line energies of the cuts separating the in- from the out-vertex (DT: 4 cuts, 24
signed functions; TBT: 9 cuts, 76 signed functions). No disconnected or non-separating
combination appears as a pole in $z$.

=== Spurious singularities

Occupations with mixed signs, such as $n_B (E_(p 1) - E_(p 5))$ in the DT, and
$z$-independent propagators that are $2 -> 2$ energy balances, such as
$E_(e 6)^2 - (E_(e 0) + E_(e 2) - E_(e 5))^2$ in the TBT, vanish on hypersurfaces of momentum
space. They cancel in the sum over terms: the Matsubara sum is finite and continuous in the
spatial momenta at every $omega_K != 0$, and a rational function of $z$ that vanishes at
infinitely many Matsubara points vanishes identically, so $R(z)$ is regular there for every
non-real $z$. Numerically the cancellation costs a relative precision of about
$10^(-16)\/delta$ at distance $delta$ from such a surface (smooth down to
$delta tilde 10^(-10)$; the review confirmed this against 60-digit arithmetic).
`engine.py` replaces any radial node whose smallest spurious denominator is below $10^(-9)$,
or whose value is non-finite, by the average at $lambda (1 plus.minus 10^(-6))$, and drops
nodes with $lambda < 10^(-6)$. There the weight $lambda^(3L-1)$ vanishes, and at
$bold(p) = 0$ all momenta vanish together so the spurious denominators are exactly
degenerate. The vacuum ($T = 0$) uses $n_B (x) = -theta(-x)$ exactly.

== Finite-$epsilon$ integrand and phase space

For $epsilon > 0$,
$J(epsilon) = kappa\/E med op("Im") integral dif^(3L) l \/ (2 pi)^(3L) R(E + i epsilon; bold(l))$
with $bold(p)$ along $z$. Loop momenta are $bold(l)_i = lambda rho_i hat(bold(n))_i$, with
$rho$ on the positive part of the unit sphere $S^(L-1)$ and unit vectors
$hat(bold(n))_i$:

- DT: $rho = (cos a, sin a)$; TBT:
  $rho = (cos a_1, sin a_1 cos a_2, sin a_1 sin a_2)$ (Jacobian $sin a_1$). The measure is
  $lambda^(3L-1) dif lambda product_i rho_i^2 dif S_(L-1)$.
- $bold(p) != 0$: the azimuth of $hat(bold(n))_1$ is fixed (factor $2 pi$), outer dimension
  $3L-2$ (DT 4, TBT 7). $bold(p) = 0$: $hat(bold(n))_1 = hat(bold(z))$ and
  $hat(bold(n))_2$ in the $x z$ plane (factors $4 pi$ and $2 pi$), outer dimension $3L-4$
  (DT 2, TBT 5).

Along a ray every line energy is $E_j (lambda) = sqrt(A_j lambda^2 + 2 B_j lambda + C_j)$.

=== Deterministic radial quadrature

For every outer point all threshold functions $h(lambda) = sum_j h_j E_j (lambda)$ (24 for
the DT, 76 for the TBT) and their derivatives are tabulated on 400 points in
$x = lambda\/(1 + lambda) <= 0.998$. Quadrature centers are

+ every root of $h(lambda) = E$ (sign change, 60 bisections), with width
  $w = 0.5 min(epsilon_min\/abs(h'), sqrt(2 epsilon_min\/abs(h'')))$;
+ every extremum of $h$ with $abs(h - E) <$ `extremum_window`, with width
  $0.5 sqrt(2 max(abs(h - E), epsilon_min)\/abs(h''))$ (near-tangent and double roots).

$[0, 1)$ in $x$ is split at the midpoints between sorted centers. On each panel
$x = c + w sinh t$, and $t$ is integrated with `nsub` equal sub-intervals of `ng`-point
Gauss–Legendre rules. The $sinh$ map resolves all scales from $w$ to the panel size, so the
same nodes serve the whole $epsilon$ ladder. Settings: DT $10 times 6$ nodes with window
0.15; TBT $8 times 4$ with window 0.15 for the thermal $(1, 0)$ run and $6 times 4$ with window
0.05 for all other TBT runs. Against $16 times 12$–$16$ nodes, a 1600–2000-point root grid and
half the widths, the DT per-ray deviations are below $10^(-5)$ ($10^(-3)$ of the mean for
single rays at $(1.25, 0.75)$, $epsilon = 0.002$) and the summed ones below
$3 times 10^(-5)$; for the TBT at $6 times 4$/0.05 they are below $3 times 10^(-4)$ of the mean
per ray and $3 times 10^(-5)$ summed. The review repeated this at four production settings:
the mean bias is below $4 times 10^(-6)$ of the statistical error. Both are far below the
statistical errors.

The outer dimensions are integrated with vegas: adaptive training iterations are discarded,
then `adapt=False` iterations are combined. The integrand returns $J(epsilon_i)$ for the whole
ladder at the same points, so vegas gives their full covariance. The first component, which
drives the adaptation, is the fitted $epsilon -> 0$ linear combination.

== Extrapolation $epsilon -> 0$

$I_R (z)$ is analytic in the upper half plane and continues analytically through the real axis
away from thresholds, so

$ J(epsilon) = J_0 + J_1 epsilon + J_2 epsilon^2 + dots, quad
  J_1 = kappa/E op("Re") I'_R (E), quad J_2 = -kappa/(2 E) op("Im") I''_R (E), $

with radius equal to the distance from $E$ to the nearest non-analytic point. Stationary
points of a signed cut energy $sum_i s_i E_i$ at fixed total momentum have all particles at a
common velocity, so (pseudo-)thresholds sit at $E_c^2 - p^2 = (sum_i s_i m_i)^2$: DT
$M_c in {0, 0.25, 0.5, 0.95, 1.2, 1.45}$, TBT additionally ${0.7, 1.7}$. The distances are
0.05 at $(1, 0)$ (from $2H - L$), 0.15 (DT) and 0.10 (TBT) at $(1.6, 0)$, 0.040 at
$(1.25, 0.75)$, 0.25 at $(0.5, 0.75)$ and 1.0 at $(0.5, 1.5)$. Timelike ladders use
$epsilon = 0.002, 0.003, 0.004, 0.006, 0.008$ ($epsilon_max\/"distance" <= 0.2$), spacelike
ones $0.004, 0.006, 0.008, 0.012, 0.016$.

The primary estimate is the unweighted quadratic least-squares fit, with its statistical error
propagated through the full covariance. The extrapolation uncertainty is the largest shift of
$J_0$ among the cubic fit, the quadratic fit without the largest $epsilon$, and the quadratic
fit without the two smallest; each shift is stored with its own correlated statistical sigma
(`shift_sigma`). Where the shift is comparable to that sigma the quoted extrapolation
uncertainty is itself statistical and therefore conservative. The correlated GLS fit is stored
but not used: neighbouring $epsilon$ are correlated at 0.98–1.000, so tiny truncation residuals
give huge $chi^2$ (e.g. $8.5 times 10^5\/2$ for DT $(1.6, 0)$) and an unstable $J_0$.

At the hardest DT point $(0.5, 0.75)$ an independent run with ten $epsilon$ from 0.0005 to
0.016 gives neighbour slopes $dif J\/dif epsilon = -1.44(33), -1.33(21), -1.24(13),
-1.19(8), -1.16(5), -1.14(3) times 10^(-3)$ from $epsilon = 0.0025$ to $0.014$: no structure
beyond a smooth quadratic (correlated quadratic fit $chi^2\/"dof" = 6.9\/7$), and
$J_0 = 3.9386(226) times 10^(-4)$ against $3.8945(160) times 10^(-4)$ from the five-point
ladder (1.6 sigma). The five-point fit shows a 2.3-sigma cubic shift there, which enters its
extrapolation uncertainty.

== Statistics and run strategy

`stat` is the propagated vegas error of the fitted $J_0$; the final-iteration $Q$ is recorded.
Independent runs of one point (different seeds and/or $epsilon$ ladders) are combined with
weights $1\/("stat"^2 + "extrap"^2)$; the propagated statistical error is multiplied by the PDG
scale factor $max(1, sqrt(chi^2\/"dof"))$, and the extrapolation uncertainties are averaged
with the same weights. The total uncertainty adds both in quadrature.

With the radial integral deterministic, the $epsilon$-dependence of the noise comes from
codimension-one sets of outer directions on which roots of two threshold functions with
opposite radial slopes coincide (principal-value structures, variance $tilde 1\/epsilon$),
and from tangencies. At $bold(p) != 0$ the thermal pseudo-thresholds have stationary
configurations at finite momenta (all cut particles co-moving with velocity $bold(p)\/E_c$),
inside the integration domain instead of at $lambda = 0$. For the TBT these collinear regions
dominate the variance: at $(1.25, 0.75)$ with uniform sampling, 64% of the variance comes from
$cos theta_p > 0.8$ and 51% from $cos theta_q > 0.8$, with a per-sample relative standard
deviation of 17–31. Consequently:

- DT: all points are at or below 0.19% except the spacelike $(0.5, 0.75)$, where three runs
  (ladders 0.004–0.016, 0.0005–0.016, 0.002–0.008; seeds 12345, 777, 4242) combine to 0.38%
  with scale factor 1.8 ($chi^2 = 6.6\/2$).
- TBT spacelike: the noise grows like $epsilon^(-0.7)$, and the standard spacelike ladder gave
  7% and 27%. Since the nearest singularity is far ($E = p$: distances 1.0 and 0.25), second
  runs use ladders 0.02–0.08 at $(0.5, 1.5)$ ($epsilon_max\/d = 0.08$) and 0.01–0.04 at
  $(0.5, 0.75)$ ($epsilon_max\/d = 0.16$), both with the quadratic fit.
- TBT $(1.25, 0.75)$: the $2H - L$ pseudo-threshold is 0.040 away, so $epsilon$ cannot grow
  and the result is statistics-limited.

== Validation

Before any new point was computed, the engine passed four gates:

+ $R(i omega_K)$ against direct Matsubara sums (above).
+ One-loop controls with the same code path (`bubble:m1:m2`, $kappa = 1$) against the
  analytic values of the #link("../../scalar_bubble/benchmark.typ")[scalar bubble]: rest
  decay, absorption, boosted decay and spacelike Landau damping.
+ The five published DT targets at $bold(p) = 0$.
+ Vacuum Lorentz invariance at $bold(p) != 0$: $Gamma_"vac" (E, p) = Gamma_"vac" (M, 0) M\/E$
  for timelike and $0$ for spacelike momenta. At the spacelike points the finite-$epsilon$
  Lorentzian tails extrapolate to $10^(-11)$–$10^(-10)$, below $10^(-6)$ of the thermal values.

#table(
  columns: 5,
  table.header([Case], [Value], [Uncertainty], [Target], [Pull]),
  [DT $(1, 0)$, $T = 0.5$], [`3.0577230e-4`], [`8.0e-9`], [`3.0576519e-4`], [+0.88],
  [DT $(1, 0)$, $T = 0.25$], [`1.4647342e-4`], [`3.3e-9`], [`1.4647685e-4`], [−1.05],
  [DT $(1, 0)$, vacuum], [`1.0980913e-4`], [`2.2e-9`], [`1.0980938e-4`], [−0.11],
  [DT $(1.6, 0)$, $T = 0.5$], [`5.2821716e-6`], [`1.0e-9`], [`5.2820877e-6`], [+0.08],
  [DT $(1.6, 0)$, vacuum], [`4.3285861e-5`], [`1.2e-9`], [`4.3284455e-5`], [+1.22],
  [Bubble $(0.5, 0.5)$, $(1, 2)$, $T = 1$], [`1.03648171e-2`], [`3.3e-8`], [`1.03647739e-2`], [+1.30],
  [Bubble $(0.4, 0.7)$, $(2.5, 1.5)$, $T = 0.6$], [`9.0940364e-3`], [`1.9e-8`], [`9.0940167e-3`], [+1.02],
  [Bubble $(0.3, 1.8)$, $(1, 0)$, $T = 0.6$], [`6.8805095e-3`], [`1.4e-9`], [`6.8805088e-3`], [+0.48],
  [Bubble $(0.4, 0.7)$, $(2, 0)$, $T = 0.6$], [`1.21045022e-2`], [`1.4e-9`], [`1.21045016e-2`], [+0.38],
  [DT vacuum $(1.25, 0.75)$], [`8.785328e-5`], [`2.8e-8`], [`8.784730e-5`], [+0.22],
  [DT vacuum $(1.76706, 0.75)$], [`3.921147e-5`], [`1.8e-8`], [`3.919357e-5`], [+0.99],
  [DT vacuum $(0.5, 0.75)$], [`-1.70e-10`], [`1.8e-10`], [0], [−0.92],
  [DT vacuum $(0.5, 1.5)$], [`-1.11e-11`], [`9.2e-12`], [0], [−1.21],
  [TBT vacuum $(1.25, 0.75)$], [`1.232308e-6`], [`3.7e-9`], [`1.228355e-6`], [+1.02],
  [TBT vacuum $(0.5, 1.5)$, review], [`-6.6e-14`], [`8.2e-14`], [0], [−0.81],
)

Bubble rows give the masses, then $(E, p)$. The independent review confirmed the routing,
masses and vertex types against both dot files, found two automorphisms fixing the external
vertices of each graph (matching $kappa = 1\/2$), and checked the large-$abs(z)$ decay of
$R$. It also evaluated the vacuum without the engine: below every threshold except $L L$,
Cutkosky factorizes $op("Im") I_"vac" = sqrt(1 - 4 m_L^2\/s)\/(16 pi) med V(s)$, with $V$ the
Euclidean vertex left of the $L L$ cut integrated over Feynman parameters. This reproduces
the published DT triangle $C_0 = 0.012747$ and gives the TBT vacuum at $(1, 0)$ as
$1.53705(2) times 10^(-6)$, against the engine's $1.53544(134) times 10^(-6)$ (0.1%,
1.2 sigma). Review reruns with other seeds and $epsilon$ ladders give
$5.28261(184) times 10^(-6)$ for DT $(1.6, 0)$ (0.29 sigma from the target),
$1.464807(34) times 10^(-4)$ for DT $(1, 0)$ at $T = 0.25$ (1.13 sigma) and
$3.941(15) times 10^(-4)$ for DT $(0.5, 0.75)$ (0.16 sigma from the result below).

== Results

$Gamma_rho$ at $T = 0.5$, with stat and extrapolation uncertainties (1 sigma, independent)
and the relative total (both in quadrature); the vacuum is given where computed. Campaign
points carry their #link("../reference.json")[`reference.json`] key.

#table(
  columns: 7,
  table.header([Graph, $(E, p)$], [Key], [$Gamma_rho$], [Stat], [Extrap.], [Rel. total], [Vacuum]),
  [DT $(1, 0)$], [`dt_rest1`], [`3.0577230e-4`], [`8.0e-9`], [`4.6e-10`], [`2.6e-5`], [`1.098091(22)e-4`],
  [DT $(1.6, 0)$], [`dt_rest16`], [`5.2821716e-6`], [`1.0e-9`], [`8.0e-12`], [`1.9e-4`], [`4.32859(12)e-5`],
  [DT $(0.5, 1.5)$], [`dt_space`], [`6.6752087e-5`], [`1.25e-7`], [`9.0e-9`], [`1.9e-3`], [0],
  [DT $(0.5, 0.75)$], [], [`3.9372175e-4`], [`1.45e-6`], [`3.2e-7`], [`3.8e-3`], [0],
  [DT $(1.25, 0.75)$], [], [`2.0619648e-4`], [`2.37e-7`], [`5.2e-9`], [`1.1e-3`], [`8.7853(28)e-5`],
  [TBT $(1, 0)$], [`tbt_rest1`], [`5.220132e-6`], [`1.35e-7`], [`2.6e-8`], [`2.6e-2`], [`1.5354(13)e-6`],
  [TBT $(1.6, 0)$], [`tbt_rest16`], [`-1.915432e-6`], [`5.7e-9`], [`4.0e-11`], [`3.0e-3`], [`-6.8799(92)e-7`],
  [TBT $(0.5, 1.5)$], [`tbt_space`], [`1.134794e-6`], [`2.66e-8`], [`4.2e-9`], [`2.4e-2`], [0],
  [TBT $(0.5, 0.75)$], [], [`6.21474e-6`], [`1.38e-6`], [`5.9e-7`], [`2.4e-1`], [],
  [TBT $(1.25, 0.75)$], [], [`3.381027e-6`], [`1.16e-7`], [`4.5e-8`], [`3.7e-2`], [`1.2323(37)e-6`],
)

The DT rest values reproduce the published targets, which
#link("../reference.json")[`reference.json`] uses with zero uncertainty. Three TBT entries
combine or correct the production runs:

- $(1.6, 0)$: the production run $(-1.92017 plus.minus 0.00360) times 10^(-6)$ and the
  review rerun $(-1.90850 plus.minus 0.00435) times 10^(-6)$ differ by 2.07 sigma. Their weighted mean
  $-1.915432 times 10^(-6)$ has the propagated error $2.77 times 10^(-9)$, which the PDG scale
  factor 2.07 raises to $5.7 times 10^(-9)$ (0.30%).
- $(1.25, 0.75)$: the production run $3.395(163)(82) times 10^(-6)$ and the review rerun
  $3.370(162)(16) times 10^(-6)$ agree (0.10 sigma) and combine to
  $3.381(116)(45) times 10^(-6)$.
- $(1, 0)$: the final vegas iterations have $Q = 0.04$ (iterations from 4.85(21) to
  5.50(20) $times 10^(-6)$), a heavy-tailed outer integrand, so the vegas error is inflated
  1.5 times to $1.35 times 10^(-7)$.

The production runs predate the final code and carry no `code_sha256`. The review reruns
used the final code, differing from the delivered files only by lint fixes and the hash
output, and reproduce them within 2.1 sigma. #link("results.json")[`results.json`] lists every
run with its seed, $epsilon$ ladder, sample counts, quadrature settings and final $Q$.

== Limitations

- DT $(0.5, 0.75)$ reaches 0.38%, not 0.3%.
- TBT precision: 0.09–0.3% for the vacuum points and the thermal $(1.6, 0)$, 2.6% at
  $(1, 0)$, 3.7% at $(1.25, 0.75)$, 2.4% at $(0.5, 1.5)$ and 24% at $(0.5, 0.75)$. The
  limit is the outer-angle variance described above, not the Matsubara sums or the radial
  quadrature. Dedicated channels (ray origins at the co-moving stationary configurations, or a
  second deterministic variable) would be the next step.
- The extrapolation assumes analyticity of $I_R (E + i epsilon)$ within the stated radius.
  The (pseudo-)threshold positions cover stationary points of single cuts; thermal anomalous
  (multi-cut) Landau singularities were not enumerated. The fine $epsilon$ scan of DT
  $(0.5, 0.75)$ and the cubic and subset fits show no sign of nearer singularities at the
  quoted precision.

== Reproduction

The scripts need Python 3 with `numpy`, `vegas` and `gvar` (the runs used Python 3.13.12,
numpy 2.4.4, vegas 6.2.1 and gvar 13.1.7). `python3 check_matsubara.py` checks the Matsubara
sums in under a minute. The driver pins BLAS threads to one; on a loaded machine the run
scripts also exported `MALLOC_MMAP_THRESHOLD_=33554432 MALLOC_TRIM_THRESHOLD_=4294967296
MALLOC_TOP_PAD_=268435456`, so that freed memory is reused instead of being faulted in again. Each run
writes its JSON and appends its log next to `--output`; `--temperature 0` gives the vacuum.
The campaign points used `--nproc 48` and the default `--degree 2 --width-factor 0.5`:

```sh
R="python3 thermal_reference.py --temperature 0.5 --nproc 48"
TL="--eps 0.002 0.003 0.004 0.006 0.008"
SL="--eps 0.004 0.006 0.008 0.012 0.016"
DT="--graph dt --ng 10 --nsub 6 --extremum-window 0.15 --nitn 6"
DT="$DT --neval-train 30000 --nitn-train 4 --batch 9600 --seed 12345"
TBT="--graph tbt --ng 6 --nsub 4 --extremum-window 0.05 --nitn 5"
TBT="$TBT --neval-train 10000 --nitn-train 6 --batch 4800"
# dt_rest1, dt_rest16, dt_space
$R $DT --energy 1 --momentum 0 $TL --neval 100000 --output dt_rest1.json
$R $DT --energy 1.6 --momentum 0 $TL --neval 200000 --output dt_rest16.json
$R $DT --energy 0.5 --momentum 1.5 $SL --neval 100000 --output dt_space.json
# tbt_rest1: 8x4 nodes, window 0.15
$R $TBT --ng 8 --extremum-window 0.15 --neval-train 12000 --seed 12345 \
  --energy 1 --momentum 0 $TL --neval 48000 --output tbt_rest1.json
# tbt_rest16: production run and review rerun (--nproc 16 in the review)
$R $TBT --energy 1.6 --momentum 0 $TL --neval 40000 --seed 12345 \
  --output tbt_rest16_a.json
$R $TBT --energy 1.6 --momentum 0 --eps 0.0025 0.0035 0.005 0.007 0.01 \
  --neval 30000 --nitn-train 5 --seed 2718 --output tbt_rest16_b.json
# tbt_space: standard spacelike ladder and large-eps ladder
$R $TBT --energy 0.5 --momentum 1.5 $SL --neval 40000 --seed 12345 \
  --output tbt_space_a.json
$R $TBT --energy 0.5 --momentum 1.5 --eps 0.02 0.03 0.04 0.06 0.08 \
  --neval 40000 --seed 426499 --output tbt_space_b.json
```

Later options override earlier ones in `$TBT`. Runs of one
point are combined as described under statistics. On the shared machine (load around 300 on
384 threads) the DT runs took 1.5–3.5 minutes and the TBT runs 13–21 minutes each.
