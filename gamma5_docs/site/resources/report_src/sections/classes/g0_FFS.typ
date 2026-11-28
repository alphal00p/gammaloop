#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Yukawa vertices $overline(Psi) Psi H$ <sec:g0_FFS>

External legs $overline(Psi)_(f_1)$, $Psi_(f_2)$ and a scalar $phi$, all momenta incoming, for the six Yukawa vertices of @eq:V-yuk:
$ (f_1, f_2, phi) = (Q, U, H^dagger), quad (Q, D, H), quad (L, E, H), \
  "conjugates:" quad (U, Q, H), quad (D, Q, H^dagger), quad (E, L, H^dagger) . $ <eq:g0_FFS-legs>
The sum of the diagrams below is $i thin Gamma^((1))_(overline(Psi)_(f_1) Psi_(f_2) phi)$, a matrix in spinor space. Every diagram has three propagators and a numerator of degree two in the loop momentum, so the superficial degree is $dg(G) = 4 + 2 - 2 dot 3 = 0$: the pole is independent of the external momenta, which are set to zero ($Top(0)$ sets them to zero anyway). The three couplings $y_u$, $y_d$, $y_e$ and all three gauge factors are covered.

As in @sec:g0_FF, the Feynman rule of the Yukawa vertex $overline(Psi)_(f_1) Psi_(f_2) phi_j$ of @eq:V-yuk is written $-i thin bold(Y)^j_(f_1 f_2) P_y$; in this subsection $bold(Y)^j equiv bold(Y)^j_(f_1 f_2)$, the coupling matrix from the colour and $"SU"(2)_L$ indices of $f_2$ to those of $f_1$, $j$ the $"SU"(2)_L$ index of the scalar field operator $phi_j in {H_j, H^dagger_j}$, and $P_y = P_(overline(c)_(f_1)) = P_(c_(f_2))$ ($PR$ for the first three vertices, $PL$ for the conjugates). For example $(bold(Y)^j_(Q U))^(alpha i)_(space beta) = y_u epsilon_(i j) delta^alpha_beta$ for $overline(Q) U H^dagger$ and $(bold(Y)^j_(Q D))^(alpha i)_(space beta) = y_d delta_(i j) delta^alpha_beta$ for $overline(Q) D H$. The tree vertex function is $Gamma^((0))_(overline(Psi)_(f_1) Psi_(f_2) phi_i) = -bold(Y)^i P_y$. The generator on the scalar leg is $Tg^a_phi = Tg^a_H$ for $phi = H$ and $Tg^a_(H^dagger) = -(Tg^a_H)^*$ for $phi = H^dagger$, so that rigid invariance of @eq:LYuk reads
$ Tg^a_(f_1) bold(Y)^i - bold(Y)^i Tg^a_(f_2) = bold(Y)^j (Tg^a_phi)_(j i) . $ <eq:g0_FFS-rigid>
Contracting @eq:g0_FFS-rigid once more with the generators gives $sum_a [Tg^a_(f_1) bold(Y)^j - bold(Y)^j Tg^a_(f_2)] (Tg^a_phi)_(j i) = C_H bold(Y)^i$, i.e.
$ 2 sum_a Tg^a_(f_1) bold(Y)^i Tg^a_(f_2) = (C_(f_1) + C_(f_2) - C_H) bold(Y)^i equiv 2 t_k thin bold(Y)^i , $ <eq:g0_FFS-casimir>
with $C_f$ of @eq:invariants ($C_(H^dagger) = C_H$); $t_k$ is defined by this equation for the vertex $k$.

*Diagrams.* At one loop the vertex carrying the external scalar has exactly two internal lines. This allows only the Yukawa vertices and the cubic $H^dagger H A$ vertex: the seagull $H^dagger H A A$ and the quartic vertex would need three internal lines, and there is no $overline(Psi) Psi H H$ or $overline(Psi) Psi A H$ vertex. The one-loop 1PI graphs are therefore the triangles of @fig:g0_FFS-diagrams: gauge exchange between the fermion legs (a), gauge exchange between a fermion leg and the scalar leg (b1, b2), and scalar exchange between the fermion legs (c).

#figure(
  diagrow(
    oneloop(
      ((legs: (("S", $phi$, 0),)), (legs: (("F", $overline(Psi)_(f_1)$, -1),)), (legs: (("F", $Psi_(f_2)$, 1),))),
      (("F", $Psi_(f_1)$, 1), ("V", $A^a$, 0), ("F", $Psi_(f_2)$, 1)), caption: [(a)], scale: 0.85),
    oneloop(
      ((legs: (("F", $overline(Psi)_(f_1)$, -1),)), (legs: (("S", $phi$, 0),)), (legs: (("F", $Psi_(f_2)$, 1),))),
      (("S", $H$, 0), ("V", $A^a$, 0), ("F", $Psi_(f_2)$, 1)), caption: [(b1)], scale: 0.85),
    oneloop(
      ((legs: (("F", $Psi_(f_2)$, 1),)), (legs: (("F", $overline(Psi)_(f_1)$, -1),)), (legs: (("S", $phi$, 0),))),
      (("F", $Psi_(f_1)$, 1), ("V", $A^a$, 0), ("S", $H$, 0)), caption: [(b2)], scale: 0.85),
    oneloop(
      ((legs: (("S", $H^dagger$, 0),)), (legs: (("F", $overline(Q)$, -1),)), (legs: (("F", $U$, 1),))),
      (("F", $D$, 1), ("S", $H$, -1), ("F", $Q$, 1)), caption: [(c), $overline(Q) U H^dagger$], scale: 0.85),
  ),
  caption: [Topologies of $Gamma^((1))_(overline(Psi)_(f_1) Psi_(f_2) phi)$. In (b1) the gauge boson joins the $Psi_(f_2)$ leg to the scalar line, in (b2) the $overline(Psi)_(f_1)$ leg. (c) is drawn for the up-type vertex; the internal scalar is $H$ or $H^dagger$ according to the vertices.],
) <fig:g0_FFS-diagrams>

The internal fermions of (c) follow from @eq:V-yuk. For $overline(Q) U H^dagger$: the vertex that absorbs $U$ can only be $overline(Q) U H^dagger$ (internal $Q$, internal scalar from $H^dagger$); the vertex that emits $overline(Q)$ must then contain $H$, i.e. it is $overline(Q) D H$ (internal $D$); the middle vertex $overline(D) Q H^dagger$ carries the external $H^dagger$. A $y_u^3$ graph does not exist, because its two outer vertices would both contain $H^dagger$. For $overline(Q) D H$ the same argument gives the outer vertices $overline(Q) D H$, $overline(Q) U H^dagger$ and the middle vertex $overline(U) Q H$. For $overline(L) E H$ both outer vertices would be $overline(L) E H$ with two fields $H$: there is no graph (c).

#figure(
  table(
    columns: (auto, auto, auto, auto, 1fr),
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [sym. factor], [sign], [group factor]),
    [(a)], [$Psi_(f_1)$, $Psi_(f_2)$, $A^a$ ($a$ over all 12)], [$1$], [$+1$], [$sum_a Tg^a_(f_1) bold(Y)^i Tg^a_(f_2) = t_k bold(Y)^i$],
    [(b1)], [$Psi_(f_2)$, $A^a$, scalar], [$1$], [$+1$], [$sum_a bold(Y)^j (Tg^a_phi)_(j i) Tg^a_(f_2) = (t_k - C_(f_2)) bold(Y)^i$],
    [(b2)], [$Psi_(f_1)$, $A^a$, scalar], [$1$], [$+1$], [$sum_a Tg^a_(f_1) bold(Y)^j (Tg^a_phi)_(j i) = (C_(f_1) - t_k) bold(Y)^i$],
    [(c), $overline(Q) U H^dagger$], [$Psi_Q$, $Psi_D$, $H$], [$1$], [$+1$], [$y_u y_d^2 epsilon_(j i') = - y_d^2 (bold(Y)^j_(Q U))^(i')$],
    [(c), $overline(Q) D H$], [$Psi_Q$, $Psi_U$, $H$], [$1$], [$+1$], [$y_d y_u^2 (epsilon epsilon)_(i' j) = - y_u^2 (bold(Y)^j_(Q D))^(i')$],
    [(c), $overline(L) E H$], [none], [], [], [],
  ),
  caption: [Field assignments of @fig:g0_FFS-diagrams. No closed loop, no identical lines. The group factors of (b1) and (b2) follow from @eq:g0_FFS-rigid and @eq:g0_FFS-casimir. In (c), $i'$ is the $"SU"(2)_L$ index of $overline(Q)$ and $j$ that of the external scalar; colour gives $delta^alpha_beta$. The conjugate vertices have the same factors with (b1) and (b2) exchanged in the roles of the legs.],
) <tab:g0_FFS-diagrams>

*Integrands.* All external momenta vanish. In (a) the loop momentum $k$ flows along both fermion lines in the direction of the arrow, the gauge boson carries $-k$ from the $Psi_(f_2)$ vertex to the $overline(Psi)_(f_1)$ vertex. In (b1) $k$ flows along $Psi_(f_2)$ into the Yukawa vertex, in (b2) along $Psi_(f_1)$ out of it; the gauge boson and the scalar close the loop. In (c) $k$ flows along both fermion lines, the scalar carries $-k$. Writing the outgoing end of the fermion line on the left,
$ G_"(a)" = intk (i gbar^mu P_(c_(f_1)) Tg^a_(f_1)) (i sl(k))/(k^2) (-i bold(Y)^i P_y) (i sl(k))/(k^2) (i gbar^nu P_(c_(f_2)) Tg^a_(f_2)) (-i g_(mu nu))/(k^2) , $ <eq:g0_FFS-int-a>
$ G_"(b1)" = intk (-i bold(Y)^j P_y) (i sl(k))/(k^2) (i gbar^nu P_(c_(f_2)) Tg^a_(f_2)) (-i g_(nu mu))/(k^2) (-i k^mu (Tg^a_phi)_(j i)) i/(k^2 - m^2) , $ <eq:g0_FFS-int-b1>
$ G_"(b2)" = intk (i gbar^mu P_(c_(f_1)) Tg^a_(f_1)) (i sl(k))/(k^2) (-i bold(Y)^j P_y) (-i g_(mu nu))/(k^2) (i k^nu (Tg^a_phi)_(j i)) i/(k^2 - m^2) . $ <eq:g0_FFS-int-b2>
The scalar-vertex factors come from @eq:V-ssv. For $phi = H^dagger$ in (b1) the internal $H$ enters the scalar vertex with momentum $k$ and the external $H^dagger$ leaves with $0$, giving $i k^mu (Tg^a_H)_(i j) = -i k^mu (Tg^a_(H^dagger))_(j i)$; for $phi = H$ the internal line leaves the scalar vertex with $-k$, giving $-i k^mu (Tg^a_H)_(j i)$. Both are $-i k^mu (Tg^a_phi)_(j i)$, and in (b2) the momentum is reversed. For (c), up-type, with the vertices of @eq:V-yuk ($overline(Q)^(alpha i') D H_n$, $overline(D) Q_l H^dagger_j$, $overline(Q)^l U_beta H^dagger_(n')$), the internal $Q$ line with the single $"SU"(2)_L$ index $l$ at both ends, and the $H$ propagator $i delta_(n n') \/ (k^2 - m^2)$,
$ G_"(c)" = intk (-i y_d delta_(i' n) PR) (i sl(k))/(k^2) (-i y_d delta_(l j) PL) (i sl(k))/(k^2) (-i y_u epsilon_(l n') PR) i/(k^2 - m^2) thin delta_(n n') delta^alpha_beta , $ <eq:g0_FFS-int-c>

*$Rop$-operation.* At one loop $Rbar G = G$. _Prefactors_: $i dot i dot (-i) dot i dot i dot (-i) = -1$ in (a), $(-i)(i)(i)(-i)(-i)(i) = +1$ in (b1), $(i)(i)(-i)(-i)(i)(i) = -1$ in (b2) and $(-i)^3 i^3 = +1$ in (c). _Dirac algebra_: every fermion propagator stands between two physical vertices, so by @eq:phys-prop and @eq:sandwich
$ "(a)": quad gbar^mu P_(c_(f_1)) sl(k) P_y sl(k) gbar^nu P_(c_(f_2)) g_(mu nu) = gbar^mu slb(k) slb(k) gbar_mu P_y = 4 bar(k)^2 P_y, quad
  "(b1)": quad P_y sl(k) gbar^nu P_(c_(f_2)) k_nu = bar(k)^2 P_y, \
  "(b2)": quad gbar^mu P_(c_(f_1)) sl(k) P_y thin k_mu = bar(k)^2 P_y, quad
  "(c)": quad PR sl(k) PL sl(k) PR = bar(k)^2 PR , $ <eq:g0_FFS-dirac>
using $gbar^mu gbar_mu = 4$ and $gbar^nu k_nu = slb(k)$. The hatted parts of $sl(k)$ drop out between the projectors; only $gbar$ and $bar(k)$ remain. _Expansion_: $Top(0)$ keeps the degree-zero term, i.e. every denominator becomes $k^2 - MUV^2$ and $m^2$ is dropped. With @eq:tensor-red-split ($bar(k)^2 -> 4 k^2 \/ d$, an unsplit pairing without $eta$) and @eq:logpole,
$ Top(0) intk (bar(k)^2)/(k^2 thin k^2 (k^2 - m^2)) = 4/d intk (k^2)/((k^2 - MUV^2)^3) = ilf 1/ep + O(ep^0) . $ <eq:g0_FFS-T0>
No trace $gmh^mu_(space mu)$ occurs, so $N_1 = 0$ for all four topologies. With the group factors of @tab:g0_FFS-diagrams, and $N_0$ written as a multiple of the tree vertex function $Gamma^((0)) = -bold(Y)^i P_y$,
$ N_0^"(a)" = -4 sum_a Tg^a_(f_1) bold(Y)^i Tg^a_(f_2) P_y = 4 t_k thin Gamma^((0)), quad
  N_0^"(b1)" = sum_a bold(Y)^j (Tg^a_phi)_(j i) Tg^a_(f_2) P_y = (C_(f_2) - t_k) Gamma^((0)), \
  N_0^"(b2)" = - sum_a Tg^a_(f_1) bold(Y)^j (Tg^a_phi)_(j i) P_y = (C_(f_1) - t_k) Gamma^((0)), quad
  N_0^"(c)" = y'_k^2 thin Gamma^((0)) , $ <eq:g0_FFS-N>
where $y'_k$ (a symbol of this subsection) is the coupling of the partner vertex: $y'_u = y_d$, $y'_d = y_u$, $y'_e = 0$. The sign of $N_0^"(c)"$ comes from $epsilon_(j i') = -epsilon_(i' j)$ (up) and $(epsilon epsilon)_(i' j) = -delta_(i' j)$ (down): the internal scalar line links an $epsilon$ vertex to a $delta$ vertex. For a gluon $Tg^a_phi = 0$, so (b1) and (b2) vanish, consistent with $t_k = C_(f_1) = C_(f_2)$ on the colour part. The conjugate vertices give the same $N_0$ per gauge factor, with the values of (b1) and (b2) interchanged.

_No evanescent pole._ The pole has dimension zero, has no external momentum and no external index, and sits between $P_(overline(c)_(f_1))$ and $P_(c_(f_2))$ with only barred matrices inside: it is proportional to $P_y$. Components with an external sterile partner vanish for the reason given in @sec:g0_FF.

*Result.*
#result[
  $ Gamma^((1))_(overline(Psi)_(f_1) Psi_(f_2) phi_i) |_"pole" = lf 1/ep N_0, \
    N_0 = [-2 sum_a Tg^a_(f_1) bold(Y)^i Tg^a_(f_2) - (C_(f_1) + C_(f_2)) bold(Y)^i - y'_k^2 bold(Y)^i] P_y = z_k thin Gamma^((0))_(overline(Psi)_(f_1) Psi_(f_2) phi_i) , $ <eq:g0_FFS-pole>
  with
  $ z_k = 2 t_k + C_(f_1) + C_(f_2) + y'_k^2 = 2 (C_(f_1) + C_(f_2)) - C_H + y'_k^2 = C_H + 4 t_k + y'_k^2 , $ <eq:g0_FFS-nu>
  identical for a vertex and its conjugate; $z_k$ is a symbol of this subsection. In the last form $4 t_k$ is (a) and $C_H$ is the sum of (b1) and (b2). The counterterm vertex function is $-lf z_k Gamma^((0)) \/ ep$. Since $Gamma^((0))$ is the vertex function of $y_k Q^+_(y_k) = -y_k integral (O_k + O_k^dagger)$ (@eq:Ok, @tab:Qbasis),
  $ delta S^((1)) supset - lf 1/ep integral dif^d x thin [ - z_k y_k (O_k + O_k^dagger) ] = - lf 1/ep thin z_k y_k thin Q^+_(y_k), quad k = u, d, e . $ <eq:g0_FFS-ct>
  $O_k$ is continued to $d$ dimensions unchanged (it contains only $P_(c_f) Psi_f$, $overline(Psi)_f P_(overline(c)_f)$). There is no $Q^-_(y_k)$ and no evanescent operator. With $y_k -> y_k + delta y_k$, $psi_f -> (1 + delta Z_f \/ 2) psi_f$ and $H -> (1 + delta Z_H \/ 2) H$ in $S_0$, the coefficient of $Q^+_(y_k)$ is
  $ delta y_k + y_k / 2 (delta Z_(f_1) + delta Z_(f_2) + delta Z_H) = - lf 1/ep thin z_k y_k , $ <eq:g0_FFS-rel>
  as in the assembled counterterm of @sec:st-solution. Here $delta Z_f = -lf (C_f + cal(Y)_f) \/ ep$ is the coefficient of $Q^f_K$ (@sec:g0_FF; $cal(Y)_f$ and $tilde(Y) = 3 y_u^2 + 3 y_d^2 + y_e^2$ are the Yukawa sums of @sec:st-solution), and $delta Z_H = lf (2 C_H - tilde(Y)) \/ ep$ is the coefficient of the #box[$d$-dimensional] $Q_(S S)$ (@sec:g0_SS; the evanescent Higgs poles $E_H$, $tilde(E)_H$ found there do not enter @eq:g0_FFS-rel). Hence, for $delta y_k$ of @sec:st-solution,
  $ delta y_k = lf 1/ep thin y_k [ - 3/2 (C_(f_1) + C_(f_2)) + 1/2 (cal(Y)_(f_1) + cal(Y)_(f_2) + tilde(Y)) - y'_k^2 ] . $ <eq:g0_FFS-dy>
]

#figure(
  table(
    columns: 6,
    align: (center, left, left, left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$k$], [(a)], [(b1)], [(b2)], [(c)], [$z_k$]),
    [$u$], [$16/3 gs^2 + 4/9 gy^2$], [$1/3 gy^2$], [$3/4 gw^2 - 1/12 gy^2$], [$y_d^2$], [$16/3 gs^2 + 3/4 gw^2 + 25/36 gy^2 + y_d^2$],
    [$d$], [$16/3 gs^2 - 2/9 gy^2$], [$1/6 gy^2$], [$3/4 gw^2 + 1/12 gy^2$], [$y_u^2$], [$16/3 gs^2 + 3/4 gw^2 + 1/36 gy^2 + y_u^2$],
    [$e$], [$2 gy^2$], [$1/2 gy^2$], [$3/4 gw^2 - 1/4 gy^2$], [$0$], [$3/4 gw^2 + 9/4 gy^2$],
  ),
  caption: [Standard-Model values of $N_0 \/ Gamma^((0))$ per topology for the vertices $overline(Q) U H^dagger$, $overline(Q) D H$, $overline(L) E H$ (@eq:g0_FFS-N); for the conjugates (b1) and (b2) are interchanged. Per gauge factor $t_k$ is $4/3 gs^2$ (quarks) plus $gy^2 Y_(f_1) Y_(f_2)$: $t_u = 4/3 gs^2 + 1/9 gy^2$, $t_d = 4/3 gs^2 - 1/18 gy^2$, $t_e = 1/2 gy^2$; the $W$ does not couple to $f_2$.],
) <tab:g0_FFS-values>

#figure(
  table(
    columns: 2,
    align: (center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$k$], [$16 pi^2 ep thin delta y_k \/ y_k$]),
    [$u$], [$-4 gs^2 - 9/8 gw^2 - 17/24 gy^2 + 9/4 y_u^2 + 3/4 y_d^2 + 1/2 y_e^2$],
    [$d$], [$-4 gs^2 - 9/8 gw^2 - 5/24 gy^2 + 3/4 y_u^2 + 9/4 y_d^2 + 1/2 y_e^2$],
    [$e$], [$-9/8 gw^2 - 15/8 gy^2 + 3/2 y_u^2 + 3/2 y_d^2 + 5/4 y_e^2$],
  ),
  caption: [Yukawa counterterms @eq:g0_FFS-dy in Feynman gauge, $overline("MS")$.],
) <tab:g0_FFS-dy>

#check[*Auxiliary mass, Higgs mass and $eta$.* $MUV$ and $m^2$ are of degree two and do not survive $Top(0)$; $N_0$ is a pure number times $Gamma^((0))$. No hatted contraction occurs, so $Kop = Kop_eta$ and $N_1 = 0$ (@sec:m-K). The coefficients are real and equal for $O_k$ and $O_k^dagger$, as required by hermiticity; no $Q^-_(y_k)$ is generated.]

#check[*Gluon part and anticommuting limit.* The numerators @eq:g0_FFS-dirac contain no hatted object, so @eq:g0_FFS-pole coincides with the result for an anticommuting $g5$. For $gs$ alone only (a) contributes, $z_k = 4 t_k = 4 C_F gs^2$ with $C_F = 4\/3$, and @eq:g0_FFS-dy gives $(delta y_k) \/ y_k = - lf 3 C_F gs^2 \/ ep$, the QCD mass renormalisation $delta Z_m$. With $beta_(y_k) = 2 y_k dot (16 pi^2 ep thin delta y_k \/ y_k) \/ (16 pi^2)$ at one loop, @tab:g0_FFS-dy reproduces the standard one-loop Standard-Model values, e.g. $16 pi^2 beta_(y_u) = y_u (9/2 y_u^2 + 3/2 y_d^2 + y_e^2 - 8 gs^2 - 9/4 gw^2 - 17/12 gy^2)$ and the gauge parts $-5/12 gy^2$ ($d$), $-9/4 gw^2 - 15/4 gy^2$ ($e$). The gauge part of @eq:g0_FFS-dy is $-3/2 (C_(f_1) + C_(f_2))$: the $C_H$ of the vertex (@eq:g0_FFS-nu) cancels against that of $delta Z_H$.]

#check[*Arithmetic.* @eq:g0_FFS-N and @tab:g0_FFS-values were reproduced for all six vertices and every $"SU"(2)_L$ component of the external scalar (the 32 multiplet-level diagrams counted in @tab:classes) in two independent ways: by contractions of the Feynman rules with explicit Dirac matrices and $"SU"(3) times "SU"(2) times "U"(1)$ generators, without using @eq:g0_FFS-rigid, and by the symbolic BMHV algebra with #box[$d$-dimensional] external momenta, which confirms that $N_0$ has no momentum-dependent or hatted part (@sec:supporting).]

#check[*Literature.* For a single gauge group the gauge part $C_H + 4 t_k$ of @eq:g0_FFS-nu is eq. (5.14) of @Belusca-Maito:2020ala at $xi = 1$, whose gauge terms $xi thin C_2(S) + (3 + xi) T_(overline(R)) Y T_R \/ Y$ (relative to the tree vertex) become $C_H + 4 t_k$ with $T_(overline(R)) Y T_R -> sum_a Tg^a_(f_1) bold(Y)^i Tg^a_(f_2)$. The fermion self-energy entering @eq:g0_FFS-rel is eq. (5.10) of the same reference at $xi = 1$.]

#check[*Relation to other classes.* The pole @eq:g0_FFS-ct combines with @sec:g0_FF and @sec:g0_SS into the invariant combination $delta y_k Q^+_(y_k) + sum_f (delta Z_f \/ 2) cal(O)^f_psi + (delta Z_H \/ 2) cal(O)_H + dots$ of @eq:closed, which is part of the assembled counterterm of @sec:st-solution. $bop_d Q^+_(y_k) = 0$, because $O_k$ contains only physical components and is BRST invariant in $d$ dimensions; @eq:g0_FFS-ct therefore enters the pole check of @sec:g1_cFFS only through the Yukawa parts of $bop_d$ of the fermion and Higgs counterterms of @sec:g0_FF, @sec:g0_FFV and @sec:g0_SS. The finite Yukawa restoring counterterms are determined in @sec:st-solution.]

#general[At $L$ loops the Yukawa pole remains proportional to $P_y$ with barred matrices inside the projectors, so no evanescent Yukawa counterterm arises at any order in the ghost-number-zero sector: the external legs attach to physical components (no counterterm contains $chi_f$), and at superficial degree $0$ all Lorentz indices of the pole are contracted internally, so the string between the projectors reduces to a number. The forest of @eq:m-Rbar subtracts the one-loop counterterms of @sec:g0_FF, @sec:g0_FFV, @sec:g0_SS and of this class; insertions of $E^((r))_A$, $E_H$ and of $Delta S^((1))$ into internal lines produce $eta$-dependent finite terms, kept as in @sec:m-beyond. From three loops on, in bookkeeping B, the ghost-number-0 subgraphs with $eta$-dependent poles are subtracted with $delta S^((L))_eta$ and the remainder $E^((L))$ is inserted as a vertex (@eq:m-LPl).]
