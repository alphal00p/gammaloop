#import "../../macros.typ": *
#import "../../diagrams.typ": *

// local shorthands (this file only)
#let Yt = $tilde(Y)$
#let br(x, y) = $⟨#x, #y⟩$

=== Breaking of $omega H^dagger H$, $omega H^dagger H A$ and higher <sec:g1_cSS_cSSV>

External legs $gh^a (q)$, $H_j (p_1)$, $H^dagger_i (p_2)$ and, where present, $A^c_nu (p_3)$, $A^e_lambda (p_4)$, all momenta incoming and $sum p = 0$; $H_j$ carries $p_1$ and $H^dagger_i$ carries $p_2$, as in @eq:Nhat. The sum of the diagrams below is $i thin (Dbrk_1 dot Gamma)^((1))$ for these legs. Superficial degrees of divergence: $dg(G) = 2$ for $gh H^dagger H$, $1$ for $gh H^dagger H A$, $0$ for $gh H^dagger H A A$ and $gh (H^dagger H)^2$. Ghost and gauge field run over all twelve gauge indices; the three Yukawa couplings contribute.

In this subsection @eq:LYuk is written as $cal(L)_"Yuk" = - sum_k y_k [ overline(Psi)_(f') thin Phi_k P_(c_f) Psi_f + overline(Psi)_f thin Phi_k^dagger P_(c_(f')) Psi_(f') ]$ with
$ (f, f', Phi_k) = (D, Q, H), quad (E, L, H), quad (Q, U, tilde(H)^dagger) quad "for" k = d, e, u, quad Phi_k = Lambda_(k j) H_j . $ <eq:g1_cSS_cSSV-Phi>
$f$ enters the $H$ vertex and $f'$ leaves it, as in @sec:g0_SS. $Lambda_(k j)$ is a constant matrix from the colour and $"SU"(2)_L$ indices of $f$ to those of $f'$: $(Lambda_(d j))_(alpha l, beta) = delta_(l j) delta_(alpha beta)$, $(Lambda_(u j))_(alpha, beta l) = epsilon_(l j) delta_(alpha beta)$, and $Lambda_(e j) = delta_(l j)$; $tilde(H)^dagger_l = epsilon_(l j) H_j$. In the notation of @sec:g0_FFS, $y_k Lambda_(k j) = bold(Y)^j$ for the vertex $(f_1, f_2, phi) = (f', f, H)$, i.e. $(Q, D, H)$, $(L, E, H)$, $(U, Q, H)$. For matrices $X$ on $f'$ and $Z$ on $f$ define the doublet matrix
$ br(X, Z)_(i j) = Tr(Lambda_(k i)^dagger thin X thin Lambda_(k j) thin Z), quad quad br(1, 1)_(i j) = N^c_k delta_(i j), quad quad N^c_u = N^c_d = 3, quad N^c_e = 1 . $ <eq:g1_cSS_cSSV-br>
Rigid invariance of the Yukawa terms is the covariance
$ Tg^a_(f') Lambda_(k j) - Lambda_(k j) Tg^a_f = Lambda_(k l) (Tg^a_H)_(l j), quad "hence" quad br(X Tg^c, Z) = br(X, Tg^c Z) + br(X, Z) Tg^c_H, $ <eq:g1_cSS_cSSV-cov>
with the doublet matrix product on the right. $N^c_k$ is the colour multiplicity $N^c_f = N^c_(f')$ of the two multiplets of the pair $k$. The Yukawa sum $Yt = sum_k N^c_k y_k^2 = 3 y_u^2 + 3 y_d^2 + y_e^2$ is that of @sec:g0_SS and @eq:sts-yuk.

*Diagrams.* The $Dbrk_1$ vertex has two fermion legs, so every one-loop graph is a fermion loop that carries the $Dbrk_1$ vertex, the Yukawa vertices of the external Higgs legs and the gauge vertices @eq:V-ffv of the external gauge legs. There is no Higgs propagator. @fig:g1_cSS_cSSV-diagrams shows the topologies, @tab:g1_cSS_cSSV-diagrams the assignments.

#figure(
  stack(dir: ttb, spacing: 4mm,
  diagrow(
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("S", $H^dagger_i$, -1),)), (legs: (("S", $H_j$, 1),))),
      (("F", $f'$, 1), ("F", $f$, 1), ("F", $f'$, 1)), caption: [(a) triangle, $Dbrk_1$ on $f'$]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("S", $H^dagger_i$, -1),)), (legs: (("S", $H_j$, 1),)), (legs: (("V", $A^c_nu$, 0),))),
      (("F", $f'$, 1), ("F", $f$, 1), ("F", $f'$, 1), ("F", $f'$, 1)), caption: [(b) box, $A$ next to $Dbrk_1$]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("S", $H^dagger_i$, -1),)), (legs: (("V", $A^c_nu$, 0),)), (legs: (("S", $H_j$, 1),))),
      (("F", $f'$, 1), ("F", $f$, 1), ("F", $f$, 1), ("F", $f'$, 1)), caption: [(c) box, $A$ opposite]),
  ),
  diagrow(
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("S", $H^dagger$, -1),)), (legs: (("V", $A$, 0),)), (legs: (("V", $A$, 0),)), (legs: (("S", $H$, 1),))),
      (("F", none, 1), ("F", none, 1), ("F", none, 1), ("F", none, 1), ("F", none, 1)), caption: [(d) $gh H^dagger H A A$]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("S", $H^dagger$, -1),)), (legs: (("S", $H$, 1),)), (legs: (("S", $H^dagger$, -1),)), (legs: (("S", $H$, 1),))),
      (("F", none, 1), ("F", none, 1), ("F", none, 1), ("F", none, 1), ("F", none, 1)), caption: [(e) $gh (H^dagger H)^2$]),
  )),
  caption: [Topologies, one placement each. The $Dbrk_1$ vertex (crossed circle) sits on the line $f'$ or $f$ of the Yukawa pair; in (b) the gauge vertex is on the same line as $Dbrk_1$, in (c) on the other line. (d) and (e) stand for all placements of the vertices.],
) <fig:g1_cSS_cSSV-diagrams>

#figure(
  table(
    columns: 5,
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines, placement], [symmetry factor], [sign], [group factor]),
    [(a1)], [$f, f', f'$; $Dbrk_1$ on $f'$], [$1$], [$-1$ (fermion loop)], [$y_k^2 br(Tg^a, 1)$],
    [(a2)], [$f', f, f$; $Dbrk_1$ on $f$], [$1$], [$-1$], [$y_k^2 br(1, Tg^a)$],
    [(b1)], [$Dbrk_1$ and $A$ on $f'$, two orders], [$1$], [$-1$], [$y_k^2 br(Tg^a Tg^c, 1)$, $y_k^2 br(Tg^c Tg^a, 1)$],
    [(b2)], [$Dbrk_1$ and $A$ on $f$, two orders], [$1$], [$-1$], [$y_k^2 br(1, Tg^a Tg^c)$, $y_k^2 br(1, Tg^c Tg^a)$],
    [(c1)], [$Dbrk_1$ on $f'$, $A$ on $f$], [$1$], [$-1$], [$y_k^2 br(Tg^a, Tg^c)$],
    [(c2)], [$Dbrk_1$ on $f$, $A$ on $f'$], [$1$], [$-1$], [$y_k^2 br(Tg^c, Tg^a)$],
    [(d), (e)], [all placements, all $k$], [$1$], [$-1$], [vanish, see below],
  ),
  caption: [Field assignments for each Yukawa pair $k = u, d, e$ of @eq:g1_cSS_cSSV-Phi; $Tg^a = Tg^a_(f')$ in the first and $Tg^a_f$ in the second slot of @eq:g1_cSS_cSSV-br. The $Dbrk_1$ vertex carries the two terms of @eq:V-delta; both are part of each row. All loops are oriented by the Yukawa vertices, so every symmetry factor is $1$. Multiplet-level graphs with a non-zero group factor: $3$ ($gh^W H^dagger H$), $6$ ($gh^B H^dagger H$); $12$ ($gh^G H^dagger H G$), $6$ ($gh^W H^dagger H W$), $9$ each ($gh^W H^dagger H B$, $gh^B H^dagger H W$), $18$ ($gh^B H^dagger H B$).],
) <tab:g1_cSS_cSSV-diagrams>

*Integrands.* Every vertex other than $Dbrk_1$ is of the form $P_(overline(c)_("out")) Gamma P_(c_("in"))$. On the line of $Dbrk_1$ (multiplet $f$, incoming momentum $p$, outgoing $p' = p + q$) the vertex @eq:V-delta and its two propagators therefore always appear between $P_(c_f)$ and $P_(overline(c)_f)$, and with @eq:sandwich
$ P_(c_f) (i sl(p)')/(p'^2) Tg^a_f (slh(p) P_(overline(c)_f) - slh(p)' P_(c_f)) (i sl(p))/(p^2) P_(overline(c)_f)
  = - Tg^a_f (slb(q) thin hat(p)^2 - (2 hat(p) dot hat(q) + hat(q)^2) slb(p))/(p'^2 p^2) P_(overline(c)_f) . $ <eq:g1_cSS_cSSV-E>
The two terms of @eq:V-delta give $-Tg^a slb(p)' hat(p)^2$ and $+Tg^a hat(p')^2 slb(p)$; they are the two forms of @eq:m-sterile, each multiplied by the physical propagator of the other line. The right side is exactly of first order in $q$: at zero ghost momentum the $Dbrk_1$ vertex with its two propagators vanishes.

(a) For (a1) the loop momentum $k$ runs on the line $f$ from the $H^dagger$ vertex to the $H$ vertex; $p = k + p_1$ enters $Dbrk_1$ and $p' = k - p_2$ leaves it. Reading the loop against the fermion flow from the $H^dagger$ vertex, with $-i y_k P_(c_(f'))$ and $-i y_k P_(c_f)$ from @eq:V-yuk,
$ G_"(a1)" = - intk tr[ (-i y_k P_(c_(f'))) thin bold(E)_(f') (p', p) thin (-i y_k P_(c_f)) (i sl(k))/(k^2) ] = -2 i y_k^2 br(Tg^a, 1)_(i j) F(p_1, p_2) , $ <eq:g1_cSS_cSSV-int-a1>
where $bold(E)_(f')$ is the left side of @eq:g1_cSS_cSSV-E for $f'$, $tr(slb(v) slb(k) P) = 2 bar(v) dot bar(k)$ for either projector, and
$ F(v_1, v_2) = intk ( (bar(p)' dot bar(k)) hat(p)^2 - hat(p')^2 (bar(p) dot bar(k)) )/(k^2 thin p^2 thin p'^2), quad p = k + v_1, quad p' = k - v_2 , $ <eq:g1_cSS_cSSV-F>
with $v_1$, $v_2$ two external momenta.
For (a2) $k$ runs on $f'$ from the $H$ to the $H^dagger$ vertex, $p = k + p_2$, $p' = k - p_1$ on $f$, and the same steps give $G_"(a2)" = -2 i y_k^2 br(1, Tg^a)_(i j) F(p_2, p_1)$. With $k -> -k$, $F(v_2, v_1) = -F(v_1, v_2)$, and by @eq:g1_cSS_cSSV-cov, $br(Tg^a, 1) - br(1, Tg^a) = N^c_k Tg^a_H$:
$ G_"(a1)" + G_"(a2)" = -2 i y_k^2 (br(Tg^a, 1) - br(1, Tg^a))_(i j) F(p_1, p_2) = -2 i N^c_k y_k^2 (Tg^a_H)_(i j) F(p_1, p_2) . $ <eq:g1_cSS_cSSV-a>

(b), (c) The box has the four vertices $Dbrk_1$, $H$, $H^dagger$, $A^c_nu$, two ordinary propagators $i slb(k)_l \/ k_l^2$, the gauge vertex $i gbar^nu Tg^c$ and the Yukawa factors $(-i y_k)^2$. Reading from $Dbrk_1$ against the flow, the Dirac string is $bold(E) thin gbar^nu slb(k)_1 slb(k)_2$ or $bold(E) thin slb(k)_1 slb(k)_2 gbar^nu$ when the gauge vertex is a neighbour of $Dbrk_1$ (b), and $bold(E) thin slb(k)_1 gbar^nu slb(k)_2$ when it is opposite (c).

(d), (e) The pentagons carry @eq:g1_cSS_cSSV-E and three ordinary propagators.

*$Rop$-operation.* At one loop $Rbar G = G$. $D = 1\/(k^2 - MUV^2)$ denotes the expanded denominators of @eq:uv-expansion. Since @eq:g1_cSS_cSSV-E is of order $q$, every term of $Top(dg(G))$ has one explicit factor $q$ from it; $MUV^2$ and $m^2$ have degree $2$.

_Zero-derivative functions_, $dg(G) = 0$. $Top(0)$ sets every external momentum to zero, so $p = p'$ and @eq:g1_cSS_cSSV-E vanishes identically: $Top(0) G_"(d)" = Top(0) G_"(e)" = 0$, i.e. $N_0 = N_1 = 0$ for every placement.

_Box_, $dg(G) = 1$. $Top(1)$ keeps the part of @eq:g1_cSS_cSSV-E linear in $q$ at $p = p' = k$ and all other momenta zero:
$ Top(1) bold(E) = - Tg^a thin gbar_nu u^nu thin D^2, quad u^nu = bar(q)^nu hat(k)^2 - 2 (hat(k) dot hat(q)) bar(k)^nu, quad quad
  Top(1) G = i/2 thin y_k^2 thin [dots] intk D^4 thin t , $ <eq:g1_cSS_cSSV-T1>
with $[dots]$ the group factor of @tab:g1_cSS_cSSV-diagrams and $t$ the four-dimensional trace of the string (the projector gives $1\/2$; its $g5$ part is proportional to $epsilon^(alpha beta gamma nu) u_alpha bar(k)_beta bar(k)_gamma = 0$). The scalar factor is $(-1)(-1)(-i y_k)^2 (i)(i)^2 = i y_k^2$, the first $-1$ from @eq:g1_cSS_cSSV-E. The traces are
$ t_"(b)" = 4 thin bar(k)^2 u^nu, quad quad t_"(c)" = 4 [ 2 (u dot bar(k)) bar(k)^nu - bar(k)^2 u^nu ] . $
Terms with an odd number of $hat(k)$ vanish, so $hat(q)$ drops out. With @eq:tensor-red and @eq:eta, $bar(k)^2 hat(k)^2 -> 4 (-2 eta) (k^2)^2 \/ 24$ and $bar(k)^alpha bar(k)^nu hat(k)^2 -> (-2 eta) gmb^(alpha nu) (k^2)^2 \/ 24$:
$ t_"(b)" -> - 4/3 eta thin (k^2)^2 bar(q)^nu, quad quad t_"(c)" -> 4 (- 1/6 + 1/3) eta thin (k^2)^2 bar(q)^nu = 2/3 eta thin (k^2)^2 bar(q)^nu , $ <eq:g1_cSS_cSSV-box-t>
and $intk (k^2)^2 D^4 -> ilf \/ ep$. Hence $N_0 = 0$ and
$ N_1^"(b)" = - (2 i)/3 y_k^2 [dots] thin bar(q)^nu, quad quad N_1^"(c)" = i/3 y_k^2 [dots] thin bar(q)^nu . $ <eq:g1_cSS_cSSV-N-box>
Applying @eq:g1_cSS_cSSV-cov twice gives $br({Tg^a, Tg^c}, 1) + br(1, {Tg^a, Tg^c}) = 2 (br(Tg^a, Tg^c) + br(Tg^c, Tg^a)) + N^c_k {Tg^a_H, Tg^c_H}$, so the six boxes of one Yukawa pair sum to
$ N_1^"(b+c)" = i/2 thin y_k^2 thin bar(q)^nu thin C^(a c)_k, quad quad
  C^(a c)_k = - 4/3 N^c_k {Tg^a_H, Tg^c_H} - 2 (br(Tg^a, Tg^c) + br(Tg^c, Tg^a)) . $ <eq:g1_cSS_cSSV-C>

_Triangle_, $dg(G) = 2$. With $v_plus.minus = v_1 plus.minus v_2$ the numerator of @eq:g1_cSS_cSSV-F is, after the quartic terms $plus.minus bar(k)^2 hat(k)^2$ cancel,
$ underbrace(2 bar(k)^2 (hat(k) dot hat(v)_+) - hat(k)^2 (bar(k) dot bar(v)_+), "degree 1") + underbrace(bar(k)^2 (hat(v)_1^2 - hat(v)_2^2) - 2 (bar(k) dot bar(v)_2)(hat(k) dot hat(v)_1) + 2 (bar(k) dot bar(v)_1)(hat(k) dot hat(v)_2), "degree 2") + O(v^3) . $
$Top(2)$ multiplies the first group with $Top(1) [k^2 (k + v_1)^2 (k - v_2)^2]^(-1) = D^3 - 2 (k dot v_-) D^4$ and the second with $D^3$; $MUV^2$ terms would have degree $3$ and do not occur. The reduction gives
$ -2 (k dot v_-) [2 bar(k)^2 (hat(k) dot hat(v)_+) - hat(k)^2 (bar(k) dot bar(v)_+)] -> (k^2)^2/24 [ -16 thin hat(v)_+ dot hat(v)_- - 4 eta thin bar(v)_+ dot bar(v)_- ], \
  [dots]_"degree 2" -> k^2 (hat(v)_1^2 - hat(v)_2^2) , $
(mixed $bar(k) hat(k)$ pairs vanish), and with $hat(v)_+ dot hat(v)_- = hat(v)_1^2 - hat(v)_2^2$, $bar(v)_+ dot bar(v)_- = bar(v)_1^2 - bar(v)_2^2$,
$ Top(2) F(v_1, v_2) = ilf [ 1/ep thin 1/3 (hat(v)_1^2 - hat(v)_2^2) - eta/ep thin 1/6 (bar(v)_1^2 - bar(v)_2^2) ] + dots $ <eq:g1_cSS_cSSV-F-pole>
Per topology, $N_0^"(a1)" = -(2i)/3 y_k^2 br(Tg^a, 1) (hat(p)_1^2 - hat(p)_2^2)$, $N_1^"(a1)" = i/3 y_k^2 br(Tg^a, 1) (bar(p)_1^2 - bar(p)_2^2)$, and (a2) the same with $br(Tg^a, 1) -> -br(1, Tg^a)$. By @eq:g1_cSS_cSSV-a the sum over (a1), (a2) and $k$ is
$ N_0 = -(2 i)/3 Yt (Tg^a_H)_(i j) (hat(p)_1^2 - hat(p)_2^2), quad quad N_1 = i/3 Yt (Tg^a_H)_(i j) (bar(p)_1^2 - bar(p)_2^2) . $ <eq:g1_cSS_cSSV-N-tri>

*Result.* Translation (@sec:m-oneloop): on $gh^a (q)$, $partial gh -> -i q$. $bop_d E_H$ of @eq:b-evanescent has the $gh H^dagger H$ vertex function $i (hat(p)_2^2 - hat(p)_1^2) Tg^a_H$; $N^((r))_S$ has $-i (p_1^2 - p_2^2) Tg^a_H$; $N^((r s))_(S V)$ has $-i q^nu {Tg^a_H, Tg^c_H}$ ($a in r$, $c in s$); $N_(S G)$ has $-i q^nu gs^2 delta^(a c)$.
#result[
  $ (Dbrk_1 dot Gamma)^((1))_(gh^a H^dagger_i H_j) |_"pole" = lf 1/ep N_0, quad
    N_0 = -(2 i)/3 Yt (Tg^a_H)_(i j) (hat(p)_1^2 - hat(p)_2^2) quad "(evanescent)", $ <eq:g1_cSS_cSSV-pole>
  $ (Dbrk_1 dot Gamma)^((1))|_("pole", gh H^dagger H) & = lf (2 Yt)/(3 ep) thin bop_d E_H \
    & = lf (2 Yt)/(3 ep) thin i integral dif^d x [ (hat(partial)_mu H)^dagger hat(partial)^mu (wm H) - hat(partial)_mu (H^dagger wm) hat(partial)^mu H ] , $ <eq:g1_cSS_cSSV-pole-op>
  $ (Dbrk_1 dot Gamma)^((1))_(gh^a H^dagger_i H_j A^c_nu) |_"pole" = 0, quad
    (Dbrk_1 dot Gamma)^((1))_(gh H^dagger H A A) |_"pole" = (Dbrk_1 dot Gamma)^((1))_(gh (H^dagger H)^2) |_"pole" = 0 . $ <eq:g1_cSS_cSSV-pole-box>
  $ & X_(1; gh^a H^dagger_i H_j) = lf thin Pphys N_1, quad Pphys N_1 = i/3 Yt (Tg^a_H)_(i j) (p_1^2 - p_2^2), \
    & X_1 supset - lf Yt/3 (N^((W))_S + N^((B))_S) , $ <eq:g1_cSS_cSSV-X2>
  $ X_(1; gh^a H^dagger_i H_j A^c_nu) = lf thin Pphys N_1, quad Pphys N_1 = i/2 thin q^nu sum_k y_k^2 thin C^(a c)_(k, i j) , $ <eq:g1_cSS_cSSV-X3>
  with $C^(a c)_k$ of @eq:g1_cSS_cSSV-C, the vertex function of
  $ X_1 supset lf sum_k y_k^2 integral dif^4 x [ & 2/3 N^c_k thin H^dagger {partial_mu wm, Am^mu} H \
    & + Tr(Phi_k^dagger (partial_mu wm^(f')) Phi_k Am^(f mu)) + Tr(Phi_k^dagger Am^(f' mu) Phi_k thin partial_mu wm^f) ] , $ <eq:g1_cSS_cSSV-X3-op>
  $ X_1 supset lf [ & 2/3 Yt N^((W W))_(S V) + (y_d^2 - 1/3 y_e^2)(N^((W B))_(S V) + N^((B W))_(S V)) \
    & + (10/3 y_u^2 + 4/3 y_d^2 + 8/3 y_e^2) N^((B B))_(S V) + (y_u^2 + y_d^2) N_(S G) ] , $ <eq:g1_cSS_cSSV-X3-basis>
  $ X_(1; gh H^dagger H A A) = X_(1; gh (H^dagger H)^2) = 0 . $ <eq:g1_cSS_cSSV-X0>
]
In @eq:g1_cSS_cSSV-X3-op $wm$, $Am$ without label are in the representation of $H$, $wm^f = Tg^a_f gh^a$, $Am^f_mu = Tg^a_f A^a_mu$, and $Tr$ runs over the indices of $f$. No $N_(S C)$ (commutator), $N_(S V E)$, $N_(S 2)$, $N_(S E)$, $N_(S m)$ or $N_cal(J)$ structure occurs, and the mixed gluon--electroweak functions vanish ($Tr T^A = 0$). For a pair with $f'$ the doublet ($k = d, e$) the coefficient of $N^((W B))_(S V)$ is $N^c_k (2/3 + Y_f) y_k^2$, for $f$ the doublet ($k = u$) it is $N^c_k (2/3 - Y_(f')) y_k^2$; the coefficient of $N^((B B))_(S V)$ is $N^c_k (2/3 + 4 Y_f Y_(f')) y_k^2$.

#figure(
  table(
    columns: 7,
    align: (left, center, center, center, center, center, center),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else if y == 3 { (bottom: 0.3pt) } else { none },
    table.header([pair], [pole: $bop_d E_H$], [$N^((W))_S, N^((B))_S$], [$N^((W W))_(S V)$], [$N^((W B))_(S V), N^((B W))_(S V)$], [$N^((B B))_(S V)$], [$N_(S G)$]),
    [$d$: $(D, Q)$], [$2 y_d^2$], [$-y_d^2$], [$2 y_d^2$], [$y_d^2$], [$4/3 y_d^2$], [$y_d^2$],
    [$u$: $(Q, U)$], [$2 y_u^2$], [$-y_u^2$], [$2 y_u^2$], [$0$], [$10/3 y_u^2$], [$y_u^2$],
    [$e$: $(E, L)$], [$2/3 y_e^2$], [$-1/3 y_e^2$], [$2/3 y_e^2$], [$-1/3 y_e^2$], [$8/3 y_e^2$], [$0$],
    [sum], [$2/3 Yt$], [$-1/3 Yt$], [$2/3 Yt$], [$y_d^2 - 1/3 y_e^2$], [$10/3 y_u^2 + 4/3 y_d^2 + 8/3 y_e^2$], [$y_u^2 + y_d^2$],
  ),
  caption: [Standard-Model values. Column 2: coefficient of $bop_d E_H$ in $16 pi^2 ep (Dbrk_1 dot Gamma)^((1))|_"pole"$, @eq:g1_cSS_cSSV-pole-op. Columns 3--7: coefficients in $16 pi^2 X_1$, @eq:g1_cSS_cSSV-X2 and @eq:g1_cSS_cSSV-X3-basis. The hypercharges enter through $Y_D = -1/3$, $Y_Q = 1/6$, $Y_U = 2/3$, $Y_E = -1$, $Y_L = -1/2$; the $"SU"(2)_L$ structure is $Tg^I_H = gw tau^I \/ 2$ in $N_S$, $N_(S V)$. The mixed functions with one gluon index and one electroweak index vanish.],
) <tab:g1_cSS_cSSV-sm>

#check[*Cancellation of the pole.* By the check of @sec:m-oneloop, @eq:g1_cSS_cSSV-pole-op must equal $-bop_d delta S^((1))$ restricted to the fermion-loop (Yukawa) part. The Yukawa parts of $16 pi^2 ep thin Gamma^((1))|_"pole"$ for $H^dagger H$, $H^dagger H A$ and $H^dagger H A A$ are $Yt (bar(p)^2 + 2/3 hat(p)^2)$, $Yt (bar(p) + bar(p)')^mu Tg^a_H$ and $Yt thin gmb^(mu nu) {Tg^a_H, Tg^b_H}$ (@sec:g0_SS, @sec:g0_SSV_SSVV), the vertex functions of $Yt [(bar(D)_mu H)^dagger bar(D)^mu H + 2/3 E_H] = Yt [(D_mu H)^dagger D^mu H + tilde(E)_H + 2/3 E_H]$. By $delta S^((1)) = -Kop Gamma^((1))$, the Yukawa part of $delta S^((1))$ in @sec:g0_SS and @sec:g0_SSV_SSVV is $-lf (1\/ep) Yt [(D_mu H)^dagger D^mu H + tilde(E)_H + 2/3 E_H]$. Since $(D H)^dagger D H$ and $tilde(E)_H$ are $bop_d$-closed, $-bop_d delta S^((1)) supset lf (2 Yt)/(3 ep) bop_d E_H$: this is @eq:g1_cSS_cSSV-pole-op, and its vanishing $gh H^dagger H A$ component ($bop_d E_H$ contains no gauge field) is @eq:g1_cSS_cSSV-pole-box. The identity was also verified diagram class by diagram class with the fermion-loop poles of $Gamma^((1))_(H^dagger H)$, $Gamma^((1))_(H^dagger H A)$ and $Gamma^((1))_(H^dagger H A A)$ recomputed from the Feynman rules and $bop_d$ applied leg by leg (@eq:m-b-rules).]

#check[*Literature.* @eq:g1_cSS_cSSV-X2 is $X_1 supset -lf Yt/3 thin bop Q_(S S)$ (@eq:btable-higgs), the image in our conventions of Eq. (6.22) of @Belusca-Maito:2020ala. @eq:g1_cSS_cSSV-X3-op equals the $partial wm$ part of
$ - bop lf sum_k y_k^2 integral dif^4 x [ & - 2/3 ( Tr_(f')(Am^(f') dot Am^(f') Phi_k Phi_k^dagger) + Tr_f (Am^f dot Am^f Phi_k^dagger Phi_k) ) \
  & + 1/3 Tr_f (Phi_k^dagger Am^(f') Phi_k dot Am^f) ] , $ <eq:g1_cSS_cSSV-lit>
the image of the $Am Am Phi Phi$ terms of the finite counterterm of @Belusca-Maito:2020ala (Eqs. (6.23), (6.46)) for a reducible representation; the two forms are related by @eq:g1_cSS_cSSV-cov. The vanishing of the functions with $dg(G) = 0$ agrees with the absence of such terms there.]

#check[*Rigid invariance, absence of $m^2$ and of $MUV$.* At $q = 0$ the $Dbrk_1$ vertex with its propagators vanishes (@eq:g1_cSS_cSSV-E). Consequently every result vanishes at $q = 0$: $p_1^2 - p_2^2 = -q dot (p_1 - p_2)$, $hat(p)_1^2 - hat(p)_2^2$, $q^nu$, and the zero-derivative functions @eq:g1_cSS_cSSV-X0. A term $i m^2 integral H^dagger wm H$ ($N^((r))_(S m)$) would be nonzero at $q = 0$ and is therefore excluded; at one loop it also cannot arise, because no graph contains a Higgs propagator and the fermions are massless, so $m^2$ never enters. $MUV$ does not appear at all: the explicit factor $q$ uses one unit of $dg(G)$, and $MUV^2$ would exceed it.]

#check[*Chirality and group structure.* $g5$ drops out of every trace (two or four barred matrices, $epsilon$ contracted with two equal momenta), so left- and right-handed lines enter alike and only through the generators. $C^(a c)_k$ is symmetric in $(a, c)$: no commutator structure $N_(S C)$ occurs, and the coefficients of $N^((W B))_(S V)$ and $N^((B W))_(S V)$ agree. The identity @eq:g1_cSS_cSSV-E was verified with an explicit $8 times 8$ representation of the Clifford algebra with four barred and two hatted directions for both chiralities; the reductions @eq:g1_cSS_cSSV-box-t and @eq:g1_cSS_cSSV-F-pole, and the group sums of @tab:g1_cSS_cSSV-sm with explicit $"SU"(3) times "SU"(2) times "U"(1)$ matrices, were redone independently; the same reduction reproduces $N_0$ and $N_1$ of @sec:g1_cV. The whole class was also computed by an independent numerical implementation (@sec:supporting): all cyclic orders of the labelled legs in one flavour space of all multiplets, symbolic traces cross-checked with explicit matrices, one triangle redone with Feynman parameters. It reproduces every entry of @tab:g1_cSS_cSSV-sm, the vanishing pole of every box, and the vanishing at $p = 0$ of the numerator of every graph with $dg(G) = 0$.]

#check[*Image of #bop.* $X_1$ of this class lies in $bop("span"{Q_alpha})$. With @eq:btable-higgs,
$ -X_1 |_(gh H^dagger H, thin gh H^dagger H A) = bop [ lf ( & Yt/3 (Q^((W))_J + Q^((B))_J) - Yt thin Q_(W W) - (13/3 y_u^2 + 7/3 y_d^2 + 3 y_e^2) Q_(B B) \
  & - (y_u^2 + 2 y_d^2) Q_(W B) - 1/2 (y_u^2 + y_d^2) Q_(G G) ) ] . $
This combination is the Higgs part of the normalised solution of @sec:st-solution (@eq:sts-r-J, @eq:sts-r-higgs, @tab:sts-r): it contains no $Q_(S S)$ (removed by (N1)), no $Q^((r))_D$ and no $Q^((W))_sigma$. Within these two classes alone it is not unique: $bop Q_(S S) = N^((W))_S + N^((B))_S$ and $bop Q^((W))_sigma supset N^((W))_S - N^((W B))_(S V)$ could replace part of it. $Q_(S S)$ is excluded by (N1), and $r[Q^((W))_sigma] = 0$ is fixed in @sec:st-solution by the class $gh gh srcH^dagger H$ (@eq:sts-r-sigma); it is confirmed here by the absence of $N^((W W))_(S V E)$ components in $X_1$ and by $X_1 = 0$ in the class $gh H^dagger H A A$ (@eq:g1_cSS_cSSV-X0), where $bop Q^((W))_sigma$ has the components $N^((W B W))_(S V V)$. The structure parallels @sec:g1_cV, whose numerator is @eq:g1_cSS_cSSV-E closed with one gauge vertex.]

#general[At one loop the functions with $dg(G) = 0$ vanish because @eq:g1_cSS_cSSV-E is of order $q$. The general statement is rigid invariance: at constant ghost, #sop rotates only the physical components, so $Dbrk_1 |_(q = 0)$ is the variation of $S_0$ under a rotation of the sterile partners alone, a linear change of variables. Its insertion into 1PI functions without external $chi_f$ vanishes, and the counterterms, which contain no $chi_f$ and are rigidly invariant, add nothing at $q = 0$. Hence $Dbrk^([l]) dot Gamma$ of @eq:m-Delta-l vanishes at zero ghost momentum at every order, and $X_l$ has no zero-derivative component in these classes and no $N^((r))_(S m)$ term. At two loops the hatted parts kept here enter through @eq:m-LP: the evanescent pole @eq:g1_cSS_cSSV-pole-op cancels against $bop_d delta S^((1))$, and only $Pphys N_1$ has been used.]
