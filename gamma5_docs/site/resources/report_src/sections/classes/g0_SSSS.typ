#import "../../macros.typ": *
#import "../../diagrams.typ": *

// local shorthand (this file only): the Bose-symmetric SU(2) structure of the quartic vertex
#let Sij = $(delta_(i j) delta_(i' j') + delta_(i j') delta_(i' j))$

=== Higgs four-point function $(H^dagger H)^2$ <sec:g0_SSSS>

External legs $H^dagger_i, H_j, H^dagger_(i'), H_(j')$, all at zero momentum; $i, j, i', j'$ are $"SU"(2)_L$ indices, and $l, n$ denote internal $"SU"(2)_L$ indices (the letter $k$ is the loop momentum, $q$ the momentum of an internal $H$ line). The sum of the diagrams below is $i thin Gamma^((1))_(H^dagger_i H_j H^dagger_(i') H_(j'))$; the tree level is $Gamma^((0))_(H^dagger_i H_j H^dagger_(i') H_(j')) = -2 lambda #Sij$ (@eq:V-ssss divided by $i$, with its indices $k, l$ renamed $i', j'$). Every diagram has superficial degree $dg(G) = 0$, so the pole is a constant and zero external momenta lose nothing. Throughout, $Tg^a equiv Tg^a_H$ (@eq:Tg-explicit), $a$ runs over $W^1, W^2, W^3, B$, and repeated indices are summed.

*Diagrams.* At one loop every vertex has two internal legs. With four external scalars the possible vertices are: the quartic vertex @eq:V-ssss with two external legs, the seagull @eq:V-ssv with two external scalars, the cubic vertex @eq:V-ssv with one external scalar, and the Yukawa vertices @eq:V-yuk with one external scalar. This gives the six topologies of @fig:g0_SSSS-diagrams. Ghosts and gluons do not couple to $H$. A three-gauge vertex in the loop leaves its third leg to end on a cubic vertex with two external scalars, which makes the graph one-particle reducible.

#figure(
  stack(dir: ttb, spacing: 5mm,
  diagrow(
    oneloop(
      ((legs: (("S", $H^dagger_i$, -1), ("S", $H_j$, 1))), (legs: (("S", $H_(j')$, 1), ("S", $H^dagger_(i')$, -1)))),
      (("S", $H$, 1), ("S", $H$, 1)), caption: [(a) $lambda lambda$ bubble (t)]),
    oneloop(
      ((legs: (("S", $H^dagger_i$, -1),)), (legs: (("S", $H_j$, 1),)), (legs: (("S", $H^dagger_(i')$, -1),)), (legs: (("S", $H_(j')$, 1),))),
      (("F", $U$, -1), ("F", $Q$, -1), ("F", $U$, -1), ("F", $Q$, -1)), caption: [(b) fermion box]),
    oneloop(
      ((legs: (("S", $H^dagger_i$, -1), ("S", $H_j$, 1))), (legs: (("S", $H_(j')$, 1), ("S", $H^dagger_(i')$, -1)))),
      (("V", $A^a$, 0), ("V", $A^b$, 0)), caption: [(c) seagull bubble]),
  ),
  diagrow(
    oneloop(
      ((legs: (("S", $H^dagger_i$, -1), ("S", $H_j$, 1))), (legs: (("S", $H_(j')$, 1),)), (legs: (("S", $H^dagger_(i')$, -1),))),
      (("V", $A^b$, 0), ("S", $H$, 1), ("V", $A^a$, 0)), caption: [(d) gauge triangle]),
    oneloop(
      ((legs: (("S", $H_j$, 1),)), (legs: (("S", $H^dagger_i$, -1),)), (legs: (("S", $H^dagger_(i')$, -1),)), (legs: (("S", $H_(j')$, 1),))),
      (("S", $H$, 1), ("V", $A^b$, 0), ("S", $H$, -1), ("V", $A^a$, 0)), caption: [(e) gauge box]),
    oneloop(
      ((legs: (("S", $H^dagger_i$, -1), ("S", $H_j$, 1))), (legs: (("S", $H_(j')$, 1),)), (legs: (("S", $H^dagger_(i')$, -1),))),
      (("S", $H$, -1), ("V", $A^a$, 0), ("S", $H$, -1)), caption: [(f) $lambda$-gauge triangle]),
  )),
  caption: [Topologies of $Gamma^((1))_(H^dagger H H^dagger H)$, each drawn for one field assignment. Arrows on scalar lines give the flow from $H$ to $H^dagger$; (b) is drawn for the $y_u^4$ chain, (e) is the uncrossed box, (f) is case (i) of @tab:g0_SSSS-diagrams.],
) <fig:g0_SSSS-diagrams>

#figure(
  { set par(justify: false); table(
    columns: (auto, 1fr, auto, auto, 1.5fr),
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram $X$], [assignment (number of graphs)], [sym. $s_X$], [sign], [group factor $c_X$]),
    [(a) s], [$H^dagger_i H^dagger_(i')$ at one vertex (1)], [$1\/2$], [$+1$], [$(delta_(i l) delta_(i' n) + delta_(i n) delta_(i' l))(delta_(l j) delta_(n j') + delta_(l j') delta_(n j)) = 2 #Sij$],
    [(a) t], [$H^dagger_i H_j$ at one vertex (1)], [$1$], [$+1$], [$(delta_(i j) delta_(l n) + delta_(i n) delta_(l j))(delta_(i' j') delta_(n l) + delta_(i' l) delta_(n j')) = 4 delta_(i j) delta_(i' j') + delta_(i j') delta_(i' j)$],
    [(a) u], [$H^dagger_i H_(j')$ at one vertex (1)], [$1$], [$+1$], [$4 delta_(i j') delta_(i' j) + delta_(i j) delta_(i' j')$],
    [(b) $y_u^4$], [$U, Q, U, Q$; cyclic orders $(H^dagger_i H_j H^dagger_(i') H_(j'))$, $(H^dagger_i H_(j') H^dagger_(i') H_j)$ (2)], [$1$], [$-1$], [$3 epsilon_(l i) epsilon_(n j) epsilon_(n i') epsilon_(l j') = 3 delta_(i j') delta_(i' j)$, resp. $3 delta_(i j) delta_(i' j')$],
    [(b) $y_d^4$], [$D, Q, D, Q$ (2)], [$1$], [$-1$], [$3 delta_(i j) delta_(i' j')$, resp. $3 delta_(i j') delta_(i' j)$],
    [(b) $y_e^4$], [$E, L, E, L$ (2)], [$1$], [$-1$], [$delta_(i j) delta_(i' j')$, resp. $delta_(i j') delta_(i' j)$],
    [(b) $y_u^2 y_d^2$], [$U, Q, D, Q$; $H^dagger_x, H_y, H_z, H^dagger_w$ with ${x,w} = {i,i'}$, ${y,z} = {j,j'}$ (4)], [$1$], [$-1$], [$3 epsilon_(l x) epsilon_(n y) delta_(n z) delta_(l w) = 3 epsilon_(w x) epsilon_(z y)$],
    [(c)], [$A^a, A^b$; seagulls $(i j)(i' j')$ and $(i j')(i' j)$ (2)], [$1\/2$], [$+1$], [${Tg^a, Tg^b}_(i j) {Tg^a, Tg^b}_(i' j')$ and $j <-> j'$],
    [(d)], [$A^a, A^b, H$; seagull on $(i j)$, $(i' j')$, $(i j')$ or $(i' j)$ (4)], [$1$], [$+1$], [${Tg^a, Tg^b}_(i j) (Tg^a Tg^b)_(i' j')$ and permutations],
    [(e)], [$H, A, H, A$; uncrossed and crossed, two pairings (4)], [$1$], [$+1$], [$(Tg^b Tg^a)_(i j) (Tg^b Tg^a)_(i' j')$, $(Tg^b Tg^a)_(i j) (Tg^a Tg^b)_(i' j')$, and $j <-> j'$],
    [(f) (i)], [$H^dagger_x H_y$ at the quartic vertex, line $H_w -> H^dagger_z$ (4)], [$1$], [$+1$], [$delta_(x y) (Tg^a Tg^a)_(z w) + Tg^a_(z y) Tg^a_(x w)$],
    [(f) (ii)], [$H^dagger_i H^dagger_(i')$ at the quartic vertex (1)], [$1$], [$+1$], [$Tg^a_(i j) Tg^a_(i' j') + Tg^a_(i' j) Tg^a_(i j')$],
    [(f) (iii)], [$H_j H_(j')$ at the quartic vertex (1)], [$1$], [$+1$], [$Tg^a_(i j) Tg^a_(i' j') + Tg^a_(i j') Tg^a_(i' j)$],
  ) },
  caption: [Field assignments of @fig:g0_SSSS-diagrams: 29 graphs. Here an internal gauge line is summed over the gauge factors $W$, $B$ inside the group factor $c_X$. @tab:classes counts 63 diagrams because there every assignment of $W$ or $B$ to the internal gauge lines is a separate diagram: (c) $2 -> 6$ ($W W$, $W B$, $B B$ per pairing), (d) $4 -> 16$, (e) $4 -> 16$, (f) $6 -> 12$, while (a) ($3$) and (b) ($10$) are unchanged. In (a) and (b) $l, n$ are the $"SU"(2)_L$ indices of the internal lines (in (b) of the two $Q$ propagators), the factor $3$ is the colour trace $N^c_Q = 3$, and the chains are read in Dirac order starting with the vertex $overline(Q) U H^dagger$ (resp. $overline(Q) D H$, $overline(L) E H$). In (f) the four assignments of (i) are $(x y z w) = (i j i' j'), (i' j' i j), (i j' i' j), (i' j i j')$. $Tg^a_(i j)$ is the $(i j)$ element of $Tg^a$.],
) <tab:g0_SSSS-diagrams>

*Integrands.* All external momenta vanish, so every internal line carries $plus.minus k$. In (a) and (c) $k$ flows along both lines; in (b) along the fermion arrow; in (d) along the internal $H$ (from the vertex of $H_(j')$ to that of $H^dagger_(i')$) and back through the two gauge lines; in (e) along the $H$ line from $H_j$ to $H^dagger_i$, so that the other $H$ line carries $-k$ (uncrossed) or $+k$ (crossed); in (f) from the quartic vertex through one $H$ line, the gauge line and the other $H$ line. The cubic vertex @eq:V-ssv with one external leg at zero momentum is $i thin q^mu Tg^a$, $q$ the momentum of the internal $H$ line measured along its arrow. With $c_X$ and $s_X$ the group factor and the symmetry factor of @tab:g0_SSSS-diagrams:
$ G_"(a)" = s_X thin c_X thin intk (-2 i lambda)^2 ((i)/(k^2 - m^2))^2 , $ <eq:g0_SSSS-int-a>
$ G_"(b)" = - c_X thin intk tr[ (-i y PR) (i sl(k))/(k^2) (-i y PL) (i sl(k))/(k^2) (-i y PR) (i sl(k))/(k^2) (-i y PL) (i sl(k))/(k^2) ] , $ <eq:g0_SSSS-int-b>
with $y^4 -> y_u^4, y_d^4, y_e^4, y_u^2 y_d^2$; the projectors follow from @eq:V-yuk ($PR$ at $overline(Q) U H^dagger$, $overline(Q) D H$, $overline(L) E H$, $PL$ at their conjugates), so all chains have the same Dirac structure.
$ G_"(c)" = 1/2 intk (i g^(mu nu) {Tg^a, Tg^b}_(i j)) (i g^(alpha beta) {Tg^a, Tg^b}_(i' j')) (-i g_(mu alpha))/(k^2) (-i g_(nu beta))/(k^2) , $ <eq:g0_SSSS-int-c>
$ G_"(d)" = intk (i g^(mu nu) {Tg^a, Tg^b}_(i j)) (i k^alpha Tg^a_(i' l)) (i)/(k^2 - m^2) (i k^beta Tg^b_(l j')) (-i g_(mu alpha))/(k^2) (-i g_(nu beta))/(k^2) , $ <eq:g0_SSSS-int-d>
$ G_"(e)" = intk (i k^nu Tg^b_(i l)) (i)/(k^2 - m^2) (i k^mu Tg^a_(l j)) thin (i k'_nu Tg^b_(i' n)) (i)/(k'^2 - m^2) (i k'_mu Tg^a_(n j')) thin (-i)^2/(k^2 k^2) , quad k' = -k , $ <eq:g0_SSSS-int-e>
and for the crossed box $a <-> b$ on the second line and $k' = k$. In (f), $theta_1, theta_2 = plus.minus 1$ according to whether the routing runs along or against the arrow of each $H$ line. In case (i) the charge flows through the quartic vertex, so both lines have the same orientation relative to the routing; in cases (ii) and (iii) both lines point into (resp. out of) the quartic vertex:
$ G_"(f)" = intk (-2 i lambda) [delta delta] thin (i theta_1 k^mu Tg^a) (i theta_2 k^nu Tg^a) ((i)/(k^2 - m^2))^2 (-i g_(mu nu))/(k^2) , quad theta_1 theta_2 = cases(+1 & "(i)", -1 & "(ii), (iii)") , $ <eq:g0_SSSS-int-f>
where $[delta delta]$ and the $Tg^a$ indices combine into the group factor $c_X$ of @tab:g0_SSSS-diagrams; for (i) with $(x y z w) = (i j i' j')$ it is $(delta_(i j) delta_(l n) + delta_(i n) delta_(l j)) Tg^a_(i' l) Tg^a_(n j')$.

*$Rop$-operation.* At one loop $Rbar G = G$. Power counting: $dg(G) = 4 - 4 = 0$ for (a) and (c), $4 + 4 - 8 = 0$ for (b) and (e), $4 + 2 - 6 = 0$ for (d) and (f). With $n$ propagators the numerators are homogeneous of degree $2 n - 4$ in $k$, so $Top(0)$ keeps only the leading term of @eq:uv-expansion: every denominator becomes $1\/(k^2 - MUV^2)$; $m^2$ and $MUV^2$ in the numerators of @eq:uv-expansion have degree two and drop out. Every integral is then of the form @eq:logpole, $Kop intk (k^2)^a \/ (k^2 - MUV^2)^(a+2) = ilf thin 1\/ep$.

_Numerators._ (a): $(-2 i lambda)^2 i^2 = 4 lambda^2$. (b): $(-1)(-i)^4 i^4 = -1$ and, by @eq:sandwich, $PR sl(k) PL = slb(k) PL$, $PL sl(k) PR = slb(k) PR$, $slb(k) slb(k) = bar(k)^2$:
$ tr[PR slb(k) PL slb(k) PR slb(k) PL slb(k)] = (bar(k)^2)^2 tr PR = 2 (bar(k)^2)^2 , $ <eq:g0_SSSS-trace>
the $g5$ part being $(bar(k)^2)^2 tr g5 = 0$. The four-momentum rule of @eq:tensor-red-split gives $(bar(k)^2)^2 -> 24 (k^2)^2 \/ (d (d + 2)) = (k^2)^2 (1 + O(ep))$, with $gmb^mu_(space mu) = 4$ and no hatted object. (c): $(1\/2) i^2 (-i)^2 g^(mu nu) g^(alpha beta) g_(mu alpha) g_(nu beta) = d\/2 = 2 + O(ep)$. (d): $i^4 (-i)^2 g^(mu nu) g_(mu alpha) g_(nu beta) k^alpha k^beta = -k^2$. (e): $i^6 (-i)^2 (k^2)^2 = (k^2)^2$ for both boxes ($k'_mu k'_nu = k_mu k_nu$). (f): $(-2 i lambda) i^4 (-i) theta_1 theta_2 k^2 = -2 lambda theta_1 theta_2 k^2$. The factors of $d$ are unsplit and enter only the finite part; the residues are
$ N_0^"(a)" = 4 lambda^2 s_X c_X, quad N_0^"(b)" = -2 y^4 c_X, quad N_0^"(c)" = 2 c_X, \
  N_0^"(d)" = -c_X, quad N_0^"(e)" = c_X, quad N_0^"(f)" = -2 lambda theta_1 theta_2 c_X . $ <eq:g0_SSSS-N-graph>
No contraction of hatted objects occurs, the external momenta vanish, and $N_1 = 0$ for every graph.

_Group algebra._ With $Tg^I = gw tau^I\/2$, $Tg^B = gy\/2$ and $tau^I_(i j) tau^I_(i' j') = 2 delta_(i j') delta_(i' j) - delta_(i j) delta_(i' j')$,
$ Tg^a Tg^a = C_H bb(1), quad
  Tg^a_(i j) Tg^a_(i' j') = gw^2/2 delta_(i j') delta_(i' j) + (gy^2 - gw^2)/4 delta_(i j) delta_(i' j'), \
  {Tg^a, Tg^b}_(i j) {Tg^a, Tg^b}_(i' j') = (3/4 gw^4 - 1/2 gw^2 gy^2 + 1/4 gy^4) delta_(i j) delta_(i' j') + gw^2 gy^2 delta_(i j') delta_(i' j) , $ <eq:g0_SSSS-group>
with $C_H = 3/4 gw^2 + 1/4 gy^2$ (@tab:invariants); the three terms of the last line come from $(a, b) in (W, W)$, $(W, B) + (B, W)$, $(B, B)$. Since ${Tg^a, Tg^b}$ is symmetric in $a b$, $(Tg^a Tg^b)_(i' j')$ and $(Tg^b Tg^a)_(i' j')$ may be replaced by $1/2 {Tg^a, Tg^b}_(i' j')$ in (d) and (e). Two index structures of this subsection collect the sums over pairings,
$ cal(U)^((1))_(i j i' j') = Tg^a_(i j) Tg^a_(i' j') + Tg^a_(i j') Tg^a_(i' j) = (gw^2 + gy^2)/4 #Sij, \
  cal(U)^((2))_(i j i' j') = {Tg^a, Tg^b}_(i j) {Tg^a, Tg^b}_(i' j') + {Tg^a, Tg^b}_(i j') {Tg^a, Tg^b}_(i' j) = 1/4 (3 gw^4 + 2 gw^2 gy^2 + gy^4) #Sij , $ <eq:g0_SSSS-Omega>
both summed over $a$ (and $b$). Summing the assignments of @tab:g0_SSSS-diagrams:
$ N_0^"(a)" = 4 lambda^2 [ #Sij + (5 delta_(i j) delta_(i' j') + 5 delta_(i j') delta_(i' j)) ] = 24 lambda^2 #Sij , $ <eq:g0_SSSS-N-a>
$ N_0^"(b)" & = -2 (3 y_u^4 + 3 y_d^4 + y_e^4) #Sij - 6 y_u^2 y_d^2 (epsilon_(i' i) + epsilon_(i i')) (epsilon_(j' j) + epsilon_(j j')) \
  & = -2 (3 y_u^4 + 3 y_d^4 + y_e^4) #Sij , $ <eq:g0_SSSS-N-b>
$ N_0^"(c)" = 2 thin cal(U)^((2))_(i j i' j'), quad N_0^"(d)" = - cal(U)^((2))_(i j i' j'), quad N_0^"(e)" = 1/2 thin cal(U)^((2))_(i j i' j') , $ <eq:g0_SSSS-N-cde>
$ N_0^"(f)" & = -2 lambda [2 C_H #Sij + 2 cal(U)^((1))_(i j i' j')] + 2 lambda [2 cal(U)^((1))_(i j i' j')] \
  & = -4 lambda C_H #Sij = -lambda (3 gw^2 + gy^2) #Sij . $ <eq:g0_SSSS-N-f>
In (a) the s channel gives $s_X c_X = #Sij$ and the channels t, u give $5 #Sij$ together. In (c) the two seagull pairings give $2 c_X$ each, i.e. $2 cal(U)^((2))$. In (d) the replacement $(Tg^a Tg^b)_(i' j') -> 1/2 {Tg^a, Tg^b}_(i' j')$ makes the four seagull positions add to $-1/2 cal(U)^((2)) - 1/2 cal(U)^((2))$. In (e), for each pairing the uncrossed and the crossed box add to $sum_(a b) (Tg^b Tg^a)_(i j) {Tg^a, Tg^b}_(i' j') = 1/2 {Tg^a, Tg^b}_(i j) {Tg^a, Tg^b}_(i' j')$, so the two pairings give $1/2 cal(U)^((2))$. In (f) the first bracket is case (i), the second cases (ii) and (iii); the $cal(U)^((1))$ terms cancel and only the Casimir of the line passing through the quartic vertex survives. The individual mixed boxes give $plus.minus 6 y_u^2 y_d^2 epsilon_(i i') epsilon_(j j') = plus.minus 6 y_u^2 y_d^2 (delta_(i j) delta_(i' j') - delta_(i j') delta_(i' j))$, antisymmetric under $i <-> i'$; the four assignments cancel in pairs. The single pairings in (c)--(e) and the channels t, u in (a) are not proportional to $#Sij$; only their sums are.

*Result.*
#result[
  $ Gamma^((1))_(H^dagger_i H_j H^dagger_(i') H_(j')) |_"pole" = lf 1/ep N_0, \
    N_0 = [ 24 lambda^2 - 2 (3 y_u^4 + 3 y_d^4 + y_e^4) - 4 lambda C_H ] #Sij + 3/2 thin cal(U)^((2))_(i j i' j') , $ <eq:g0_SSSS-pole>
  with $Tg^a Tg^a = C_H bb(1)$ and $cal(U)^((2))_(i j i' j') = sum_(a b) ({Tg^a, Tg^b}_(i j) {Tg^a, Tg^b}_(i' j') + {Tg^a, Tg^b}_(i j') {Tg^a, Tg^b}_(i' j))$ of @eq:g0_SSSS-Omega. Explicitly,
  $ N_0 = [ 24 lambda^2 - 2 (3 y_u^4 + 3 y_d^4 + y_e^4) - lambda (3 gw^2 + gy^2) + 3/8 (3 gw^4 + 2 gw^2 gy^2 + gy^4) ] \
    times #Sij . $ <eq:g0_SSSS-pole-sm>
   The counterterm vertex function is $-lf N_0 \/ ep$. Since $Q_lambda = integral (H^dagger H)^2$ has the vertex function $2 #Sij$,
  $ delta S^((1)) supset - lf 1/ep thin c_lambda integral dif^d x thin (H^dagger H)^2 = - lf 1/ep thin c_lambda thin Q_lambda , $ <eq:g0_SSSS-ct>
  $ c_lambda = 12 lambda^2 - (3 y_u^4 + 3 y_d^4 + y_e^4) - lambda/2 (3 gw^2 + gy^2) + 3/16 (3 gw^4 + 2 gw^2 gy^2 + gy^4) . $ <eq:g0_SSSS-clambda>
  $c_lambda$ (a symbol of this subsection) is $c_lambda = 12 lambda^2 - Y_4 - 2 lambda C_H + 3/16 G_4$ with $Y_4$, $G_4$ of @sec:st-solution, and $-lf c_lambda \/ ep = -delta lambda - 2 lambda delta Z_H$ is the coefficient of $Q_lambda$ in the assembled counterterm. $Q_lambda$ carries no Lorentz index, so its #box[$d$-dimensional] continuation is unique. There is no evanescent operator and no $m^2$ dependence.
]

#figure(
  table(
    columns: 4,
    align: (left, left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([coupling structure], [topologies], [$N_0 \/ #Sij$], [contribution to $c_lambda$]),
    [$lambda^2$], [(a): s $4$, t $+$ u $20$], [$24 lambda^2$], [$12 lambda^2$],
    [$y_u^4$, $y_d^4$, $y_e^4$], [(b)], [$-6 y_u^4 - 6 y_d^4 - 2 y_e^4$], [$-3 y_u^4 - 3 y_d^4 - y_e^4$],
    [$y_u^2 y_d^2$], [(b), four graphs], [$0$], [$0$],
    [$gw^4$], [(c) $3\/2$, (d) $-3\/4$, (e) $3\/8$], [$9\/8 thin gw^4$], [$9\/16 thin gw^4$],
    [$gw^2 gy^2$], [(c) $1$, (d) $-1\/2$, (e) $1\/4$], [$3\/4 thin gw^2 gy^2$], [$3\/8 thin gw^2 gy^2$],
    [$gy^4$], [(c) $1\/2$, (d) $-1\/4$, (e) $1\/8$], [$3\/8 thin gy^4$], [$3\/16 thin gy^4$],
    [$lambda gw^2$], [(f): (i) $-4$, (ii) $+$ (iii) $+1$], [$-3 lambda gw^2$], [$-3\/2 thin lambda gw^2$],
    [$lambda gy^2$], [(f): (i) $-2$, (ii) $+$ (iii) $+1$], [$-lambda gy^2$], [$-1\/2 thin lambda gy^2$],
  ),
  caption: [Standard-Model values of @eq:g0_SSSS-pole and @eq:g0_SSSS-clambda (one family, Feynman gauge). Colour enters only through the multiplicity $N^c_Q = 3$ in the quark boxes; $gs$ does not appear.],
) <tab:g0_SSSS-values>

#check[*Auxiliary mass, Higgs mass, Bose symmetry.* The integrands are logarithmic, $Top(0)$ keeps only the leading term, and $N_0$ contains neither $MUV$ nor $m^2$. Rigid $"SU"(2)_L times "U"(1)_Y$ invariance allows $delta_(i j) delta_(i' j')$ and $delta_(i j') delta_(i' j)$; Bose symmetry at zero momenta ($i <-> i'$) leaves $#Sij$ only. The antisymmetric structure $epsilon_(i i') epsilon_(j j')$ of the single $y_u^2 y_d^2$ boxes cancels, in agreement with $tilde(H)^dagger H = 0$: there is no invariant quartic other than $(H^dagger H)^2$.]

#check[*BMHV.* By @eq:phys-prop the fermion numerators are four-dimensional and the $g5$ trace vanishes (@eq:g0_SSSS-trace); the bosonic numerators are #box[$d$-dimensional] and unsplit. The result therefore coincides with that of an anticommuting $g5$, and $N_1 = 0$. Since $dg(G) = 0$ the pole is momentum independent, so no hatted external momentum can occur at any external kinematics: the class has no evanescent counterterm.]

#check[*Beta function.* With $lambda_0 = mu^(2 ep) (lambda + delta lambda)$ and $H_0 = Z_H^(1\/2) H$ (@sec:st-solution), $-lambda_0 (H_0^dagger H_0)^2$ gives the coefficient $-delta lambda - 2 lambda delta Z_H$ of $Q_lambda$ in $delta S^((1))$. With @eq:g0_SSSS-ct and $16 pi^2 ep thin delta Z_H = 2 C_H - tilde(Y)$ of @sec:g0_SS ($tilde(Y) = 3 y_u^2 + 3 y_d^2 + y_e^2$), $16 pi^2 ep thin delta lambda = c_lambda + 2 lambda (tilde(Y) - 2 C_H)$, the value of @sec:st-solution, and
$ 16 pi^2 beta_lambda & = 2 (c_lambda + 2 lambda (tilde(Y) - 2 C_H)) \
  & = 24 lambda^2 - 3 lambda (3 gw^2 + gy^2) + 4 lambda (3 y_u^2 + 3 y_d^2 + y_e^2) \
  & quad - 2 (3 y_u^4 + 3 y_d^4 + y_e^4) + 3/8 (2 gw^4 + (gw^2 + gy^2)^2) , $
the standard one-loop $beta_lambda$ for $V supset lambda (H^dagger H)^2$ and one family. The $lambda$-gauge term $-3 lambda (3 gw^2 + gy^2) = -12 lambda C_H$ is gauge independent only after adding $Z_H$: the vertex part gives $-4 lambda C_H$, the field renormalisation $-8 lambda C_H$.]

#check[*Slavnov--Taylor operator.* $bop Q_lambda = 0$ (@eq:btable-higgs), so @eq:g0_SSSS-ct produces no ghost-number-1 vertex and does not enter the pole checks of @sec:g1. $Q_lambda$ is a closed direction of @eq:closed on its own (the renormalisation of $lambda$), and it also enters $cal(O)_H$ with the weight $-4 lambda$. The Slavnov--Taylor identity therefore puts no constraint on $c_lambda$, and $Q_lambda$ is excluded from $Delta S^((1))_"ST"$ by the normalisation condition (N1) of @sec:normalisation; the solution is in @sec:st-solution.]

#general[At higher orders the class still has superficial degree $0$ and a single structure, because $(H^dagger H)^2$ is the only rigidly invariant, Bose-symmetric scalar quartic of dimension four, and it has no evanescent partner. The BMHV scheme enters in three ways: the four-dimensional fermion numerators against #box[$d$-dimensional] denominators (the $O(ep)$ terms of $(bar(k)^2)^2 -> 24 (k^2)^2 \/ (d (d+2))$ now multiply subdivergences); $Rbar$ inserts the one-loop counterterms, including the evanescent $E_H$ and $E^((r))_A$ and the finite $Delta S^((1))$, into the two-loop graphs; and graphs with two fermion loops and two $g5$ traces give products $epsilon epsilon$ that are reduced with @eq:epseps. The coefficient of $Q_lambda$ is then read off as above.]
