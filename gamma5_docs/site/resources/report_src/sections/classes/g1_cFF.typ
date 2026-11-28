#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Breaking of $omega overline(Psi) Psi$ <sec:g1_cFF>

External legs $gh^a (q)$, $overline(Psi)_f$ (outgoing momentum $p' = p + q$, i.e. incoming $-p'$) and $Psi_f (p)$, as in @eq:V-delta, for every multiplet $f in {Q, U, D, L, E}$ and every ghost factor $r = G, W, B$ with $Tg^a_f != 0$; the sum of the diagrams below is $i thin (Dbrk_1 dot Gamma)^((1))_(gh^a overline(Psi)_f Psi_f)$, a matrix in spinor and multiplet space with the ghost as the leftmost Grassmann factor. External momenta and spinor structures are $d$-dimensional and $Psi_f$ is not projected, so the pole keeps the sterile components of the legs at a $Dbrk_1$ vertex; $Pphys$ ($hat(p) = hat(p)' = 0$) removes them. Every graph has three propagators and a numerator of degree three, so its superficial degree of divergence is $dg(G) = 4 + 3 - 2 dot 3 = 1$.

*Diagrams.* The ghost leg either sits on the $Dbrk_1$ vertex, which is then internal to the fermion line (topologies (a), (b)), or it enters at the ghost--gauge vertex @eq:V-ghost; the antighost there is contracted with the ghost of a $Dbrk_1$ vertex placed at an external fermion leg (topology (c), two graphs). There is no $agh gh A A$, $agh gh H$ or $overline(Psi) Psi A A$ vertex, so these are all one-loop 1PI graphs.

#figure(
  diagrow(
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("F", $overline(Psi)_f$, -1),)), (legs: (("F", $Psi_f (p)$, 1),))),
      (("F", $k + q$, 1), ("V", $A^e$, 0), ("F", $k$, 1)), caption: [(a)]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("F", $overline(Psi)_f$, -1),)), (legs: (("F", $Psi_f (p)$, 1),))),
      (("F", $Psi_(f')$, 1), ("S", $H$, 0), ("F", $Psi_(f')$, 1)), caption: [(b)]),
    oneloop(
      ((legs: (("G", $gh^a (q)$, 1),)), (legs: (("F", $overline(Psi)_f$, -1),)), (legs: (("F", $Psi_f (p)$, 1),), mark: "delta")),
      (("V", $A^e$, 0), ("F", $k$, -1), ("G", $gh^b$, -1)), caption: [(c1)]),
    oneloop(
      ((legs: (("G", $gh^a (q)$, 1),)), (legs: (("F", $overline(Psi)_f$, -1),), mark: "delta"), (legs: (("F", $Psi_f (p)$, 1),))),
      (("G", $gh^b$, 1), ("F", $k$, -1), ("V", $A^e$, 0)), caption: [(c2)]),
  ),
  caption: [Topologies of $(Dbrk_1 dot Gamma)^((1))_(gh overline(Psi) Psi)$. Crossed circle: the vertex @eq:V-delta. In (a), (b) the fermion line enters the $Dbrk_1$ vertex with momentum $k$ and leaves it with $k + q$; in (b) the internal fermion is a Yukawa partner $Psi_(f')$. In (c1) [(c2)] the $Dbrk_1$ vertex sits at the incoming [outgoing] fermion leg and the internal fermion carries $k$.],
) <fig:g1_cFF-diagrams>

#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto),
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [symmetry factor], [sign], [group factor]),
    [(a)], [$Psi_f$, $Psi_f$, $A^e$ ($e$ over all 12)], [$1$], [$+1$], [$sum_e Tg^e_f Tg^a_f Tg^e_f = (C_f - 1/2 C_"ad"^((a))) Tg^a_f$],
    [(b)], [$Psi_(f')$, $Psi_(f')$, $H$; $f' = U, D$ for $Q$; $Q$ for $U, D$; $E$ for $L$; $L$ for $E$], [$1$], [$+1$], [$cal(M)^a_f = sum_(f', j) bold(Y)^j_(f f') Tg^a_(f') (bold(Y)^j_(f f'))^dagger$],
    [(c1)], [$Psi_f$, $A^e$, $gh^b$; $Dbrk_1$ at the $Psi_f$ leg], [$1$], [$+1$], [$fg^(b e a) Tg^e_f Tg^b_f = -i/2 C_"ad"^((a)) Tg^a_f$],
    [(c2)], [$Psi_f$, $A^e$, $gh^b$; $Dbrk_1$ at the $overline(Psi)_f$ leg], [$1$], [$+1$], [$fg^(b e a) Tg^b_f Tg^e_f = +i/2 C_"ad"^((a)) Tg^a_f$],
  ),
  caption: [Field assignments of @fig:g1_cFF-diagrams. No closed loop occurs. (c1), (c2) exist only for $a$ in $"SU"(3)_c$ or $"SU"(2)_L$, since $fg^(b e a) = 0$ for $a = B$. In (b) the $Dbrk_1$ vertex carries the generator of the partner $f'$. $bold(Y)^j_(f f')$ is the coupling matrix $bold(Y)^j$ of @sec:g0_FFS for the vertex $overline(Psi)_f Psi_(f') phi_j$ of @eq:V-yuk, i.e. that vertex without the factor $-i$ and the projector ($phi_j in {H_j, H^dagger_j}$), so that the vertex absorbing $Psi_f$ is $-i (bold(Y)^j_(f f'))^dagger P_(c_f)$; $cal(M)^a_f$ is defined here for this subsection only.],
) <tab:g1_cFF-diagrams>

The sign of (c) follows from Wick's theorem on $[gh^b overline(Psi)_1 M Psi_1] [(partial agh^e) fg A gh^a] [overline(Psi)_2 V Psi_2]$: the contraction $gh^b agh$ is adjacent up to the even bilinear $overline(Psi)_1 M Psi_1$ and gives the propagator $i delta \/ ell^2$ of @eq:propagators; the external $gh^a$ moves to the left past even bilinears only; the fermion contraction $Psi_2 overline(Psi)_1$ (c1) or $Psi_1 overline(Psi)_2$ (c2) is that of an open line. The result has the form $gh^a overline(Psi) [dots] Psi$ with sign $+1$, as for (a) and (b).

*Integrands.* Write $Pc = P_(c_f)$, $Pcb = P_(overline(c)_f)$. In (a) and (b) the boson carries $p - k$ from the vertex that absorbs $Psi_f (p)$ to the vertex that emits $overline(Psi)_f$; in @eq:V-delta the fermion momenta are $k$ (in) and $k + q$ (out):
$ G_"(a)" = intk (i gbar^mu Pc Tg^e_f) (i (sl(k) + sl(q)))/((k+q)^2) Tg^a_f (slh(k) Pcb - (slh(k) + slh(q)) Pc) (i sl(k))/(k^2) (i gbar^nu Pc Tg^e_f) (-i g_(mu nu))/((k - p)^2) , $ <eq:g1_cFF-int-a>
$ G_"(b)" = sum_(f', j) intk & (-i bold(Y)^j_(f f') P_(c_(f'))) (i (sl(k) + sl(q)))/((k+q)^2) Tg^a_(f') (slh(k) P_(overline(c)_(f')) - (slh(k) + slh(q)) P_(c_(f'))) \
  & times (i sl(k))/(k^2) (-i (bold(Y)^j_(f f'))^dagger Pc) (i)/((k - p)^2 - m^2) , $ <eq:g1_cFF-int-b>
with $P_(c_(f')) = Pcb$. In (c1) the internal fermion leaves the $Dbrk_1$ vertex with $k$, the ghost arrives there with $ell = k - p$ (outgoing at the antighost leg of @eq:V-ghost), and $A^e$ carries $p' - k$ to the gauge vertex; in (c2) the internal fermion enters the $Dbrk_1$ vertex with $k$, the ghost arrives with $ell = p' - k$ and $A^e$ carries $k - p$:
$ G_"(c1)" = intk (i gbar^nu Pc Tg^e_f) (i sl(k))/(k^2) Tg^b_f (slh(p) Pcb - slh(k) Pc) (i)/((k-p)^2) (- fg^(b e a) (k - p)^mu) (-i g_(mu nu))/((k - p')^2) , $ <eq:g1_cFF-int-c1>
$ G_"(c2)" = intk Tg^b_f (slh(k) Pcb - slh(p)' Pc) (i sl(k))/(k^2) (i gbar^nu Pc Tg^e_f) (i)/((p' - k)^2) (- fg^(b e a) (p' - k)^mu) (-i g_(mu nu))/((k - p)^2) . $ <eq:g1_cFF-int-c2>
The numerical prefactors are $-i$ in (a), $+i$ in (b) and $+1$ in (c1), (c2).

*$Rop$-operation.* At one loop $Rbar G = G$. _Dirac algebra._ With @eq:sandwich, $Pc sl(v) = slb(v) Pcb + slh(v) Pc$ and $Pcb slh(v) Pc = 0$; the sterile insertions are those of @eq:m-sterile. For (a),
$ Pc (sl(k) + sl(q)) (slh(k) Pcb - (slh(k) + slh(q)) Pc) sl(k) gbar_mu Pc = [hat(k)^2 (slb(k) + slb(q)) - (hat(k) + hat(q))^2 slb(k)] gbar_mu Pc , $
the first term from the sterile line entering $Dbrk_1$, the second from the sterile line leaving it. With $gbar^mu slb(v) gbar_mu = -2 slb(v)$ the terms cubic in $k$ cancel, and (b) is the same string without $gbar^mu dots gbar_mu$ and with $f'$ in place of $f$ (the projector on the right is again $Pc$):
$ G_"(a)" = 2 i thin (C_f - 1/2 C_"ad"^((a))) Tg^a_f intk (hat(k)^2 slb(q) - 2 (hat(k) dot hat(q)) slb(k) - hat(q)^2 slb(k)) / (k^2 (k + q)^2 (k - p)^2) Pc , $ <eq:g1_cFF-G-a>
$ G_"(b)" = i thin cal(M)^a_f intk (hat(k)^2 slb(q) - 2 (hat(k) dot hat(q)) slb(k) - hat(q)^2 slb(k)) / (k^2 (k + q)^2 ((k - p)^2 - m^2)) Pc . $ <eq:g1_cFF-G-b>
For (c), $gbar^nu (k-p)_nu = (slb(k) - slb(p))$ and $Pc sl(k) slh(p) Pcb = slb(k) slh(p) Pcb$, $Pc sl(k) slh(k) Pc = hat(k)^2 Pc$; similarly in (c2):
$ G_"(c1)" = - i/2 C_"ad"^((a)) Tg^a_f intk ((slb(k) - slb(p)) slb(k) slh(p) Pcb - hat(k)^2 (slb(k) - slb(p)) Pc) / (k^2 (k - p)^2 (k - p')^2) , $ <eq:g1_cFF-G-c1>
$ G_"(c2)" = i/2 C_"ad"^((a)) Tg^a_f intk (hat(k)^2 (slb(p)' - slb(k)) Pc - slh(p)' slb(k) (slb(p)' - slb(k)) Pc) / (k^2 (k - p')^2 (k - p)^2) . $ <eq:g1_cFF-G-c2>
In (c) the terms with $Pcb$ on the right [with $slh(p)'$ on the left] have the sterile component as the external $Psi_f$ [$overline(Psi)_f$] leg and a physical internal line; the terms with $hat(k)^2$ have physical external legs and the sterile component on the internal line.

_Expansion._ $Top(1)$ keeps the terms of degree $<= 1$ in $p$, $p'$, $q$, $m$, $MUV$. With $D = 1 \/ (k^2 - MUV^2)$ and @eq:uv-expansion, $1\/((k - v)^2) = D + 2 (k dot v) D^2 + O(2)$; $MUV^2$ and $m^2$ are of degree two and never enter. The integrals needed are
$ hat(k)^2 D^3 -> ilf (-2 eta \/ d) / ep -> ilf (- eta / 2) 1/ep, quad quad
  bar(k)^2 D^3 -> ilf (4 \/ d) / ep -> ilf 1/ep, \
  hat(k)^2 slb(k) thin 2 (k dot P) thin D^4 -> ilf (-(4 eta)/(d (d + 2))) 1/ep slb(P) -> ilf (- eta/6) 1/ep slb(P) , $ <eq:g1_cFF-ints>
by @eq:tensor-red, @eq:tensor-red-split and @eq:logpole: in the four-index reduction only the pairing $g^(alpha beta) g^(gamma delta)$ with $gmh_(alpha beta)$ survives, $gmh_(alpha beta) g^(alpha beta) = -2 eta$, and the other two pairings contract $gmh$ with $gbar$. In (a) and (b) every term of the numerator carries an explicit external factor, so only $D^3$ is needed: $hat(k)^2 slb(q)$ gives the tag, $(hat(k) dot hat(q)) slb(k)$ vanishes ($hat(k)^alpha bar(k)^beta -> 0$), and $hat(q)^2 slb(k)$ is of degree $2$. In (c1), $(slb(k) - slb(p)) slb(k) slh(p) -> bar(k)^2 slh(p)$ and $-hat(k)^2 (slb(k) - slb(p)) = -hat(k)^2 slb(k) + hat(k)^2 slb(p)$ with $P = p + p'$ from the two denominators; (c2) is the mirror image. Hence
$ Top(1) G_"(c1)" & = - i/2 C_"ad"^((a)) Tg^a_f thin ilf 1/ep [ slh(p) Pcb + eta (-1/3 slb(p) + 1/6 slb(p)') Pc ] + dots , \
  Top(1) G_"(c2)" & = i/2 C_"ad"^((a)) Tg^a_f thin ilf 1/ep [ slh(p)' Pc + eta (1/6 slb(p) - 1/3 slb(p)') Pc ] + dots . $
The same coefficients follow from Feynman parameters: after the shift $k = l + x_2 v_2 + x_3 v_3$, only $hat(l)^2$ survives $Pphys$, and the average of each $x_i$ over the simplex is $1\/3$ (e.g. (c1): $(slb(k) - slb(p)) -> -2/3 slb(p) + 1/3 slb(p)'$, times $-hat(l)^2 -> eta \/ 2$). Per topology, in the form @eq:N0N1,
$ N_0^"(a)" = N_0^"(b)" = 0, quad
  N_1^"(a)" = -i (C_f - 1/2 C_"ad"^((a))) Tg^a_f slb(q) Pc, quad
  N_1^"(b)" = -i/2 thin cal(M)^a_f slb(q) Pc , $ <eq:g1_cFF-N-ab>
$ N_0^"(c)" = N_0^"(c1)" + N_0^"(c2)" = -i/2 C_"ad"^((a)) Tg^a_f (slh(p) Pcb - slh(p)' Pc), \
  N_1^"(c1)" = i/12 C_"ad"^((a)) Tg^a_f (2 slb(p) - slb(p)') Pc, quad
  N_1^"(c2)" = i/12 C_"ad"^((a)) Tg^a_f (slb(p) - 2 slb(p)') Pc, quad
  N_1^"(c)" = -i/4 C_"ad"^((a)) Tg^a_f slb(q) Pc . $ <eq:g1_cFF-N-c>
$N_0^"(c)"$ is $-i/2 thin C_"ad"^((a))$ times the Feynman rule $i V_(Dbrk_1)$ of @eq:V-delta, i.e. $1/2 thin C_"ad"^((a)) V_(Dbrk_1)$; its two terms have an external sterile leg and are evanescent. All $N_1$ are four-dimensional, $Pphys N_1 = N_1$.

*Result.* Summing the three topologies, with $C_"ad"^((B)) = 0$:
#result[
  $ (Dbrk_1 dot Gamma)^((1))_(gh^a overline(Psi)_f Psi_f) |_"pole" = lf 1/ep N_0, quad quad
    N_0 = -i/2 C_"ad"^((a)) Tg^a_f (slh(p) P_(overline(c)_f) - slh(p)' P_(c_f)) quad "(evanescent)" , $ <eq:g1_cFF-pole>
  the vertex function of the local functional
  $ (Dbrk_1 dot Gamma)^((1))|_("pole", gh overline(Psi) Psi) = lf 1/ep sum_f sum_(r = G, W) 1/2 C_"ad"^((r)) thin Dbrk_1^((f r)), \
    Dbrk_1^((f r)) = integral dif^d x thin gh^a_((r)) [ overline(Psi)_f Tg^a_f ghat^mu P_(overline(c)_f) partial_mu Psi_f + (partial_mu overline(Psi)_f) ghat^mu P_(c_f) Tg^a_f Psi_f ] , $ <eq:g1_cFF-pole-op>
  where $Dbrk_1^((f r))$ is the part of @eq:Delta1 with multiplet $f$ and a ghost of factor $r$ (vertex function $-i Tg^a_f (slh(p) Pcb - slh(p)' Pc)$).
  $ X_(1; gh^a overline(Psi)_f Psi_f) = lf thin Pphys N_1, quad quad
    Pphys N_1 = -i [ (C_f - 1/4 C_"ad"^((a))) Tg^a_f + 1/2 cal(M)^a_f ] slb(q) P_(c_f) , $ <eq:g1_cFF-X>
  $ X_1 & supset - lf sum_f integral dif^4 x thin overline(psi)_f gamma^mu (partial_mu [ C_f wm - 1/4 sum_(r) C_"ad"^((r)) wm^((r)) + 1/2 sum_(f', j) bold(Y)^j_(f f') wm_(f') (bold(Y)^j_(f f'))^dagger ]) psi_f \
    & = lf sum_f sum_(r in f) x_(f r) thin N^(f r)_F , $ <eq:g1_cFF-X-op>
  $ x_(f r) = x^"(a)"_(f r) + x^"(b)"_(f r) + x^"(c)"_(f r) = -(C_f - 1/2 C_"ad"^((r))) - 1/2 nu_(f r) - 1/4 C_"ad"^((r)) , $ <eq:g1_cFF-x>
  with $nu_(f r)$ defined by $cal(M)^a_f = nu_(f r) Tg^a_f$ for $a$ in the factor $r$ (@tab:g1_cFF-sm).
]
Translation: $N^(f r)_F = integral overline(psi)_f gamma^mu (partial_mu wm^((r))) psi_f = - integral (partial_mu gh^a) overline(psi)_f gamma^mu Tg^a_f psi_f$; with the ghost first and $partial -> -i q$ on it, its vertex function is $+i slb(q) thin Tg^a_f P_(c_f)$, in agreement with $bop Q^(f r)_V = -N^(f r)_F$ and example 2 of @sec:dictionary. In $wm_(f') = Tg^a_(f') gh^a$ the generators of the partner act on the indices it shares with $f$. The pole @eq:g1_cFF-pole-op is not in the basis of @sec:evanescent (which lists ghost-number-0 operators); it is the evanescent part of $bop_d Q^(f r)_R$: by @eq:m-bYR, $bop_d Q^(f r)_R |_(gh overline(Psi) Psi) = Dbrk_1^((f r)) + N^(f r)_F$ with $Q^(f r)_R$ continued to $d$ dimensions as in @eq:Lext-explicit.

#figure(
  table(
    columns: 7,
    align: (center, center, center, left, left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$f$], [$r$], [$nu_(f r)$], [$x^"(a)"_(f r)$], [$x^"(b)"_(f r)$], [$x^"(c)"_(f r)$], [$x_(f r)$ (coefficient of $N^(f r)_F$ in $16 pi^2 X_1$)]),
    [$Q$], [$G$], [$y_u^2 + y_d^2$], [$1/6 gs^2 - 3/4 gw^2 - 1/36 gy^2$], [$-(y_u^2 + y_d^2)\/2$], [$-3/4 gs^2$], [$-7/12 gs^2 - 3/4 gw^2 - 1/36 gy^2 - 1/2 (y_u^2 + y_d^2)$],
    [$Q$], [$W$], [$0$], [$-4/3 gs^2 + 1/4 gw^2 - 1/36 gy^2$], [$0$], [$-1/2 gw^2$], [$-4/3 gs^2 - 1/4 gw^2 - 1/36 gy^2$],
    [$Q$], [$B$], [$4 y_u^2 - 2 y_d^2$], [$-4/3 gs^2 - 3/4 gw^2 - 1/36 gy^2$], [$-2 y_u^2 + y_d^2$], [$0$], [$-4/3 gs^2 - 3/4 gw^2 - 1/36 gy^2 - 2 y_u^2 + y_d^2$],
    [$U$], [$G$], [$2 y_u^2$], [$1/6 gs^2 - 4/9 gy^2$], [$-y_u^2$], [$-3/4 gs^2$], [$-7/12 gs^2 - 4/9 gy^2 - y_u^2$],
    [$U$], [$B$], [$y_u^2 \/ 2$], [$-4/3 gs^2 - 4/9 gy^2$], [$-y_u^2 \/ 4$], [$0$], [$-4/3 gs^2 - 4/9 gy^2 - 1/4 y_u^2$],
    [$D$], [$G$], [$2 y_d^2$], [$1/6 gs^2 - 1/9 gy^2$], [$-y_d^2$], [$-3/4 gs^2$], [$-7/12 gs^2 - 1/9 gy^2 - y_d^2$],
    [$D$], [$B$], [$-y_d^2$], [$-4/3 gs^2 - 1/9 gy^2$], [$y_d^2 \/ 2$], [$0$], [$-4/3 gs^2 - 1/9 gy^2 + 1/2 y_d^2$],
    [$L$], [$W$], [$0$], [$1/4 gw^2 - 1/4 gy^2$], [$0$], [$-1/2 gw^2$], [$-1/4 gw^2 - 1/4 gy^2$],
    [$L$], [$B$], [$2 y_e^2$], [$-3/4 gw^2 - 1/4 gy^2$], [$-y_e^2$], [$0$], [$-3/4 gw^2 - 1/4 gy^2 - y_e^2$],
    [$E$], [$B$], [$y_e^2$], [$-gy^2$], [$-y_e^2 \/ 2$], [$0$], [$-gy^2 - 1/2 y_e^2$],
  ),
  caption: [Standard-Model values of @eq:g1_cFF-x, with $C_f$ and $C_"ad"^((r))$ of @tab:invariants. $nu_(f r)$: for $f = Q$, $cal(M)^a_Q = y_u^2 Tg^a_U + y_d^2 Tg^a_D$ (times $bb(1)$ in $"SU"(2)_L$); for $f = U, D$, $cal(M)^a_f = y_(u,d)^2 Tr_2 Tg^a_Q$ (trace over the $"SU"(2)_L$ index); for $L$, $cal(M)^a_L = y_e^2 Tg^a_E$; for $E$, $cal(M)^a_E = y_e^2 Tr_2 Tg^a_L$. The pole @eq:g1_cFF-pole-op has the coefficient $C_"ad"^((r)) \/ 2 = 3/2 gs^2$ ($r = G$, quarks), $gw^2$ ($r = W$, doublets), $0$ ($r = B$) of $Dbrk_1^((f r))$, in units of $1 \/ (16 pi^2 ep)$.],
) <tab:g1_cFF-sm>

#check[*Cancellation of the pole.* By @sec:m-oneloop, $sum_G Kop Top(1) G + i V_(bop_d delta S^((1))) = 0$, where $V$ is the vertex function of @sec:dictionary and $i V$ the Feynman rule. The loops give $ilf N_0 \/ ep = lf (1\/ep) (1/2 C_"ad"^((a))) thin i V_(Dbrk_1)$, with $i V_(Dbrk_1) = Tg^a_f (slh(p) Pcb - slh(p)' Pc)$ the Feynman rule @eq:V-delta. The only evanescent $gh overline(Psi) Psi$ part of $V_(bop_d delta S^((1)))$ is $z^((r))_R V_(Dbrk_1)$, from the source counterterm $z^((r))_R Q^(f r)_R$ (@eq:m-bYR). The one-loop graph of $Gamma_(overline(srcF) gh Psi)$ is (c2) with the $Dbrk_1$ vertex replaced by $-Tg^b Pc$ of @eq:V-src-matter. Its value is
$ -i/2 C_"ad"^((r)) Tg^a_f intk (Pc sl(k) (slb(p)' - slb(k)) Pc) / (k^2 (p' - k)^2 (k - p)^2) , $
the numerator gives $-bar(k)^2 Pc$ at degree zero, and the pole is $Gamma^((1))|_"pole" = lf (1\/ep) (1/2 C_"ad"^((r))) Gamma^((0))_(overline(srcF) gh Psi)$ with $Gamma^((0))_(overline(srcF) gh Psi) = i Tg^a_f Pc$, in agreement with @sec:g0_ghost_sources. Hence $z^((r))_R = -lf (1\/ep) C_"ad"^((r)) \/ 2$, and $i z^((r))_R V_(Dbrk_1)$ cancels the loops exactly. The four-dimensional pole is zero in @eq:g1_cFF-pole. Correspondingly the coefficients of $N^(f r)_F$ in $16 pi^2 ep thin bop_d delta S^((1))$ cancel: $-(C_f + cal(Y)_f)$ from $Q^f_K$ (@sec:g0_FF), $+(C_f + C_"ad"^((r)) + cal(Y)_f)$ from $Q^(f r)_V$ (@sec:g0_FFV, $bop Q^(f r)_V = -N^(f r)_F$), $-1/2 C_"ad"^((r))$ from $Q^(f r)_R$, and $-1/2 C_"ad"^((r))$ from $Q^((r))_(rho partial)$, which enters $16 pi^2 ep thin delta S^((1))$ with $+1/2 C_"ad"^((r))$ (@sec:g0_ghost_sources) and has $bop Q^((r))_(rho partial) supset - sum_f N^(f r)_F$ (@sec:btable). Here $cal(Y)_f$ is the Yukawa coefficient of the fermion self-energy pole (@eq:sts-yuk): $cal(Y)_Q = (y_u^2 + y_d^2)\/2$, $cal(Y)_U = y_u^2$, $cal(Y)_D = y_d^2$, $cal(Y)_L = y_e^2\/2$, $cal(Y)_E = y_e^2$. With $delta S^((1))$ recomputed from its own graphs (fermion self-energy, $-i q_mu$ times the fermion--gauge vertex, $Gamma_(gh srcA)$, $Gamma_(overline(srcF) gh Psi)$, $Gamma_(gh overline(Psi) srcF)$), both the evanescent and the four-dimensional part of the identity hold for each of the ten pairs $(f, r)$, Yukawa terms included.]

#check[*Rigid invariance and locality.* At $q = 0$ the ghost is constant and $Pphys X$ must vanish (@eq:dict-ii). (a) and (b) are proportional to $slb(q)$ separately; (c1) and (c2) are not, $N_1^"(c1)" + N_1^"(c2)" |_(q = 0) = (i/12 - i/12) C_"ad"^((a)) Tg^a_f slb(p) Pc = 0$. The sum is local, independent of $MUV$ and of $m^2$ (both of degree two, beyond $dg(G) = 1$). $g5$ enters only through $Pc$: left- and right-handed multiplets have the same form, and no $epsilon$ tensor appears.]

#check[*Known limits.* For a simple gauge group without scalars, @eq:g1_cFF-x is $x = -(C_f - 1/4 C_"ad")$, with $-(C_f - 1/2 C_"ad")$ from the gauge exchange and $-1/4 C_"ad"$ from the two internal-ghost graphs, which also carry the whole pole $1/2 C_"ad" Dbrk_1$. This is the Feynman-gauge value of @Belusca-Maito:2020ala and of @Lazopoulos:2026two (there $(5 + xi)\/6 thin C_f - xi thin C_"ad"\/4$). The Yukawa term $-1/2 nu_(f r)$ agrees with the scalar-exchange term of @Belusca-Maito:2020ala, rewritten for Dirac multiplets: $cal(M)^a_f$ is the sum over the Yukawa partners of $Y Tg^a_(f') Y^dagger$.]

#check[*Independent evaluations.* The Dirac identities behind @eq:g1_cFF-G-a to @eq:g1_cFF-G-c2 were evaluated with an explicit $8 times 8$ representation of the six-dimensional Clifford algebra (four barred, two hatted directions) for both chiralities, and the raw products of Feynman rules were compared with the symbolic BMHV algebra at $d = 6, 8$. $N_0$ and $Pphys N_1$ were recomputed with the expansion operator of @sec:integrals for all ten pairs $(f, r)$; they agree per graph with @eq:g1_cFF-N-ab and @eq:g1_cFF-N-c, do not change under shifts of the loop momentum, and do not depend on $MUV$. The group factors of @tab:g1_cFF-diagrams and the values $nu_(f r)$ were evaluated numerically with explicit $"SU"(3)$, $"SU"(2)$ matrices and $epsilon_(i j)$ (@sec:supporting). In the solution of @sec:st-solution, $N^(f r)_F$ is produced by $bop Q^f_K$, $bop Q^(f r)_V$, $bop Q^(f r)_R$ and $bop Q^((r))_(rho partial)$ (@eq:btable-fermion, @eq:btable-ghost); $Q^(f r)_R$ and $Q^((r))_(rho partial)$ also produce structures of @sec:g1_cFFV and @sec:g1_ccRF, where the same coefficients must reappear.]

#general[The pole proportional to $Dbrk_1$ itself is the one-loop case of a general statement: a source counterterm $z^((L))_R Q_R$ of $delta S^((L))$ has $bop_d$ variation $z^((L))_R (Dbrk_1 + N_F)$ (@eq:m-bYR), so the pole of $(Dbrk^([l]) dot Gamma)_(gh overline(Psi) Psi)$ contains $Dbrk_1$ at every order. At two loops the vertex $z^((1))_R Dbrk_1$ of $bop_d delta S^((1))$ is inserted into one-loop graphs (class (c) of @sec:m-beyond); next to an internal sterile line it produces $eta \/ ep^2$ terms like the tree $Dbrk_1$. A finite source counterterm in $Delta S^((1))$ enters in the same way, through $Evan bop_d Delta S^((1))$ of @eq:m-LP. In this class $Evan [ilf N_1] = 0$ at one loop.]
