#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Fermion--gauge vertices $overline(Psi) Psi A$ <sec:g0_FFV>

External legs $Psi_f (p)$ and $A^a_mu (q)$ incoming and $overline(Psi)_f (p')$ outgoing with $p' = p + q$, at $p = p' = q = 0$; $mu$ is a #box[$d$-dimensional] index. The sum of the diagrams below is $i thin Gamma^((1))_(overline(Psi)_f Psi_f A^a_mu)$, a matrix in Dirac, colour and $"SU"(2)_L$ indices with $overline(Psi)_f$ on the left. The superficial degree is $dg(G) = 4 + 2 - 2 dot 3 = 0$ for every diagram, so the pole is momentum independent and zero external momenta suffice. The class covers the ten pairs (multiplet, gauge factor) with $Tg^a_f != 0$: $Q G$, $Q W$, $Q B$, $U G$, $U B$, $D G$, $D B$, $L W$, $L B$, $E B$. The tree-level vertex function is $Gamma^((0))_(overline(Psi)_f Psi_f A^a_mu) = gbar^mu Pc Tg^a_f$ (@eq:V-ffv divided by $i$), with $Pc = P_(c_f)$, $Pcb = P_(overline(c)_f)$.

*Diagrams.* Four topologies (@fig:g0_FFV-diagrams). In (c) and (d) the internal fermion is a Yukawa partner $f'$ of $f$ and the internal scalar is $cal(H) in {H, tilde(H)}$, the doublet for which the Yukawa term of @eq:LYuk reads $overline(f) thin cal(H) thin f'$ or its conjugate. In this subsection $cal(H)$ denotes the doublet whose charge flows _along_ the fermion line, from the $Psi_f$ vertex to the $overline(Psi)_f$ vertex (a doublet field, $H$ or $tilde(H)$, not the field operator $phi_j in {H_j, H^dagger_j}$ that labels the Yukawa vertices in @sec:g0_FF and @sec:g0_FFS); it carries the coupling-weighted generators $Tg^a_H$ of @eq:Tg-explicit or
$ (Tg^a_(tilde(H)))_(i j) = - epsilon_(i k) (Tg^a_H)_(l k) epsilon_(j l), quad "i.e." quad Tg^I_(tilde(H)) = gw tau^I/2, quad Tg^B_(tilde(H)) = -gy/2 . $ <eq:g0_FFV-THt>
$Tg_(tilde(H))$ is the generator of $tilde(H)_i = epsilon_(i j) H^*_j$, and $Tg^a_f = Tg^a_(f') + Tg^a_(cal(H))$ on the product representation (gauge invariance of @eq:LYuk).

#figure(
  diagrow(
    oneloop(
      ((legs: (("V", $A^a_mu$, 0),)), (legs: (("F", $overline(Psi)_f$, -1),)), (legs: (("F", $Psi_f$, 1),))),
      (("F", $Psi_f$, 1), ("V", $A^b$, 0), ("F", $Psi_f$, 1)), caption: [(a)]),
    oneloop(
      ((legs: (("V", $A^a_mu$, 0),)), (legs: (("F", $overline(Psi)_f$, -1),)), (legs: (("F", $Psi_f$, 1),))),
      (("V", $A^b$, 0), ("F", $Psi_f$, -1), ("V", $A^c$, 0)), caption: [(b)]),
    oneloop(
      ((legs: (("V", $A^a_mu$, 0),)), (legs: (("F", $overline(Psi)_f$, -1),)), (legs: (("F", $Psi_f$, 1),))),
      (("F", $Psi_(f')$, 1), ("S", $cal(H)$, 0), ("F", $Psi_(f')$, 1)), caption: [(c)]),
    oneloop(
      ((legs: (("V", $A^a_mu$, 0),)), (legs: (("F", $overline(Psi)_f$, -1),)), (legs: (("F", $Psi_f$, 1),))),
      (("S", $cal(H)$, 0), ("F", $Psi_(f')$, -1), ("S", $cal(H)$, 0)), caption: [(d)]),
  ),
  caption: [Topologies of the class $overline(Psi)_f Psi_f A^a_mu$. (a) gauge exchange; (b) three-gauge vertex; (c) scalar exchange, gauge boson on the internal fermion; (d) scalar exchange, gauge boson on the internal scalar.],
) <fig:g0_FFV-diagrams>

#figure(
  table(
    columns: 5,
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [symmetry factor], [sign], [group factor]),
    [(a)], [$Psi_f$, $Psi_f$, $A^b$ (all twelve $b$)], [$1$], [$+1$], [$Tg^b_f Tg^a_f Tg^b_f = (C_f - 1/2 C_"ad"^((a))) Tg^a_f$],
    [(b)], [$Psi_f$, $A^b$, $A^c$ ($b, c$ in the factor of $a$)], [$1$], [$+1$], [$fg^(a b c) Tg^b_f Tg^c_f = i/2 C_"ad"^((a)) Tg^a_f$],
    [(c)], [$Psi_(f')$, $Psi_(f')$, $cal(H)$], [$1$], [$+1$], [$y^2 thin [Tg^a_(f')]_f$, @tab:g0_FFV-yuk],
    [(d)], [$Psi_(f')$, $cal(H)$, $cal(H)$], [$1$], [$+1$], [$y^2 thin [Tg^a_(cal(H))]_f$, @tab:g0_FFV-yuk],
  ),
  caption: [Field assignments, symmetry factors, signs and group factors. No closed loop, hence no fermion-loop sign. In (b) the ordered pair $(b, c)$ is summed over; the two gauge lines end on different vertices, so there is no symmetry factor.],
) <tab:g0_FFV-diagrams>

The Yukawa diagrams exist for the following $(f; f', cal(H), y)$. $[X]_f$ denotes the matrix $X$ of the partner's indices carried onto the indices of $f$ by the two Yukawa vertices and the scalar propagator:

#figure(
  table(
    columns: 5,
    align: (center, center, left, left, center),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$f$], [$(f', cal(H), y)$], [$[Tg^a_(f')]_f$ in (c)], [$[Tg^a_(cal(H))]_f$ in (d)], [$n$]),
    [$Q$], [$(U, tilde(H), y_u)$, $(D, H, y_d)$], [$Tg^a_U bb(1)_2$, $Tg^a_D bb(1)_2$], [$Tg^a_(tilde(H))$, $Tg^a_H$], [$1$],
    [$L$], [$(E, H, y_e)$], [$Tg^a_E bb(1)_2$], [$Tg^a_H$], [$1$],
    [$U$], [$(Q, H, y_u)$], [$Tr_2 Tg^a_Q$], [$Tr_2 Tg^a_H$], [$2$],
    [$D$], [$(Q, tilde(H), y_d)$], [$Tr_2 Tg^a_Q$], [$Tr_2 Tg^a_(tilde(H))$], [$2$],
    [$E$], [$(L, tilde(H), y_e)$], [$Tr_2 Tg^a_L$], [$Tr_2 Tg^a_(tilde(H))$], [$2$],
  ),
  caption: [Yukawa partners. $Tr_2$ is the trace over the $"SU"(2)_L$ index of the internal doublet; $Tg^a_U bb(1)_2$ acts on the colour index of $Q$ and trivially on its $"SU"(2)_L$ index. Colour indices are $delta^alpha_beta$ at every Yukawa vertex. $n$ is the number of $"SU"(2)_L$ components running in the loop.],
) <tab:g0_FFV-yuk>

The entries follow from @eq:V-yuk: for $f = Q$, $f' = U$ the $"SU"(2)_L$ factor is $epsilon_(i j) epsilon_(i' j) = delta_(i i')$ in (c) and $epsilon_(i j') (Tg^a_H)_(j j') epsilon_(i' j) = -(Tg^a_(tilde(H)))_(i i')$ in (d), the minus sign being compensated by the orientation of the scalar momentum (below); for $f = U$, $f' = Q$ it is $epsilon_(i j) (Tg^a_Q)_(i' i) epsilon_(i' j) = Tr_2 Tg^a_Q$ in (c) and $epsilon_(i j) (Tg^a_H)_(j' j) epsilon_(i j') = Tr_2 Tg^a_H$ in (d). For $f = Q, L$ with $f' = D, E$, and for $f = D, E$, the $epsilon$'s are replaced by $delta$'s. In every case $[Tg^a_(f')]_f + [Tg^a_(cal(H))]_f = n thin Tg^a_f$.

*Integrands.* At zero external momenta every internal fermion carries the loop momentum $k$ along the fermion flow. In (a) the gauge line carries $k$ from the $overline(Psi)_f$ vertex to the $Psi_f$ vertex. In (b) the line $A^b_beta$ joins the $overline(Psi)_f$ vertex (index $rho$) to the three-gauge vertex and $A^c_gamma$ the $Psi_f$ vertex (index $nu$) to it; the incoming momenta at the three-gauge vertex are $(0, k, -k)$ for $(A^a_mu, A^b_beta, A^c_gamma)$. In (c) and (d) the scalar carries $-k$ along the fermion flow. Reading against the fermion flow,
$ G_"(a)" = intk [i gbar^rho Pc Tg^b] (i sl(k))/(k^2) [i gbar^mu Pc Tg^a] (i sl(k))/(k^2) [i gbar^nu Pc Tg^b] (-i g_(nu rho))/(k^2) , $ <eq:g0_FFV-int-a>
$ G_"(b)" = intk [i gbar_beta Pc Tg^b] (i sl(k))/(k^2) [i gbar_gamma Pc Tg^c] (-i)/(k^2) (-i)/(k^2) thin fg^(a b c) [ - g^(mu beta) k^gamma + 2 g^(beta gamma) k^mu - g^(gamma mu) k^beta ] , $ <eq:g0_FFV-int-b>
$ G_"(c)" = intk [-i y Pcb] (i sl(k))/(k^2) [i gbar^mu Pcb thin Tg^a_(f')] (i sl(k))/(k^2) [-i y Pc] (i)/(k^2 - m^2) , $ <eq:g0_FFV-int-c>
$ G_"(d)" = intk [-i y Pcb] (i sl(k))/(k^2) [-i y Pc] (i)/(k^2 - m^2) [i (-2 k^mu) Tg^a_(cal(H))] (i)/(k^2 - m^2) , $ <eq:g0_FFV-int-d>
where in (b) $fg^(a b c)[dots]$ is @eq:V-vvv at $(p_1, p_2, p_3) = (0, k, -k)$, and in (c), (d) the group matrices are understood as $[Tg^a_(f')]_f$, $[Tg^a_(cal(H))]_f$ of @tab:g0_FFV-yuk. The factor $-2 k^mu$ in (d) is $(p + p')^mu$ of @eq:V-ssv: for $cal(H) = H$ the $H$ line runs along the fermion flow with $p = p' = -k$; for $cal(H) = tilde(H)$ the $H$ line runs against it with $p = p' = k$, and the extra sign of the $"SU"(2)_L$ contraction turns $Tg^a_H$ into $Tg^a_(tilde(H))$. The overall factors are $i^5 (-i) = 1$ in (a), $i^3 (-i)^2 = i$ in (b) and $(-i y)^2 i^4 = -y^2$ in (c) and (d).

*$Rop$-operation.* At one loop the forest is empty, $Rbar G = G$, and $Z(G) = -Kop Top(0) G$ (@eq:m-Z1). With $dg(G) = 0$, $Top(0)$ keeps only the leading term of @eq:uv-expansion: every denominator becomes $k^2 - MUV^2$, and the external momenta and $m^2$ (degree 1 and 2) drop out; this is why the integrands are written at zero external momenta. Every numerator below is proportional to $k^2$, and
$ Kop intk (k^2)/((k^2 - MUV^2)^3) = ilf 1/ep . $ <eq:g0_FFV-master>
_Dirac algebra._ By @eq:phys-prop each projector sandwich gives a four-dimensional numerator, $Pc sl(k) Pcb = slb(k) Pcb$, so that, with all projectors moved to the right,
$ "(a)": quad gbar^rho slb(k) gbar^mu slb(k) gbar_rho Pc = -2 slb(k) gbar^mu slb(k) Pc, quad quad
  "(c)": quad Pcb sl(k) gbar^mu Pcb sl(k) Pc = slb(k) gbar^mu slb(k) Pc, $
$ "(b)": quad gbar_beta slb(k) gbar_gamma [dots]^(mu beta gamma) Pc = (-2 bar(k)^2 gbar^mu - 4 k^mu slb(k)) Pc, quad quad
  "(d)": quad -2 k^mu Pcb sl(k) Pc = -2 k^mu slb(k) Pc . $
All contractions are between barred matrices (@eq:contr2); no $ghat$ and no $gmh^nu_(space nu)$ occurs. _Tensor reduction_ (@eq:tensor-red, @eq:tensor-red-split): $slb(k) gbar^mu slb(k) -> (k^2 \/ d) gbar_alpha gbar^mu gbar^alpha = -(2 \/ d) k^2 gbar^mu$, $bar(k)^2 -> (4 \/ d) k^2$, $k^mu slb(k) -> (k^2 \/ d) gbar^mu$. The factors $4$ are traces of $gmb$, not of $gmh$, so no $eta$ is produced. With the overall factors and group factors of @tab:g0_FFV-diagrams:
$ Top(0) G_"(a)" = 4/d (C_f - 1/2 C_"ad"^((a))) gbar^mu Pc Tg^a_f thin I, quad quad
  Top(0) G_"(b)" = i dot i/2 C_"ad"^((a)) dot (- 12/d) gbar^mu Pc Tg^a_f thin I = 6/d C_"ad"^((a)) gbar^mu Pc Tg^a_f thin I, \
  Top(0) G_"(c)" = (-y^2)(-2/d) gbar^mu Pc [Tg^a_(f')]_f thin I, quad quad
  Top(0) G_"(d)" = (-y^2)(-2/d) gbar^mu Pc [Tg^a_(cal(H))]_f thin I, $
with $I = intk k^2 \/ (k^2 - MUV^2)^3$. With @eq:g0_FFV-master and $d -> 4$ in the residue, every topology $X$ takes the form @eq:N0N1, $Top(0) G_X = ilf (1/ep N_0^X + eta/ep N_1^X) + dots$, with
$ N_0^"(a)" = (C_f - 1/2 C_"ad"^((a))) gbar^mu Pc Tg^a_f, quad
  N_0^"(b)" = 3/2 C_"ad"^((a)) gbar^mu Pc Tg^a_f, \
  N_0^"(c)" = y^2/2 gbar^mu Pc [Tg^a_(f')]_f, quad
  N_0^"(d)" = y^2/2 gbar^mu Pc [Tg^a_(cal(H))]_f, $ <eq:g0_FFV-N>
and $N_1^X = 0$ for every topology. Summed over the Yukawa partners, $N_0^"(c)" + N_0^"(d)" = cal(Y)_f thin gbar^mu Pc Tg^a_f$ by $[Tg^a_(f')]_f + [Tg^a_(cal(H))]_f = n thin Tg^a_f$, with
$ cal(Y)_f = 1/2 sum_((f', cal(H), y)) n thin y^2 : quad
  cal(Y)_Q = (y_u^2 + y_d^2)/2, quad cal(Y)_U = y_u^2, quad cal(Y)_D = y_d^2, quad cal(Y)_L = y_e^2 / 2, quad cal(Y)_E = y_e^2 . $ <eq:g0_FFV-Y>
These are the Yukawa coefficients $cal(Y)_f$ of the fermion self-energy pole in @sec:g0_FF, the Yukawa sums of @sec:st-solution.

*Result.*
#result[
  $ Gamma^((1))_(overline(Psi)_f Psi_f A^a_mu) |_"pole" = lf 1/ep N_0, quad quad
    N_0 = z_(f r) thin gbar^mu Pc Tg^a_f, quad quad z_(f r) = C_f + C_"ad"^((r)) + cal(Y)_f , $ <eq:g0_FFV-pole>
  for $a$ in the gauge factor $r$; $z_(f r)$ is a symbol of this subsection. In terms of the local operators $Q^(f r)_V$ of @tab:Qbasis, continued to $d$ dimensions as $integral dif^d x thin overline(Psi)_f gbar^mu Am^((r))_mu Pc Psi_f$ (vertex function $gbar^mu Pc Tg^a_f$),
  $ delta S^((1)) supset - lf 1/ep sum_f sum_(r in f) z_(f r) integral dif^d x thin overline(Psi)_f thin gbar^mu Am^((r))_mu thin Pc thin Psi_f = - lf 1/ep sum_f sum_(r in f) z_(f r) thin Q^(f r)_V . $ <eq:g0_FFV-ct>
  The pole has no evanescent part: no $ghat^mu$, no $hat(p)$. $N_1 = 0$. In terms of the constants of @sec:st-solution, $-lf z_(f r) \/ ep = delta Z_f + delta Z_(c_r) + 1/2 delta Z^((r))_A$ is the coefficient of $Q^(f r)_V$ in the assembled counterterm (see the check below).
]

#figure(
  table(
    columns: 6,
    align: (center, left, left, left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$f r$], [(a)], [(b)], [(c)], [(d)], [$z_(f r)$]),
    [$Q G$], [$-1/6 gs^2 + 3/4 gw^2 + 1/36 gy^2$], [$9/2 gs^2$], [$(y_u^2 + y_d^2)/2$], [$0$], [$13/3 gs^2 + 3/4 gw^2 + 1/36 gy^2 + (y_u^2 + y_d^2)/2$],
    [$Q W$], [$4/3 gs^2 - 1/4 gw^2 + 1/36 gy^2$], [$3 gw^2$], [$0$], [$(y_u^2 + y_d^2)/2$], [$4/3 gs^2 + 11/4 gw^2 + 1/36 gy^2 + (y_u^2 + y_d^2)/2$],
    [$Q B$], [$4/3 gs^2 + 3/4 gw^2 + 1/36 gy^2$], [$0$], [$2 y_u^2 - y_d^2$], [$3/2 (y_d^2 - y_u^2)$], [$4/3 gs^2 + 3/4 gw^2 + 1/36 gy^2 + (y_u^2 + y_d^2)/2$],
    [$U G$], [$-1/6 gs^2 + 4/9 gy^2$], [$9/2 gs^2$], [$y_u^2$], [$0$], [$13/3 gs^2 + 4/9 gy^2 + y_u^2$],
    [$U B$], [$4/3 gs^2 + 4/9 gy^2$], [$0$], [$y_u^2 \/ 4$], [$3 y_u^2 \/ 4$], [$4/3 gs^2 + 4/9 gy^2 + y_u^2$],
    [$D G$], [$-1/6 gs^2 + 1/9 gy^2$], [$9/2 gs^2$], [$y_d^2$], [$0$], [$13/3 gs^2 + 1/9 gy^2 + y_d^2$],
    [$D B$], [$4/3 gs^2 + 1/9 gy^2$], [$0$], [$-y_d^2 \/ 2$], [$3 y_d^2 \/ 2$], [$4/3 gs^2 + 1/9 gy^2 + y_d^2$],
    [$L W$], [$-1/4 gw^2 + 1/4 gy^2$], [$3 gw^2$], [$0$], [$y_e^2 \/ 2$], [$11/4 gw^2 + 1/4 gy^2 + y_e^2 \/ 2$],
    [$L B$], [$3/4 gw^2 + 1/4 gy^2$], [$0$], [$y_e^2$], [$-y_e^2 \/ 2$], [$3/4 gw^2 + 1/4 gy^2 + y_e^2 \/ 2$],
    [$E B$], [$gy^2$], [$0$], [$y_e^2 \/ 2$], [$y_e^2 \/ 2$], [$gy^2 + y_e^2$],
  ),
  caption: [Standard-Model values. Each entry is the coefficient of $gbar^mu Pc Tg^a_f$ in $N_0$ of the topology ($a in r$), with $C_f$ and $C_"ad"^((r))$ from @tab:invariants. In the columns (c), (d) for $r = B$ the entries are $n y^2 Y_(f') \/ (2 Y_f)$ and $n y^2 Y_(cal(H)) \/ (2 Y_f)$ summed over the partners, with $Y_(tilde(H)) = -1/2$; for $r = G, W$ the partner that does not carry $r$ gives zero.],
) <tab:g0_FFV-sm>

#check[*Abelian Ward identity.* For $r = B$, $C_"ad"^((B)) = 0$ and $z_(f B) = C_f + cal(Y)_f$: the $B$ vertex equals $Tg^B_f thin partial \/ partial p_mu$ of the fermion self-energy pole of @sec:g0_FF, which must be $lf (1\/ep) (C_f + cal(Y)_f) slb(p) Pc$ (gauge loop $C_f$, scalar loop $cal(Y)_f$). This is $Z_1 = Z_2$ for $"U"(1)_Y$. It holds separately for the Yukawa part, topology by topology in the sum (c) + (d), by $Y_(f') + Y_(cal(H)) = Y_f$ (column sums of @tab:g0_FFV-sm, e.g. $Q B$: $2 y_u^2 - y_d^2 + 3/2 (y_d^2 - y_u^2) = (y_u^2 + y_d^2)/2$).]

#check[*Relation to the renormalisation constants.* With the bare fields and couplings of @sec:st-solution ($psi_(f 0) = Z_f^(1\/2) psi_f$, $A_0 = (Z^((r))_A)^(1\/2) A$, $c_(r 0) = mu^ep Z_(c_r) c_r$ with $c_r$ the coupling of the factor $r$, @sec:bases), the invariant fermion part of $delta S^((1))$ is $sum_f [delta Z_f Q^f_K + sum_(r in f) (delta Z_f + delta Z_(c_r) + 1/2 delta Z^((r))_A) Q^(f r)_V]$, i.e. $sum_f 1/2 delta Z_f cal(O)^f_psi + sum_r (delta Z_(c_r) cal(O)^((r))_g + 1/2 delta Z^((r))_A cal(O)^((r))_A)$ restricted to $Q_K$, $Q_V$ (@eq:closed; for $r = B$, $delta Z_(c_B) + 1/2 delta Z^((B))_A = 0$). Hence
$ delta Z_f + delta Z_(c_r) + 1/2 delta Z^((r))_A = - lf (z_(f r))/ep . $
With $delta Z_f = - lf (C_f + cal(Y)_f)\/ep$ from @sec:g0_FF this gives $delta Z_(c_r) + 1/2 delta Z^((r))_A = - lf C_"ad"^((r)) \/ ep$, independent of $f$ (universality of the gauge coupling), zero for $B$. The same combination follows from @sec:g0_VV, which fixes $delta Z^((r))_A$, and @sec:g0_VVV, which fixes $delta Z^((r))_(3 A) = delta Z_(c_r) + 3/2 delta Z^((r))_A$: $delta Z_(c_r) + 1/2 delta Z^((r))_A = delta Z^((r))_(3 A) - delta Z^((r))_A = - lf C_"ad"^((r)) \/ ep$ (@sec:g0_VV, @sec:g0_VVV). The fermion and scalar loops cancel in this combination, and the relation holds. It is the Slavnov--Taylor relation between the vertex constants, $Z_(1 F) \/ Z_f = Z^((r))_1 \/ Z^((r))_A = Z_(c_r) (Z^((r))_A)^(1\/2)$, with $Z_(1 F) = Z_f Z_(c_r) (Z^((r))_A)^(1\/2)$ the constant of $Q^(f r)_V$ and $Z^((r))_1 = 1 + delta Z^((r))_(3 A)$ that of $Q^((r))_(3 f)$.]

#check[*Known limits.* For one Dirac quark in QCD in Feynman gauge, $Z_(1 F) = 1 - lf (C_F + C_A) gs^2 \/ ep$; here the gluon part of $z_(Q G)$ is $C_Q^((G)) + C_"ad"^((G)) = (4/3 + 3) gs^2$, with $C_F = 4/3$, $C_A = 3$, and (a), (b) reproduce the standard split $C_F - C_A\/2$ and $3 C_A \/ 2$. For QED ($C_f = e^2$, $C_"ad" = 0$) the vertex pole is $+lf e^2 \/ ep$ times the tree vertex. For a simple gauge group with a scalar and Yukawa couplings in the same regularisation, @Belusca-Maito:2020ala (eqs. (5.10), (5.16)) give the vertex coefficient $C_f + C_"ad" + cal(Y)_f$, with no evanescent pole. A chiral multiplet gives the Dirac result because by @eq:phys-prop the numerator algebra is four-dimensional and $g5$ only enters through the projector $Pc$ at the end of the string.]

#check[*Auxiliary mass, Higgs mass, evanescent part.* $Top(0)$ contains only $intk k^2 \/ (k^2 - MUV^2)^3$, whose pole is independent of $MUV$; $m^2$ has degree 2 and does not enter. The pole is proportional to $gbar^mu$: its $ghat^mu$ component vanishes, consistent with the statement of @sec:evanescent that fermion bilinears (and their gauge vertices) carry no evanescent one-loop pole. $N_1 = 0$ as required for ghost number $0$ (@sec:m-K).]

#general[At one loop the pole is purely four-dimensional and its operator $Q^(f r)_V$ is continued to $d$ dimensions with $gbar^mu$, $bar(A)$ and $Pc$, as the tree vertex @eq:Lferm. From two loops on the same class contains the one-loop counterterm vertices $delta S^((1))$ (including the evanescent $E^((r))_A$, $E_H$ of @eq:evanescent) as subgraph insertions; $Rbar$ subtracts them, $hat(k)^2$ from an inserted $hat(partial)^2$ produces $eta\/ep^2$, and the residue then depends on the chosen #box[$d$-dimensional] continuation of the one-loop counterterms. The vertex pole can then also acquire non-invariant, non-evanescent parts.]
