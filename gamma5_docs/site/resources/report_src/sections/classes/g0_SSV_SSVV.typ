#import "../../macros.typ": *
#import "../../diagrams.typ": *

// local shorthands (this file only)
#let Yt = $tilde(Y)$
#let DM(n) = $(k^2 - MUV^2)^(#n)$
#let Ca(x) = $C_"ad"^((#x))$
#let tl(x) = $t^(#x)$
#let Ups = $cal(B)^(a b)$
#let sm(..ds) = oneloop(..ds, scale: 0.72, leg: 0.62)

=== Higgs--gauge functions $H^dagger H A$ and $H^dagger H A A$ <sec:g0_SSV_SSVV>

External legs $H_j (p)$ and $A^a_mu (q)$ incoming and $H^dagger_i$ outgoing with momentum $p' = p + q$; for four legs $H_j (p)$, $A^a_mu (q_1)$ and $A^b_nu (q_2)$ incoming and $H^dagger_i$ outgoing with $p' = p + q_1 + q_2$. The indices $mu$, $nu$ are #box[$d$-dimensional]. The sums of the diagrams are $i thin Gamma^((1))_(H^dagger_i H_j A^a_mu)$ (superficial degree $dg(G) = 1$) and $i thin Gamma^((1))_(H^dagger_i H_j A^a_mu A^b_nu)$ ($dg(G) = 0$), matrices in $(i, j)$. The tree-level vertex functions (@eq:V-ssv divided by $i$) are
$ Gamma^((0))_(H^dagger H A^a_mu) = (p + p')^mu Tg^a_H , quad quad Gamma^((0))_(H^dagger H A^a_mu A^b_nu) = g^(mu nu) {Tg^a_H, Tg^b_H} . $
The gauge factors are $W$ and $B$; the functions with an external gluon have no pole.

*Diagrams.* @fig:g0_SSV_SSVV-diagrams shows every topology once; @tab:g0_SSV_SSVV-diagrams lists the assignments. In the fermion loops $psi$ enters the $H_j$ vertex and $psi'$ leaves it, $(psi, psi') = (D, Q), (E, L), (Q, U)$ for $y_d, y_e, y_u$ (the same convention as in @sec:g0_SS). The gauge vertex sits on the $psi'$ line (from the $H_j$ vertex to the $H^dagger_i$ vertex along the fermion arrow) or on the $psi$ line.

#figure(
  stack(dir: ttb, spacing: 4mm,
  diagrow(
    sm(((legs: (("S", $H_j$, 1),)), (legs: (("V", $A^a$, 0),)), (legs: (("S", $H^dagger_i$, -1),))),
      (("F", [$psi'$#h(2mm)], 1), ("F", [#h(2.5mm)$psi'$], 1), ("F", [$psi$#h(2mm)], 1)), caption: [(F) fermion loop]),
    sm(((legs: (("S", $H_j$, 1),)), (legs: (("V", $A^a$, 0),)), (legs: (("S", $H^dagger_i$, -1),))),
      (("S", none, 1), ("S", none, 1), ("V", [$A^c$#h(2.5mm)], 0)), caption: [(a) gauge exchange]),
    sm(((legs: (("S", $H_j$, 1),)), (legs: (("V", $A^a$, 0),)), (legs: (("S", $H^dagger_i$, -1),))),
      (("V", [$A^c$#h(2.5mm)], 0), ("V", [#h(2.5mm)$A^d$], 0), ("S", none, -1)), caption: [(b) $A A A$ vertex]),
  ),
  diagrow(
    sm(((legs: (("S", $H_j$, 1), ("V", $A^a$, 0))), (legs: (("S", [#h(2mm)$H^dagger_i$], -1),))),
      (("S", none, 1), ("V", $A^c$, 0)), caption: [(c) seagull]),
    sm(((legs: (("S", $H_j$, 1), ("S", $H^dagger_i$, -1))), (legs: (("V", $A^a$, 0),))),
      (("V", none, 0), ("V", none, 0)), caption: [(d)]),
    sm(((legs: (("S", $H_j$, 1), ("S", $H^dagger_i$, -1))), (legs: (("V", $A^a$, 0),))),
      (("S", none, 1), ("S", none, 1)), caption: [(e)]),
  ),
  diagrow(
    sm(((legs: (("S", $H$, 1),)), (legs: (("V", $A^a$, 0),)), (legs: (("V", $A^b$, 0),)), (legs: (("S", $H^dagger$, -1),))),
      (("F", none, 1), ("F", none, 1), ("F", none, 1), ("F", none, 1)), caption: [(F4) fermion box]),
    sm(((legs: (("S", $H$, 1),)), (legs: (("V", $A^a$, 0),)), (legs: (("V", $A^b$, 0),)), (legs: (("S", $H^dagger$, -1),))),
      (("S", none, 1), ("S", none, 1), ("S", none, 1), ("V", none, 0)), caption: [(i)]),
    sm(((legs: (("S", $H$, 1),)), (legs: (("V", $A^a$, 0),)), (legs: (("V", $A^b$, 0),)), (legs: (("S", $H^dagger$, -1),))),
      (("V", none, 0), ("V", none, 0), ("V", none, 0), ("S", none, -1)), caption: [(ii)]),
    sm(((legs: (("S", $H$, 1),)), (legs: (("V", $A^a$, 0),)), (legs: (("S", [#h(2mm)$H^dagger$], -1),)), (legs: (("V", $A^b$, 0),))),
      (("S", none, 1), ("S", none, 1), ("V", none, 0), ("V", none, 0)), caption: [(iii)]),
    sm(((legs: (("S", $H$, 1), ("S", $H^dagger$, -1))), (legs: (("V", $A^a$, 0),)), (legs: (("V", $A^b$, 0),))),
      (("V", none, 0), ("V", none, 0), ("V", none, 0)), caption: [(iv)]),
  ),
  diagrow(
    sm(((legs: (("S", $H$, 1), ("S", $H^dagger$, -1))), (legs: (("V", $A^a$, 0),)), (legs: (("V", $A^b$, 0),))),
      (("S", none, 1), ("S", none, 1), ("S", none, 1)), caption: [(v)]),
    sm(((legs: (("S", $H$, 1), ("V", $A^a$, 0))), (legs: (("V", $A^b$, 0),)), (legs: (("S", $H^dagger$, -1),))),
      (("S", none, 1), ("S", none, 1), ("V", none, 0)), caption: [(vi)]),
    sm(((legs: (("S", $H$, 1), ("V", $A^a$, 0))), (legs: (("V", $A^b$, 0),)), (legs: (("S", $H^dagger$, -1),))),
      (("V", none, 0), ("V", none, 0), ("S", none, -1)), caption: [(vii)]),
    sm(((legs: (("S", $H$, 1),)), (legs: (("V", [#h(2mm)$A^a$], 0), ("V", [$A^b$#h(2mm)], 0))), (legs: (("S", $H^dagger$, -1),))),
      (("S", none, 1), ("S", none, 1), ("V", none, 0)), caption: [(viii)]),
    sm(((legs: (("S", $H$, 1),)), (legs: (("V", [#h(2mm)$A^a$], 0), ("V", [$A^b$#h(2mm)], 0))), (legs: (("S", $H^dagger$, -1),))),
      (("V", none, 0), ("V", none, 0), ("S", none, -1)), caption: [(ix)]),
  ),
  diagrow(
    sm(((legs: (("S", $H$, 1), ("S", $H^dagger$, -1))), (legs: (("V", $A^a$, 0), ("V", $A^b$, 0)))),
      (("S", none, 1), ("S", none, 1)), caption: [(x)]),
    sm(((legs: (("S", $H$, 1), ("S", $H^dagger$, -1))), (legs: (("V", $A^a$, 0), ("V", $A^b$, 0)))),
      (("V", none, 0), ("V", none, 0)), caption: [(xi)]),
    sm(((legs: (("S", $H$, 1), ("V", $A^a$, 0))), (legs: (("S", $H^dagger$, -1), ("V", $A^b$, 0)))),
      (("S", none, 1), ("V", none, 0)), caption: [(xii)]),
  )),
  caption: [Topologies of $H^dagger H A$ (rows one and two) and $H^dagger H A A$ (rows three to five). Dashed lines with arrows: Higgs doublet; wavy: $W$ or $B$; solid: fermions. In (F) and (F4) the gauge bosons may sit on either fermion line; the mirror images of (c), (vi), (vii), (xii) (seagull at the $H^dagger$ vertex, or $a <-> b$) are separate diagrams.],
) <fig:g0_SSV_SSVV-diagrams>

#figure(
  {
  set par(justify: false)
  table(
    columns: 5,
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [symm.], [sign], [group factor (matrix in $i j$)]),
    [(F)$psi'$, (F)$psi$], [$psi, psi'$ for $y_d, y_e, y_u$; $A^a$ on $psi'$ or $psi$], [$1$], [$-1$], [$y^2 N^c_psi thin tl(a)_(psi')$, $y^2 N^c_psi thin tl(a)_psi$],
    [(a)], [$H, H, A^c$], [$1$], [$+1$], [$Tg^c Tg^a Tg^c = (C_H - 1/2 Ca(a)) Tg^a$],
    [(b)], [$H, A^c, A^d$], [$1$], [$+1$], [$fg^(c d a) Tg^d Tg^c = - i/2 Ca(a) Tg^a$],
    [(c), (c$'$)], [$H, A^c$], [$1$], [$+1$], [$Tg^c {Tg^a, Tg^c} = {Tg^a, Tg^c} Tg^c = (2 C_H - 1/2 Ca(a)) Tg^a$],
    [(d)], [$A^c, A^d$], [$1\/2$], [$+1$], [$fg^(c d a) {Tg^c, Tg^d} = 0$],
    [(e)], [$H, H$], [$1$], [$+1$], [$lambda (delta_(i j) Tr Tg^a + Tg^a)$; integrand odd],
    [(F4)], [as (F); $a, b$ on $psi' psi'$, $psi psi$, $psi' psi$, $psi psi'$], [$1$], [$-1$], [$y^2 N^c_psi thin tl(a) tl(b)$ products, @eq:g0_SSV_SSVV-F4],
    [(i)], [$H, H, H, A^c$], [$1$], [$+1$], [$Tg^c Tg^b Tg^a Tg^c$ $+ (a <-> b)$],
    [(ii)], [$H, A^c, A^e, A^d$], [$1$], [$+1$], [$fg^(c a e) fg^(e b d) Tg^d Tg^c$ $+ (a <-> b)$],
    [(iii)], [$H, H, A^c, A^d$], [$1$], [$+1$], [$fg^(c b d) Tg^d Tg^a Tg^c$ $+ (a <-> b)$],
    [(iv)], [$A^c, A^e, A^d$], [$1$], [$+1$], [$fg^(c a e) fg^(e b d) {Tg^c, Tg^d}$],
    [(v)], [$H, H, H$], [$1$], [$+1$], [$lambda [delta_(i j) Tr(Tg^b Tg^a) + Tg^b Tg^a] + (a <-> b)$],
    [(vi)], [$H, H, A^c$], [$1$], [$+1$], [$Tg^c Tg^b {Tg^a, Tg^c}$, ${Tg^a, Tg^c} Tg^b Tg^c$ $+ (a <-> b)$],
    [(vii)], [$H, A^c, A^d$], [$1$], [$+1$], [$fg^(c b d) Tg^d {Tg^a, Tg^c}$, $fg^(d b c) {Tg^a, Tg^c} Tg^d$ $+ (a <-> b)$],
    [(viii)], [$H, H, A^c$], [$1$], [$+1$], [$Tg^c {Tg^a, Tg^b} Tg^c$],
    [(ix)], [$H, A^c, A^d$], [$1$], [$+1$], [$fg^(x y e) fg^(z w e) Tg^d Tg^c$ (three pairings of @eq:V-vvvv)],
    [(x)], [$H, H$], [$1$], [$+1$], [$lambda [delta_(i j) Tr {Tg^a, Tg^b} + {Tg^a, Tg^b}]$],
    [(xi)], [$A^c, A^d$], [$1\/2$], [$+1$], [$fg^(x y e) fg^(z w e) {Tg^c, Tg^d}$ (three pairings)],
    [(xii)], [$H, A^c$], [$1$], [$+1$], [${Tg^b, Tg^c}{Tg^a, Tg^c}$ $+ (a <-> b)$],
  )},
  caption: [Field assignments, symmetry factors, signs and group factors ($Tg = Tg_H$ on the Higgs line; internal $c, d, e$ run over $W^I$ and $B$). Matrices are ordered against the Higgs arrow (the $H^dagger_i$ end on the left). $N^c_psi$ is the colour multiplicity of the loop: $3$ for the quark loops, $1$ for the lepton loop. In (ix) and (xi) $(x, y, z, w)$ is $(a, b, c, d)$, $(a, c, b, d)$, $(a, d, b, c)$ for the three terms of @eq:V-vvvv.],
) <tab:g0_SSV_SSVV-diagrams>

The fermion-line factors $tl(a)_(psi')$, $tl(a)_psi$ are the generators of the two lines carried onto the indices $(i, j)$ by the Yukawa vertices @eq:V-yuk. For the $Q$--$D$ and $L$--$E$ loops the doublet is $psi'$: $tl(a)_(psi') = Tg^a_(psi')$, $tl(a)_psi = Tg^a_psi bb(1)$. For the $Q$--$U$ loop the doublet is $psi = Q$ and the two vertices carry $epsilon_(l j)$ and $epsilon_(l' i)$: $tl(a)_(Q) = epsilon^T (Tg^a_Q)^T epsilon$, i.e. $tl(I)_Q = - gw tau^I \/2$, $tl(B)_Q = gy\/6$, and $tl(a)_U = Tg^a_U bb(1)$. Hypercharge and isospin conservation at the Yukawa vertices give in all three loops
$ tl(a)_(psi') - tl(a)_psi = Tg^a_H . $ <eq:g0_SSV_SSVV-tdiff>
Gluons: in (F) the colour factor is $Tr_"colour" T^A = 0$; the Higgs has no gluon coupling, so no bosonic diagram has an external gluon.

*Integrands.* Three-point functions. In (F)$psi'$ the line $psi$ carries $k$ from the $H^dagger_i$ vertex to the $H_j$ vertex, $psi'$ carries $k + p$ and $k + p'$; in (F)$psi$ $psi$ carries $k$ (into the $H_j$ vertex) and $k - q$, $psi'$ carries $k + p$. With @eq:sandwich every projector is absorbed and the numerators are four-dimensional (@eq:phys-prop); the loop factor $(-1)(-i y)^2 i^4 = y^2$ and the single remaining projector give
$ G_"(F)"^(psi') = y^2 N^c_psi tl(a)_(psi') intk (tr[slb(k) (slb(k) + slb(p)') gbar^mu (slb(k) + slb(p)) P])/(k^2 (k + p')^2 (k + p)^2), \
  G_"(F)"^psi = y^2 N^c_psi tl(a)_psi intk (tr[slb(k) gbar^mu (slb(k) - slb(q)) (slb(k) + slb(p)) P])/(k^2 (k - q)^2 (k + p)^2), $ <eq:g0_SSV_SSVV-int-F>
with $P = PR$ ($y_d$, $y_e$) or $PL$ ($y_u$). In (a) the gauge line carries $k$ from the $H_j$ vertex to the $H^dagger_i$ vertex and the Higgs line $p - k$, $p' - k$; in (b) the Higgs line carries $p - k$, the gauge lines $k$ and $k + q$; in (c) the gauge line carries $k$ and the Higgs line $p' - k$:
$ G_"(a)" = Tg^c Tg^a Tg^c intk ((2p - k) dot (2 p' - k) thin (p + p' - 2k)^mu)/(k^2 ((p - k)^2 - m^2)((p' - k)^2 - m^2)), $ <eq:g0_SSV_SSVV-int-a>
$ G_"(b)" = i fg^(c d a) Tg^d Tg^c intk ((2p - k)^rho (p + p' - k)^sigma cal(V)_(rho sigma)^(space space mu))/(k^2 (k + q)^2 ((p - k)^2 - m^2)), \
  cal(V)_(rho sigma mu) = g_(rho sigma) (2k + q)_mu - g_(sigma mu) (k + 2q)_rho + g_(mu rho) (q - k)_sigma , $ <eq:g0_SSV_SSVV-int-b>
$ G_"(c)" = - Tg^c {Tg^a, Tg^c} intk ((2p' - k)^mu)/(k^2 ((p' - k)^2 - m^2)), quad
  G_"(c')" = - {Tg^a, Tg^c} Tg^c intk ((2p - k)^mu)/(k^2 ((p - k)^2 - m^2)) . $ <eq:g0_SSV_SSVV-int-c>
$cal(V)_(rho sigma mu)$ is the bracket of @eq:V-vvv for the legs $(A^c_rho, k)$, $(A^d_sigma, -k - q)$, $(A^a_mu, q)$ (a symbol of this subsection, not a vertex function). (e) is $prop intk (2k + q)^mu \/ ((k^2 - m^2)((k + q)^2 - m^2)) = 0$ (odd under $k -> -k - q$).

Four-point functions. Since $dg(G) = 0$, $Top(0)$ sets all external momenta and $m$ to zero and every denominator to $k^2 - MUV^2$; only the value of the numerator at $p = q_1 = q_2 = 0$ is needed. With $k$ on every gauge line (in the direction from the vertex nearest $H_j$) and $-k$ on every Higgs line, the vertices are $i (-k)^rho Tg^c$ at the $H_j$ and $H^dagger_i$ vertices, $i (-2k)^mu Tg^a$ for an external $A^a$ on the Higgs line, $i (+2k)^mu$ on a closed Higgs loop, and the bracket of @eq:V-vvv with $(p_1, p_2, p_3) = (k, 0, -k)$, written $cal(V)_0$:
$ cal(V)_0^(rho' mu alpha) = g^(rho' mu) k^alpha + g^(mu alpha) k^(rho') - 2 g^(alpha rho') k^mu , quad k_(rho') cal(V)_0^(rho' mu alpha) = g^(mu alpha) k^2 - k^mu k^alpha . $ <eq:g0_SSV_SSVV-vvv0>
The numerators $n_X$ of the topologies $X$ (group factor and coupling constants stripped, all factors $i$ and $-i$ included) are
$ & n_"(i)" = - 4 k^2 k^mu k^nu, quad
  n_"(ii)" = (g^(mu alpha) k^2 - k^mu k^alpha)(g_alpha^(space nu) k^2 - k_alpha k^nu), quad
  n_"(iii)" = 2 i k^mu thin k_(rho') k_(sigma') cal(V)_0^(rho' nu sigma') = 0, \
  & n_"(iv)" = - g_(rho' sigma') cal(V)_0^(rho' mu alpha) cal(V)_(0 alpha)^(space nu sigma') = - 2 g^(mu nu) k^2 - (4 d - 6) k^mu k^nu, quad
  n_"(v)" = 8 lambda k^mu k^nu, quad
  n_"(vi)" = 2 k^mu k^nu, \
  & n_"(vii)" = - i (g^(mu nu) k^2 - k^mu k^nu), quad
  n_"(viii)" = g^(mu nu) k^2, quad
  n_"(ix)" = k_(rho') k_(sigma') thin Pi^(mu nu rho' sigma'), quad
  n_"(x)" = -2 lambda g^(mu nu), \
  & n_"(xi)" = - 1/2 g_(rho' sigma') Pi^(mu nu rho' sigma'), quad
  n_"(xii)" = - g^(mu nu) , $ <eq:g0_SSV_SSVV-n4>
each divided by $DM(ell_X)$, $ell_X$ the number of propagators of $X$. In (ix) and (xi) $Pi^(mu nu rho' sigma')$ (a symbol of this subsection) is the bracket of @eq:V-vvvv (the factor $-i$ is in the prefactor), one pairing at a time, for the legs $(A^a_mu, A^b_nu, A^c_(rho'), A^d_(sigma'))$; per pairing, $k_(rho') k_(sigma') (g^(mu nu) g^(rho' sigma') - g^(mu sigma') g^(nu rho')) = g^(mu nu) k^2 - k^mu k^nu$ for the pairings $(a c)(b d)$, $(a d)(b c)$, and $0$ for $(a b)(c d)$; $g_(rho' sigma')$ gives $(d - 1) g^(mu nu)$, $(d - 1) g^(mu nu)$ and $0$. The fermion box (F4) has the loop factor $(-1)(-i y)^2 i^2 i^4 = -y^2$ (fermion loop, two Yukawa vertices, two gauge vertices, four propagators), and with the Dirac algebra of (F) its numerators are $-y^2 N^c_psi thin 1/2 tr[slb(k) slb(k) gbar^nu slb(k) gbar^mu slb(k)]$ when both gauge vertices sit on one line (either line) and $-y^2 N^c_psi thin 1/2 tr[slb(k) gbar^mu slb(k) slb(k) gbar^nu slb(k)]$ when they sit on different lines, with $mu$ and $nu$ exchanged for the other order of the two vertices. The factor $1\/2$ is $tr[dots P] = 1/2 tr[dots]$ (the $g5$ part has only one independent vector and vanishes). The two traces are
$ 1/2 tr[slb(k) slb(k) gbar^nu slb(k) gbar^mu slb(k)] = 2 bar(k)^2 (2 bar(k)^mu bar(k)^nu - bar(k)^2 gmb^(mu nu)), quad
  1/2 tr[slb(k) gbar^mu slb(k) slb(k) gbar^nu slb(k)] = 2 (bar(k)^2)^2 gmb^(mu nu) . $ <eq:g0_SSV_SSVV-F4-traces>

*$Rop$-operation.* At one loop $Rbar G = G$ and $Z(G) = - Kop Top(dg(G)) G$ (@eq:m-Z1). For the three-point functions, $Top(1)$ keeps the numerator to first order in the external momenta and
$ Top(1) 1/((k + v_1)^2 - m^2) dots 1/((k + v_ell)^2 - m^2) = 1/DM(ell) - (2 k dot (v_1 + dots + v_ell))/DM(ell+1) $ <eq:g0_SSV_SSVV-T1>
where $v_1, dots, v_ell$ are the external momenta (or sums of them) carried by the $ell$ propagators in addition to $k$. $m^2$ and $MUV^2$ have degree two and do not enter $Top(1)$ or $Top(0)$. The integrals needed are $Kop intk (k^2)^a \/ DM(a+2) = ilf \/ep$ (@eq:logpole) with @eq:tensor-red and @eq:tensor-red-split at $d = 4$, among them $bar(k)^2 bar(k)^mu k^nu -> (k^2)^2 thin 6 gmb^(mu nu) \/ (d (d+2))$ and $(k^2) k^mu k^nu -> (k^2)^2 g^(mu nu) \/ d$.

(F) In $G^(psi')_"(F)"$, $1/2 tr[slb(a) slb(b) gbar^mu slb(c)] = 2 [(bar(a) dot bar(b)) bar(c)^mu - bar(a)^mu (bar(b) dot bar(c)) + (bar(a) dot bar(c)) bar(b)^mu]$; the terms with $g5$ are $prop epsilon^(mu alpha beta gamma)$ contracted with $k$ and at most one external momentum and vanish after reduction. The numerator is $2 bar(k)^2 bar(k)^mu + 2 bar(k)^2 (bar(p) + bar(p)')^mu + O(p^2)$ (the terms $(bar(k) dot bar(p)) bar(k)^mu$ cancel). The pole is
$ 2 (bar(p) + bar(p)')^mu - 4 thin 6/24 thin (bar(p) + bar(p)')^mu = (bar(p) + bar(p)')^mu . $
In $G^psi_"(F)"$ the numerator is $2 bar(k)^2 bar(k)^mu + 4 (bar(k) dot bar(p)) bar(k)^mu - 2 bar(k)^2 bar(p)'^mu + O(p^2)$ and the denominators give $-2 k dot (p - q)$; the pole is $(bar(p) - 2 bar(p)') - (2 bar(p) - bar(p)') = - (bar(p) + bar(p)')^mu$. Only $bar(k) dot bar(p)$ survives from $k dot p = bar(k) dot bar(p) + hat(k) dot hat(p)$ (a single $hat(k)$ integrates to zero), so no hatted structure and no $eta$ arise.

(a) Numerator $-2 k^2 k^mu + k^2 (p + p')^mu + 4 (k dot (p + p')) k^mu + O(p^2)$, denominators $+2 k dot (p + p')$: pole $(1 + 1 - 1)(p + p')^mu$.

(b) $(2p - k)^rho (p + p' - k)^sigma cal(V)_(rho sigma mu)$ has no $k^3$ term; at first order it is $k^2 q_mu - (k dot q) k_mu - 2 [(k dot p) k_mu - k^2 p_mu] - [(k dot (p + p')) k_mu - k^2 (p + p')_mu]$, with pole $3/4 q + 3/2 p + 3/4 (p + p') = 3/2 (p + p')$.

(c) $-(2 p' - k)^mu [1 + 2 k dot p' \/ (k^2 - MUV^2)]$ gives $-(2 - 1/2) p'^mu$; (c$'$) gives $-3/2 thin p^mu$.

Per topology, in units of $ilf\/ep$ and with the group factors of @tab:g0_SSV_SSVV-diagrams,
$ & N_0^"(F)" = y^2 N^c_psi (tl(a)_(psi') - tl(a)_psi)(bar(p) + bar(p)')^mu = y^2 N^c_psi (bar(p) + bar(p)')^mu Tg^a_H, quad
  N_0^"(a)" = (C_H - 1/2 Ca(a)) (p + p')^mu Tg^a_H, \
  & N_0^"(b)" = 3/4 Ca(a) (p + p')^mu Tg^a_H, quad
  N_0^"(c)" + N_0^"(c')" = - 3/2 (2 C_H - 1/2 Ca(a)) (p + p')^mu Tg^a_H, quad N_0^"(d)" = N_0^"(e)" = 0, $ <eq:g0_SSV_SSVV-N3>
and $N_1 = 0$ for all of them.

Four-point functions. With $k^mu k^nu -> k^2 g^(mu nu)\/4$ and @eq:logpole the bosonic numerators @eq:g0_SSV_SSVV-n4 give $N_0^X = r_X thin g^(mu nu) times$ (group factor of $X$), with
$ r_"(i)" = -1, quad r_"(ii)" = 3/4, quad r_"(iii)" = 0, quad r_"(iv)" = -9/2, quad r_"(vi)" = 1/2, \
  r_"(vii)" = -3/4 i, quad r_"(viii)" = 1, quad r_"(ix)" = 3/4, quad r_"(xi)" = -3/2, quad r_"(xii)" = -1, $
per diagram, and per pairing of $Pi$ in (ix), (xi) (@eq:V-vvvv), and $r_"(v)" = 2 lambda$, $r_"(x)" = -2 lambda$. In (iv) $4d - 6 = 10$ at $ep = 0$. The two $lambda$ diagrams have the same group factor $delta_(i j) Tr{Tg^a, Tg^b} + {Tg^a, Tg^b}$ and cancel. With $bar(k)^2 bar(k)^mu bar(k)^nu -> (k^2)^2 gmb^(mu nu)\/4$ and $(bar(k)^2)^2 -> (k^2)^2$ the traces @eq:g0_SSV_SSVV-F4-traces become $-(k^2)^2 gmb^(mu nu)$ and $2 (k^2)^2 gmb^(mu nu)$; with the loop factor $-y^2 N^c_psi$ the fermion box gives $+gmb^(mu nu)$ for both vertices on one line and $-2 gmb^(mu nu)$ for one on each line, so that
$ N_0^"(F4)" = y^2 N^c_psi [ {tl(a)_(psi'), tl(b)_(psi')} + {tl(a)_psi, tl(b)_psi} - 2 (tl(a)_(psi') tl(b)_psi + tl(b)_(psi') tl(a)_psi) ] gmb^(mu nu) \
  = y^2 N^c_psi {Tg^a_H, Tg^b_H} gmb^(mu nu) $ <eq:g0_SSV_SSVV-F4>
by @eq:g0_SSV_SSVV-tdiff (the factor $t$ of the singlet line is a multiple of $bb(1)$). For external gluons both lines carry $gs T^A$ on one colour trace, $tl(A)_(psi') = tl(A)_psi$, and the bracket vanishes: $N_0 = 0$ for $H^dagger H G G$; $H^dagger H G$, $H^dagger H G W$ and $H^dagger H G B$ vanish identically ($Tr_"colour" T^A = 0$).

The bosonic group factors are reduced with $[Tg^a, Tg^b] = i fg^(a b c) Tg^c$ and the Jacobi identity to three structures,
$ N_0^X = g^(mu nu) [ c^((1))_X thin C_H {Tg^a, Tg^b} + c^((2))_X thin (Ca(a) + Ca(b)) {Tg^a, Tg^b} + c^((3))_X thin Ups ], \
  Ups = sum_c Tg^c_H {Tg^a_H, Tg^b_H} Tg^c_H , $ <eq:g0_SSV_SSVV-xyz>
with the coefficients of @tab:g0_SSV_SSVV-xyz. Two rules give the reduction ($Tg = Tg_H$, sums over repeated gauge indices). Every structure constant next to a generator is a commutator, $fg^(x y e) Tg^e = -i [Tg^x, Tg^y]$. A summed pair $Tg^c dots Tg^c$ is reduced by $Tg^c Tg^c = C_H bb(1)$ when the two are adjacent and, when one generator stands between them, by
$ Tg^c Tg^a Tg^c = (C_H - 1/2 Ca(a)) Tg^a , $ <eq:g0_SSV_SSVV-TaT>
which follows from $Tg^c Tg^a = Tg^a Tg^c + i fg^(c a e) Tg^e$ and @eq:invariants. Of a product of four generators, what remains is $#Ups$. Two rows of @tab:g0_SSV_SSVV-xyz as examples. In (xii) the group factor ${Tg^b, Tg^c}{Tg^a, Tg^c}$ expands into four products, three of which reduce with @eq:g0_SSV_SSVV-TaT and $Tg^c Tg^c = C_H$; in (ii) the two structure constants become $fg^(e b d) Tg^d = -i [Tg^e, Tg^b]$ and $fg^(c a e) Tg^c = -i [Tg^a, Tg^e]$:
$ "(xii)": quad {Tg^b, Tg^c}{Tg^a, Tg^c} + (a <-> b) = (3 C_H - 1/2 (Ca(a) + Ca(b))) {Tg^a, Tg^b} + #Ups, \
  "(ii)": quad fg^(c a e) fg^(e b d) Tg^d Tg^c + (a <-> b) = - [Tg^e, Tg^b][Tg^a, Tg^e] + (a <-> b) \
  = (C_H - 1/2 (Ca(a) + Ca(b))) {Tg^a, Tg^b} - #Ups . $ <eq:g0_SSV_SSVV-red-ex>
Multiplied by $r_"(xii)" = -1$ and $r_"(ii)" = 3/4$ these give the columns (xii), $(c^((1)), c^((2)), c^((3))) = (-3, 1/2, -1)$, and (ii), $(3/4, -3/8, -3/4)$, of @tab:g0_SSV_SSVV-xyz. The other columns follow in the same way; for $a, b$ both in $B$ all structure constants vanish and $Ca(B) = 0$.

#figure(
  table(
    columns: 12,
    align: center,
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([], [(i)], [(ii)], [(iii)], [(iv)], [(v)+(x)], [(vi)], [(vii)], [(viii)], [(ix)], [(xi)], [(xii)]),
    [$c^((1))_X$], [$0$], [$3/4$], [$0$], [$-9/2$], [$0$], [$1$], [$3/2$], [$0$], [$-3/4$], [$3$], [$-3$],
    [$c^((2))_X$], [$0$], [$-3/8$], [$0$], [$9/4$], [$0$], [$-1/4$], [$0$], [$0$], [$3/8$], [$-3/2$], [$1/2$],
    [$c^((3))_X$], [$-1$], [$-3/4$], [$0$], [$9/2$], [$0$], [$1$], [$-3/2$], [$1$], [$3/4$], [$-3$], [$-1$],
  ),
  caption: [Coefficients of @eq:g0_SSV_SSVV-xyz per topology (all field assignments and $a <-> b$ included). The sums are $sum_X c^((1))_X = -2$, $sum_X c^((2))_X = 1$, $sum_X c^((3))_X = 0$.],
) <tab:g0_SSV_SSVV-xyz>

(i) and (viii) cancel; $#Ups$ drops out of the sum. No diagram produces a hatted structure or $eta$: the bosonic numerators are #box[$d$-dimensional] and never split, and in the fermion loops at most one $hat(k) dot hat(r)$ appears, which integrates to zero. Hence $N_1 = 0$ and the external indices of the fermion parts are barred.

*Result.*
#result[
  $ Gamma^((1))_(H^dagger_i H_j A^a_mu) |_"pole" = lf 1/ep N_0, quad
    N_0 = [ Yt (bar(p) + bar(p)')^mu + (Ca(a) - 2 C_H)(p + p')^mu ] (Tg^a_H)_(i j), $ <eq:g0_SSV_SSVV-pole3>
  $ Gamma^((1))_(H^dagger_i H_j A^a_mu A^b_nu) |_"pole" = lf 1/ep N_0, quad
    N_0 = [ Yt thin gmb^(mu nu) + (Ca(a) + Ca(b) - 2 C_H) g^(mu nu) ] {Tg^a_H, Tg^b_H}_(i j), $ <eq:g0_SSV_SSVV-pole4>
  with $Yt = 3 y_u^2 + 3 y_d^2 + y_e^2 = sum_"loops" N^c_psi y^2$ (@sec:st-solution). The bosonic parts are #box[$d$-dimensional] and covariant; the fermion loops produce no hatted structure. Relative to the #box[$d$-dimensional] tree structures the evanescent poles are $-Yt (hat(p) + hat(p)')^mu Tg^a_H$ and $-Yt thin gmh^(mu nu) {Tg^a_H, Tg^b_H}$, the vertex functions of $Yt thin tilde(E)_H$. The functions with an external gluon have no pole. As a local functional (terms with one and two gauge fields, $c_r in {gw, gy}$ the coupling of the factor $r$):
  $ delta S^((1)) supset - lf 1/ep integral dif^d x thin [ & Yt thin (bar(D)_mu H)^dagger bar(D)^mu H \
    & + (sum_r Ca(r) c_r partial/(partial c_r) - 2 C_H) (D_mu H)^dagger D^mu H ] |_(A, thin A A) , $ <eq:g0_SSV_SSVV-ct>
  where $bar(D)_mu = bar(partial)_mu - i gmb_mu^(space nu) Am_nu$ has a barred index and couples only to $bar(A)$. In the basis of @sec:bases ($Q_J$, $Q_(W W)$, $Q_(B B)$, $Q_(W B)$ continued with #box[$d$-dimensional] indices as in $S_0$, and $tilde(E)_H$ of @eq:evanescent, of which only the parts with gauge fields enter here):
  $ delta S^((1)) supset - lf 1/ep [ & - (Yt + Ca(W) - 2 C_H) Q^((W))_J - (Yt - 2 C_H) Q^((B))_J \
    & + (Yt + 2 Ca(W) - 2 C_H) Q_(W W) + (Yt - 2 C_H) Q_(B B) \
    & + (Yt + Ca(W) - 2 C_H) Q_(W B) + Yt thin tilde(E)_H |_(A, thin A A) ] , $ <eq:g0_SSV_SSVV-ct-basis>
  and $Q_(G G)$ does not appear. With $delta Z_H = lf (2 C_H - Yt)\/ep$ from @sec:g0_SS this is $integral [delta Z_H + sum_r (delta Z_(c_r) + 1/2 delta Z^((r))_A) c_r partial_(c_r)] (D H)^dagger D H |_(A, A A)$ plus the evanescent term, with $delta Z_(c_r)$ and $delta Z^((r))_A$ the coupling and gauge-field renormalisation constants of the factor $r$ (@sec:st-solution) and
  $ delta Z_(c_r) + 1/2 delta Z^((r))_A = - lf Ca(r)/ep, quad "i.e." quad - lf (2 gw^2)/ep "for" r = W, quad 0 "for" r = B . $ <eq:g0_SSV_SSVV-Z>
  The bosonic loops multiply the #box[$d$-dimensional] operators, the fermion loops the barred ones; the evanescent coefficients $-lf Yt\/ep$ of $tilde(E)_H$ and $-2/3 lf Yt\/ep$ of $E_H$ (@sec:g0_SS) are those of the assembled counterterm of @sec:st-solution.
]

#figure(
  table(
    columns: 4,
    align: (left, center, center, center),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([function], [tree], [coefficient of $(p + p')^mu$ or $g^(mu nu)$ (#box[$d$-dim.])], [coeff. of $(bar(p) + bar(p)')^mu$ or $gmb^(mu nu)$]),
    [$H^dagger H W^I_mu$], [$gw tau^I\/2$], [$1/2 gw^2 - 1/2 gy^2$], [$3 y_u^2 + 3 y_d^2 + y_e^2$],
    [$H^dagger H B_mu$], [$gy\/2$], [$-3/2 gw^2 - 1/2 gy^2$], [$3 y_u^2 + 3 y_d^2 + y_e^2$],
    [$H^dagger H W^I_mu W^J_nu$], [$gw^2 delta^(I J)\/2$], [$5/2 gw^2 - 1/2 gy^2$], [$3 y_u^2 + 3 y_d^2 + y_e^2$],
    [$H^dagger H W^I_mu B_nu$], [$gw gy tau^I\/2$], [$1/2 gw^2 - 1/2 gy^2$], [$3 y_u^2 + 3 y_d^2 + y_e^2$],
    [$H^dagger H B_mu B_nu$], [$gy^2\/2$], [$-3/2 gw^2 - 1/2 gy^2$], [$3 y_u^2 + 3 y_d^2 + y_e^2$],
    [$H^dagger H G$, $H^dagger H G G$, $H^dagger H G W$, $H^dagger H G B$], [], [$0$], [$0$],
  ),
  caption: [Standard-Model values: $N_0$ = tree structure ($Tg^a_H$ resp. ${Tg^a_H, Tg^b_H}$, column 2) times the coefficients in columns 3 and 4. Per topology for $W W$: (i) $-C_H$, (ii) $-3/2 gw^2$, (iv) $9 gw^2$, (vi) $1/2 (gw^2 + gy^2)$, (vii) $0$, (viii) $C_H$, (ix) $3/2 gw^2$, (xi) $-6 gw^2$, (xii) $-(gw^2 + gy^2)$; for $W B$: (i) $(gw^2 - gy^2)\/4$, (vi) $gy^2\/2$, (vii) $3/2 gw^2$, (viii) $-(gw^2 - gy^2)\/4$, (xii) $-(gw^2 + gy^2)$; for $B B$: (i) $-C_H$, (vi) $2 C_H$, (viii) $C_H$, (xii) $-4 C_H$; all others zero.],
) <tab:g0_SSV_SSVV-sm>

Coefficients in $16 pi^2 ep thin delta S^((1))$ (@eq:g0_SSV_SSVV-ct-basis): $Q^((W))_J$: $Yt + 1/2 gw^2 - 1/2 gy^2$; $Q^((B))_J$: $Yt - 3/2 gw^2 - 1/2 gy^2$; $Q_(W W)$: $-(Yt + 5/2 gw^2 - 1/2 gy^2)$; $Q_(B B)$: $-(Yt - 3/2 gw^2 - 1/2 gy^2)$; $Q_(W B)$: $-(Yt + 1/2 gw^2 - 1/2 gy^2)$; $tilde(E)_H |_(A, thin A A)$: $-Yt$; $Q_(G G)$: $0$.

#check[*Independence of $MUV$ and $m$.* $MUV^2$ and $m^2$ have degree two and do not enter $Top(1)$ and $Top(0)$; $N_0$ is independent of both identically.]

#check[*Group theory.* The sum over the bosonic topologies is proportional to the tree structure, as required by rigid invariance, and only through $C_H$ and $Ca(a)$. The reduction @eq:g0_SSV_SSVV-xyz and the totals of @eq:g0_SSV_SSVV-pole3 and @eq:g0_SSV_SSVV-pole4 were also evaluated numerically with explicit matrices for the Standard-Model doublet ($a, b in {W^I, B}$) and, as an independent test of the invariant form, for a colour-triplet scalar of $"SU"(3) times "U"(1)$, where $#Ups$ is not proportional to ${Tg^a, Tg^b}$ diagram by diagram. In both cases the sum is $(Ca(a) + Ca(b) - 2 C_H) {Tg^a, Tg^b}$.]

#check[*Known limits and Ward identity.* $Pphys N_0$ is the result with an anticommuting $g5$: the counterterm vertices are $delta Z_H + delta Z_(c_r) + 1/2 delta Z^((r))_A$ times the tree vertex for $H^dagger H A^((r))$ and $delta Z_H + sum_(t = r, s) (delta Z_(c_t) + 1/2 delta Z^((t))_A)$ times the tree vertex for $H^dagger H A^((r)) A^((s))$, with $delta Z_H$ from @sec:g0_SS. These are the standard Feynman-gauge $overline("MS")$ constants: $delta Z_H = lf (2 C_H - Yt)\/ep$, and $delta Z_(c_r) + 1/2 delta Z^((r))_A = lf [(-11/6 + 5/6) Ca(r) + (1/3 - 1/3) sum_f S^((r))_f + (1/6 - 1/6) S^((r))_H] \/ ep = - lf Ca(r)\/ep$ with the general forms in the table of constants of @sec:st-solution: the matter loops cancel in this combination. It is $Z_(c_r) (Z^((r))_A)^(1\/2) - 1$, which also relates the gauge vertices of @sec:g0_VVV and @sec:g0_FFV to the field constants (@sec:g0_FFV). For $B$ the abelian Ward identity holds for the #box[$d$-dimensional] bosonic parts: $q_mu (-2 C_H)(p + p')^mu Tg^B = -2 C_H (p'^2 - p^2) Tg^B$, the difference of the $p$-dependent bosonic poles $-2 C_H p^2$ of $Gamma^((1))_(H^dagger H) (p')$ and $Gamma^((1))_(H^dagger H) (p)$ (@sec:g0_SS), times $Tg^B$.]

#check[*Evanescent structure and relation to @sec:g0_SS and @sec:g1_cSS_cSSV.* The fermion-loop parts of the $H^dagger H$, $H^dagger H A$ and $H^dagger H A A$ poles, $Yt (bar(p)^2 + 2/3 hat(p)^2)$, $Yt (bar(p) + bar(p)')^mu Tg^a$ and $Yt gmb^(mu nu) {Tg^a, Tg^b}$, are the vertex functions of the single operator $Yt [(bar(D) H)^dagger bar(D) H + 2/3 (hat(partial) H)^dagger hat(partial) H] = Yt [(D H)^dagger D H + tilde(E)_H + 2/3 E_H]$. The two-point function fixes only the difference of the coefficients of $E_H$ and $tilde(E)_H$ (the two-point part of $tilde(E)_H$ is $-E_H$); the $H^dagger H A$ function fixes the coefficient $-lf Yt\/ep$ of the $bop_d$-closed $tilde(E)_H$, and hence $-2/3 lf Yt\/ep$ for the non-invariant $E_H$, as in @sec:g0_SS. Consequently $[(Dbrk_1 dot Gamma)^((1))_(gh H^dagger H)]_"div" = - bop_d delta S^((1)) supset 2/3 lf Yt\/ep thin bop_d E_H$, and the divergent part of $(Dbrk_1 dot Gamma)^((1))_(gh H^dagger H A)$ vanishes, since $bop_d E_H$ has no gauge field and $bop_d tilde(E)_H = 0$. Directly: $q_mu Gamma^mu_(H^dagger H B)$ and $Tg^B [Gamma_(H^dagger H) (p') - Gamma_(H^dagger H) (p)]$ differ by the evanescent $-2/3 lf Yt (hat(p)'^2 - hat(p)^2) Tg^B\/ep$, the $gh^B H^dagger H$ structure $prop (hat(p)'^2 - hat(p)^2) Tg^B$ of $bop_d E_H$ (@eq:b-evanescent).]

#general[At one loop the forest is empty and the fermion lines couple only to $bar(A)$, so the fermion-loop poles are those of $(bar(D) H)^dagger bar(D) H$ plus the evanescent $E_H$. At higher orders the bosonic subgraphs stay #box[$d$-dimensional] and BRST invariant, while $delta S^((1))$ must be inserted with its complete #box[$d$-dimensional] form: $Yt thin tilde(E)_H$ (closed but evanescent) and $2/3 Yt thin E_H$ (not closed) produce, through $hat(partial) -> hat(k)$ and the hatted components of internal gauge lines, factors $eta\/ep^2$ in two-loop graphs, i.e. two-loop poles and evanescent-induced finite parts (@sec:m-beyond).]
