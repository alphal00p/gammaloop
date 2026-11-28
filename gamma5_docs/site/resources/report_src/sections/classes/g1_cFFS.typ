#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Breaking of $omega overline(Psi) Psi H$ <sec:g1_cFFS>

External legs $gh^a (q)$, $Psi_(f_2) (p)$ and $phi (p_phi)$ incoming, $overline(Psi)_(f_1)$ outgoing with $p' = p + q + p_phi$ (incoming $-p'$), for the six Yukawa vertices
$ (f_1, f_2, phi) = (Q, U, H^dagger), quad (Q, D, H), quad (L, E, H), \
  "conjugates:" quad (U, Q, H), quad (D, Q, H^dagger), quad (E, L, H^dagger) , $ <eq:g1_cFFS-legs>
and every gauge index $a$. The sum of the diagrams below is $i thin (Dbrk_1 dot Gamma)^((1))_(gh^a overline(Psi)_(f_1) Psi_(f_2) phi_i)$, a matrix in spinor space, written with the ghost as the leftmost Grassmann factor (as in @eq:V-delta). Every graph has four propagators and a numerator of degree four, so its superficial degree of divergence is $dg(G) = 4 + 4 - 2 dot 4 = 0$: $N_0$ and $N_1$ are constants, and $Top(0)$ sets $q, p, p', p_phi$ and $m^2$ to zero and replaces every denominator by $D = 1 \/ (k^2 - MUV^2)$.

The notation is that of @sec:g0_FFS. The Yukawa vertex $overline(Psi)_(f_1) Psi_(f_2) phi_j$ of @eq:V-yuk is $-i thin bold(Y)^j P_y$, with $bold(Y)^j$ the coupling matrix from the indices of $f_2$ to those of $f_1$ and $P_y = P_(overline(c)_(f_1)) = P_(c_(f_2))$ ($PR$ for the first three vertices of @eq:g1_cFFS-legs, $PL$ for the conjugates). The generator on the scalar is $Tg^a_phi = Tg^a_H$ for $phi = H$ and $Tg^a_(H^dagger) = -(Tg^a_H)^*$ for $phi = H^dagger$; rigid invariance of @eq:LYuk is
$ Tg^a_(f_1) bold(Y)^i - bold(Y)^i Tg^a_(f_2) = bold(Y)^j (Tg^a_phi)_(j i) . $ <eq:g1_cFFS-rigid>
$P_(c_1) = P_(c_(f_1))$ and $P_(overline(c)_1) = P_(overline(c)_(f_1)) = P_y$.

*Diagrams.* A one-loop graph has exactly two internal lines at the $Dbrk_1$ vertex. _Internal ghost_ (topologies (a), (b)): $Dbrk_1$ sits at an external fermion leg, its ghost $gh^b$ ends on the antighost of the ghost--gauge vertex @eq:V-ghost that carries the external $gh^a$, and the gauge field $A^d$ of that vertex is absorbed by the other external fermion (a) or by the scalar line (b). Each has two graphs: $Dbrk_1$ at $Psi_(f_2)$ (index 1, drawn) or at $overline(Psi)_(f_1)$ (index 2). _External ghost_ (topologies (c)--(e)): $gh^a$ sits on $Dbrk_1$, whose two fermion legs are internal; the loop closes by gauge exchange between the two fermion legs with the Yukawa vertex of $phi$ on the line (c), by gauge exchange between a fermion leg and the scalar (d), or by Higgs exchange with three Yukawa vertices (e, quark vertices only); $Dbrk_1$ sits on either internal fermion segment. There is no $agh gh H$, $overline(Psi) Psi A H$ or $overline(Psi) Psi H H$ vertex, so these are all one-loop 1PI graphs.

#figure(
  stack(dir: ttb, spacing: 5mm,
  diagrow(
    oneloop(
      ((legs: (("G", $gh^a$, 1),)), (legs: (("F", $overline(Psi)_(f_1)$, -1),)), (legs: (("S", $phi$, 0),)), (legs: (("F", $Psi_(f_2)$, 1),), mark: "delta")),
      (("V", $A^d$, 0), ("F", $f_1$, -1), ("F", $chi_(f_2)$, -1), ("G", $gh^b$, -1)), scale: 0.85, caption: [(a1)]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),)), (legs: (("S", $phi$, 0),)), (legs: (("F", $overline(Psi)_(f_1)$, -1),)), (legs: (("F", $Psi_(f_2)$, 1),), mark: "delta")),
      (("V", $A^d$, 0), ("S", none, 0), ("F", $chi_(f_2)$, -1), ("G", $gh^b$, -1)), scale: 0.85, caption: [(b1)]),
  ),
  diagrow(
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("S", $phi$, 0),)), (legs: (("F", $overline(Psi)_(f_1)$, -1),)), (legs: (("F", $Psi_(f_2)$, 1),))),
      (("F", $f_2$, 1), ("F", $f_1$, 1), ("V", $A^d$, 0), ("F", $f_2$, 1)), scale: 0.85, caption: [(c1)]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("F", $overline(Psi)_(f_1)$, -1),)), (legs: (("S", $phi$, 0),)), (legs: (("F", $Psi_(f_2)$, 1),))),
      (("F", $f_2$, 1), ("S", none, 0), ("V", $A^d$, 0), ("F", $f_2$, 1)), scale: 0.85, caption: [(d1)]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("S", $phi$, 0),)), (legs: (("F", $overline(Psi)_(f_1)$, -1),)), (legs: (("F", $Psi_(f_2)$, 1),))),
      (("F", $f_3$, 1), ("F", $f_4$, 1), ("S", none, 0), ("F", $f_3$, 1)), scale: 0.85, caption: [(e1)]),
  )),
  caption: [Topologies of $(Dbrk_1 dot Gamma)^((1))_(gh overline(Psi) Psi H)$ (crossed circle: @eq:V-delta). The line marked $chi_(f_2)$ is the sterile component next to $Dbrk_1$. (a2), (b2): $Dbrk_1$ at the leg $overline(Psi)_(f_1)$, the gauge field absorbed by $Psi_(f_2)$ (a2) or by the scalar (b2). (c2), (d2), (e2): $Dbrk_1$ on the $f_1$ side of the line; in (d2) the gauge field joins the $overline(Psi)_(f_1)$ leg to the scalar. In (e) $f_3$, $f_4$ are Yukawa partners.],
) <fig:g1_cFFS-diagrams>

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([graph], [internal lines], [sym. factor], [sign], [group factor]),
    [(a1)], [$gh^b$, $A^d$, $Psi_(f_1)$, $chi_(f_2)$], [$1$], [$+1$], [$fg^(b d a) Tg^d_(f_1) bold(Y)^i Tg^b_(f_2)$],
    [(a2)], [$gh^b$, $A^d$, $Psi_(f_2)$, $chi_(f_1)$], [$1$], [$+1$], [$fg^(b d a) Tg^b_(f_1) bold(Y)^i Tg^d_(f_2)$],
    [(b1)], [$gh^b$, $A^d$, $phi$, $chi_(f_2)$], [$1$], [$+1$], [$fg^(b d a) bold(Y)^j (Tg^d_phi)_(j i) Tg^b_(f_2)$],
    [(b2)], [$gh^b$, $A^d$, $phi$, $chi_(f_1)$], [$1$], [$+1$], [$fg^(b d a) Tg^b_(f_1) bold(Y)^j (Tg^d_phi)_(j i)$],
    [(c1)], [$A^d$, $Psi_(f_2)$, $Psi_(f_2)$, $Psi_(f_1)$], [$1$], [$+1$], [$Tg^d_(f_1) bold(Y)^i Tg^a_(f_2) Tg^d_(f_2)$],
    [(c2)], [$A^d$, $Psi_(f_2)$, $Psi_(f_1)$, $Psi_(f_1)$], [$1$], [$+1$], [$Tg^d_(f_1) Tg^a_(f_1) bold(Y)^i Tg^d_(f_2)$],
    [(d1)], [$A^d$, $phi$, $Psi_(f_2)$, $Psi_(f_2)$], [$1$], [$+1$], [$bold(Y)^j (Tg^d_phi)_(j i) Tg^a_(f_2) Tg^d_(f_2)$],
    [(d2)], [$A^d$, $phi$, $Psi_(f_1)$, $Psi_(f_1)$], [$1$], [$+1$], [$Tg^d_(f_1) Tg^a_(f_1) bold(Y)^j (Tg^d_phi)_(j i)$],
    [(e1), (e2)], [$phi'$, $Psi_(f_3)$, $Psi_(f_3)$, $Psi_(f_4)$ / $phi'$, $Psi_(f_3)$, $Psi_(f_4)$, $Psi_(f_4)$], [$1$], [$+1$], [Yukawa matrices, $Tg^a$ on the $Dbrk_1$ line],
  ),
  caption: [Graphs of @fig:g1_cFFS-diagrams. $chi_(f_1)$, $chi_(f_2)$: the line is the sterile component $chi_f = P_(overline(c)_f) Psi_f$ at $Dbrk_1$; $phi$, $phi'$: scalar lines. Sums over $b$, $d$ (all twelve values) and over the scalar index $j$ are implied; (a), (b) exist only for $a$ in $"SU"(3)_c$ or $"SU"(2)_L$. No closed loop occurs. (e) exists only for the four quark vertices of @eq:g1_cFFS-legs, with the Yukawa coupling of the other quark flavour on the two internal vertices ($f_3$, $f_4$ are fixed by the chirality flips at the three Yukawa vertices); for $(L, E, H)$ and $(E, L, H^dagger)$ there is no such graph. Its group factor is not needed, since the graph vanishes at $Top(0)$ for every group structure (@eq:g1_cFFS-extghost).],
) <tab:g1_cFFS-diagrams>

The signs follow from Wick's theorem. In (c)--(e) the external $gh^a$ stands first in $Dbrk_1$ and the rest are even bilinears of an open line: sign $+1$. In (a1) the product of the four vertices, in the field order of @eq:Delta1, @eq:Lgf, @eq:LYuk, @eq:Lferm, is $[gh^b overline(Psi)_2 Psi^"ext"_2] [(partial agh^c) gh^a] [overline(Psi)_1 Psi_2] [overline(Psi)^"ext"_1 Psi_1]$; bringing it to $gh^a overline(Psi)^"ext"_1 (Psi_1 overline(Psi)_1) (Psi_2 overline(Psi)_2) Psi^"ext"_2 (gh^b agh^c)$ is an even permutation, and $(gh^b agh^c)$ gives the propagator $i delta^(b c) \/ ell^2$ of @eq:propagators. The same count gives $+1$ for (a2), (b1), (b2).

*Integrands.* _Scalar line._ For both $phi = H$ and $phi = H^dagger$, the scalar propagator joined to the vertex @eq:V-ssv with the external $phi_i$ gives, at the Yukawa vertex $-i bold(Y)^j P_y$ with the internal field $phi_j$,
$ -i thin bold(Y)^j (Tg^d_phi)_(j i) P_y thin ((ell - p_phi)^mu)/(ell^2 - m^2) , $ <eq:g1_cFFS-scalar>
where $ell$ is the momentum flowing from the Yukawa vertex to the $H^dagger H A$ vertex and $mu$, $d$ are the indices of the gauge leg. _Momentum routing._ $k$ is the momentum of the internal fermion at $Dbrk_1$, along the fermion arrow. In (a1), (b1) the ghost brings $k - p$ into $Dbrk_1$ (outgoing momentum of the antighost leg of @eq:V-ghost), in (a2), (b2) it brings $p' - k$. The gauge line carries $k - p - q$ in (a1), (b1) and $k - p' + q$ in (a2), (b2):
$ G_"(a1)" = intk & (i gbar^mu P_(c_1) Tg^d_(f_1)) (i (sl(k) + sl(p)_phi))/((k + p_phi)^2) (-i bold(Y)^i P_y) (i sl(k))/(k^2) Tg^b_(f_2) (slh(p) P_(overline(c)_(f_2)) - slh(k) P_y) \
  & times (-i g_(mu nu))/((k - p - q)^2) (- fg^(b d a) (k - p)^nu) i/((k - p)^2) , $ <eq:g1_cFFS-int-a1>
$ G_"(a2)" = intk & Tg^b_(f_1) (slh(k) P_y - slh(p)' P_(c_1)) (i sl(k))/(k^2) (-i bold(Y)^i P_y) (i (sl(k) - sl(p)_phi))/((k - p_phi)^2) (i gbar^mu P_y Tg^d_(f_2)) \
  & times (-i g_(mu nu))/((k - p' + q)^2) (fg^(b d a) (k - p')^nu) i/((k - p')^2) , $ <eq:g1_cFFS-int-a2>
$ G_"(b1)" = intk & (-i bold(Y)^j (Tg^d_phi)_(j i) P_y) ((k - p' - p_phi)^mu)/((k - p')^2 - m^2) (i sl(k))/(k^2) Tg^b_(f_2) (slh(p) P_(overline(c)_(f_2)) - slh(k) P_y) \
  & times (-i g_(mu nu))/((k - p - q)^2) (- fg^(b d a) (k - p)^nu) i/((k - p)^2) , $ <eq:g1_cFFS-int-b1>
$ G_"(b2)" = intk & Tg^b_(f_1) (slh(k) P_y - slh(p)' P_(c_1)) (i sl(k))/(k^2) (-i bold(Y)^j (Tg^d_phi)_(j i) P_y) ((p - k - p_phi)^mu)/((k - p)^2 - m^2) \
  & times (-i g_(mu nu))/((k - p' + q)^2) (fg^(b d a) (k - p')^nu) i/((k - p')^2) . $ <eq:g1_cFFS-int-b2>
The products of the numerical factors are $+1$, $-1$, $-1$, $+1$ in (a1), (a2), (b1), (b2). In (c)--(e) the $Dbrk_1$ vertex stands between two internal propagators; with the fermion entering with $k$ and leaving with $k + q$, @eq:m-sterile gives, for the neighbouring vertices with physical projections,
$ P_(c_f) (i (sl(k) + sl(q)))/((k + q)^2) thin Tg^a_f (slh(k) P_(overline(c)_f) - (slh(k) + slh(q)) P_(c_f)) (i sl(k))/(k^2) P_(overline(c)_f)
  = i thin Tg^a_f [ (i (slb(k) + slb(q)))/((k + q)^2) (hat(k)^2)/(k^2) - ((hat(k) + hat(q))^2)/((k + q)^2) (i slb(k))/(k^2) ] P_(overline(c)_f) . $ <eq:g1_cFFS-extghost>
The bracket vanishes identically at $q = 0$; every other factor of these graphs (vertices, propagators, the scalar line) is independent of $q$ at $Top(0)$.

*$Rop$-operation.* At one loop $Rbar G = G$, and $Top(0)$ is the integrand at zero external momenta and $m = 0$ with $D$ on every line. _External ghost._ By @eq:g1_cFFS-extghost, $Top(0) G_"(c)" = Top(0) G_"(d)" = Top(0) G_"(e)" = 0$ identically: $N_0 = N_1 = 0$ for all six graphs. _Internal ghost._ The terms $slh(p)$, $slh(p)'$ of @eq:V-delta (external sterile leg) carry an external momentum and are removed by $Top(0)$. The remaining Dirac strings, with $P_y = P_(overline(c)_1)$, $P_(c_1) sl(k) P_y = slb(k) P_y$ and $P_y sl(k) P_y = slh(k) P_y$, are
$ "(a1)": slb(k) P_(c_1) sl(k) P_y sl(k) slh(k) P_y = bar(k)^2 hat(k)^2 P_y, quad
  "(a2)": slh(k) P_y sl(k) P_y sl(k) slb(k) P_y = hat(k)^2 bar(k)^2 P_y, \
  "(b1)": P_y sl(k) slh(k) P_y = hat(k)^2 P_y, quad
  "(b2)": slh(k) P_y sl(k) P_y = hat(k)^2 P_y , $
where in (a) $gbar^mu$ was contracted with the momentum $k^nu$ of the ghost vertex through $g_(mu nu)$, and in (b) the scalar factor $plus.minus k^mu$ of @eq:g1_cFFS-scalar with $k^nu$ gives $plus.minus k^2$. Hence
$ Top(0) G_"(a1)" = - fg^(b d a) Tg^d_(f_1) bold(Y)^i Tg^b_(f_2) P_y intk bar(k)^2 hat(k)^2 D^4, quad
  Top(0) G_"(a2)" = - fg^(b d a) Tg^b_(f_1) bold(Y)^i Tg^d_(f_2) P_y intk bar(k)^2 hat(k)^2 D^4, \
  Top(0) G_"(b1)" = fg^(b d a) bold(Y)^j (Tg^d_phi)_(j i) Tg^b_(f_2) P_y intk k^2 hat(k)^2 D^4, quad
  Top(0) G_"(b2)" = - fg^(b d a) Tg^b_(f_1) bold(Y)^j (Tg^d_phi)_(j i) P_y intk k^2 hat(k)^2 D^4 . $ <eq:g1_cFFS-T0>
With @eq:tensor-red, @eq:tensor-red-split and @eq:logpole (in the four-index reduction of $bar(k)^2 hat(k)^2$ only the pairing $gmb_(alpha beta) gmh_(gamma delta) g^(alpha beta) g^(gamma delta) = 4 (-2 eta)$ survives),
$ intk bar(k)^2 hat(k)^2 D^4 -> ilf (- (8 eta)/(d (d+2))) 1/ep -> ilf (- eta/3) 1/ep, \
  intk k^2 hat(k)^2 D^4 -> ilf (- (2 eta)/d) 1/ep -> ilf (- eta/2) 1/ep . $ <eq:g1_cFFS-ints>
Every graph is proportional to $eta$, so $N_0 = 0$, and
$ N_1^"(a1)" = 1/3 fg^(b d a) Tg^d_(f_1) bold(Y)^i Tg^b_(f_2) P_y, quad
  N_1^"(a2)" = 1/3 fg^(b d a) Tg^b_(f_1) bold(Y)^i Tg^d_(f_2) P_y, \
  N_1^"(b1)" = -1/2 fg^(b d a) bold(Y)^j (Tg^d_phi)_(j i) Tg^b_(f_2) P_y, quad
  N_1^"(b2)" = 1/2 fg^(b d a) Tg^b_(f_1) bold(Y)^j (Tg^d_phi)_(j i) P_y . $ <eq:g1_cFFS-N>
Exchanging $b <-> d$ in (a2) and using $fg^(d b a) = -fg^(b d a)$: $N_1^"(a1)" + N_1^"(a2)" = 0$. In (b) @eq:g1_cFFS-rigid applied to $bold(Y)^j$ and $[Tg^b_phi, Tg^d_phi] = i fg^(b d x) Tg^x_phi$ (valid for $phi = H^dagger$ as well) give
$ N_1^"(b1)" + N_1^"(b2)" = 1/2 fg^(b d a) bold(Y)^k (Tg^b_phi Tg^d_phi)_(k i) P_y = i/4 thin fg^(b d a) fg^(b d x) bold(Y)^k (Tg^x_phi)_(k i) P_y = i/4 thin C_"ad"^((a)) bold(Y)^j (Tg^a_phi)_(j i) P_y . $ <eq:g1_cFFS-Nb>
For a gluon ghost $Tg^A_phi = 0$, (b1) and (b2) vanish separately, and $N_1^"(a1)" = -N_1^"(a2)" = -i/6 thin C_"ad"^((G)) Tg^A_(f_1) bold(Y)^i P_y$ is non-zero graph by graph. For a $W$ ghost one of the two generators in (a1), (a2) acts on an $"SU"(2)_L$ singlet, so both vanish, and of (b1), (b2) only the graph with $Dbrk_1$ on the doublet leg survives. $N_1$ is four-dimensional, $Pphys N_1 = N_1$.

*Result.*
#result[
  $ (Dbrk_1 dot Gamma)^((1))_(gh^a overline(Psi)_(f_1) Psi_(f_2) phi_i) |_"pole" = lf 1/ep N_0, quad quad N_0 = 0 , $ <eq:g1_cFFS-pole>
  $ X_(1; gh^a overline(Psi)_(f_1) Psi_(f_2) phi_i) & = lf thin Pphys N_1, \
    Pphys N_1 & = i/4 thin C_"ad"^((a)) bold(Y)^j (Tg^a_phi)_(j i) P_y = i/4 thin C_"ad"^((a)) (Tg^a_(f_1) bold(Y)^i - bold(Y)^i Tg^a_(f_2)) P_y . $ <eq:g1_cFFS-X>
  Only $a = W^I$ contributes ($C_"ad"^((B)) = 0$, $Tg^A_phi = 0$). As a local functional (one term per vertex of @eq:g1_cFFS-legs, ghost leftmost),
  $ X_1 & supset lf sum_("vertices") i/4 thin C_"ad"^((W)) integral dif^4 x thin gh^a_((W)) thin overline(psi)_(f_1) (Tg^a_(f_1) bold(Y)^i - bold(Y)^i Tg^a_(f_2)) P_y psi_(f_2) phi_i \
    & = lf gw^2/2 sum_(k = u, d, e) ( - N^(O_k W)_(Y L) + N^(O_k^dagger W)_(Y R) ) , $ <eq:g1_cFFS-X-op>
  with $C_"ad"^((W)) = 2 gw^2$. Equivalently, with @eq:btable-yuk-GW and @eq:btable-higgs,
  $ X_1 |_(gh overline(psi) psi H) = lf gw^2/2 ( [bop Q^(Q W)_R]_"Yuk" + [bop Q^(L W)_R]_"Yuk" ) = - lf gw^2/2 thin [bop Q^((W))_sigma]_"Yuk" . $ <eq:g1_cFFS-X-b>
]
Translation: $N^(O W)_(Y L) = i y_k integral O[overline(psi) -> overline(psi) wm^((W))] = - i integral gh^a overline(psi)_(f_1) Tg^a_(f_1) bold(Y)^i P_y psi_(f_2) phi_i$, since $wm = Tg^a gh^a$ and $gh^a$ anticommutes with $overline(psi)$; its vertex function (ghost first) is $-i Tg^a_(f_1) bold(Y)^i P_y$. Likewise $N^(O^dagger W)_(Y R)$ has the vertex function $-i bold(Y)^i Tg^a_(f_2) P_y$. Here $bold(Y)^i$ contains $y_k$, which in @tab:Nbasis is the explicit factor $y_k$. For $O_k$ the doublet is $f_1$ ($Tg^I_(f_2) = 0$), for $O_k^dagger$ it is $f_2$ ($Tg^I_(f_1) = 0$).

#figure(
  table(
    columns: 7,
    align: (center, center, left, left, left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([vertex], [$r$], [$N_1^"(a1)"$], [$N_1^"(a2)"$], [$N_1^"(b1)"$], [$N_1^"(b2)"$], [coefficient in $16 pi^2 X_1$]),
    [$O_u$, $O_d$], [$G$], [$-i/2 gs^2 Tg^A_(f_1) bold(Y)^i$], [$+i/2 gs^2 Tg^A_(f_1) bold(Y)^i$], [$0$], [$0$], [$0 dot N^(O G)_(Y L)$],
    [$O_u^dagger$, $O_d^dagger$], [$G$], [$-i/2 gs^2 Tg^A_(f_1) bold(Y)^i$], [$+i/2 gs^2 Tg^A_(f_1) bold(Y)^i$], [$0$], [$0$], [$0 dot N^(O^dagger G)_(Y L)$],
    [$O_u$, $O_d$, $O_e$], [$W$], [$0$], [$0$], [$0$], [$i/2 gw^2 Tg^I_(f_1) bold(Y)^i$], [$-gw^2\/2 dot N^(O_k W)_(Y L)$],
    [$O_u^dagger$, $O_d^dagger$, $O_e^dagger$], [$W$], [$0$], [$0$], [$-i/2 gw^2 bold(Y)^i Tg^I_(f_2)$], [$0$], [$+gw^2\/2 dot N^(O_k^dagger W)_(Y R)$],
    [all six], [$B$], [$0$], [$0$], [$0$], [$0$], [$0 dot N^(O)_(Y B)$],
  ),
  caption: [Standard-Model values of @eq:g1_cFFS-N (the projector $P_y$ on the right is omitted) and the coefficients of the basis monomials of @tab:Nbasis in $16 pi^2 X_1$. $C_"ad"^((G)) = 3 gs^2$, $C_"ad"^((W)) = 2 gw^2$, $C_"ad"^((B)) = 0$. The external-ghost graphs (c)--(e) give $0$ everywhere. In the $G$ rows $Tg^A_(f_1) bold(Y)^i = bold(Y)^i Tg^A_(f_2)$ ($H$ is a colour singlet).],
) <tab:g1_cFFS-sm>

#check[*Cancellation of the pole.* By @sec:m-oneloop, $sum_G Kop Top(0) G + i V_(bop_d delta S^((1)))[gh overline(Psi) Psi phi] = 0$, with $V$ the vertex function of @sec:dictionary. Here $N_0 = 0$, so the $gh overline(Psi) Psi H$ vertex of $bop_d delta S^((1))$ must vanish. Of the monomials of $delta S^((1))$ only the source terms $Q^(f r)_R$, $Q^((r))_sigma$ have such a vertex ($Q^+_(y_k)$ is closed, and $Q^f_K$, $Q^(f r)_V$ and the Higgs and gauge monomials have none in this class); it comes from part (iii) of @eq:dict: the coefficient of the source times the variation of @eq:LYuk under the BRST transformation of that field alone. By @sec:g0_ghost_sources all fermion and Higgs source terms of factor $r$ have the same pole coefficient, $-lf (1\/ep) C_"ad"^((r)) \/ 2$. The variations of the three fields of a Yukawa monomial then add up to $sop$ of the monomial, which vanishes in $d$ dimensions ($O_k$ contains only physical components and no $gamma$ matrix). Hence $V_(bop_d delta S^((1))) = 0$ in this class, in agreement with @eq:g1_cFFS-pole.]

#check[*Abelian ghost and normalisation.* $X_1$ has no component along $gh^B$, in particular none along $N^(O_k)_(Y B) - N^(O_k^dagger)_(Y B)$: at one loop $gh^B$ can only sit on $Dbrk_1$ with both fermion legs internal (topologies (c)--(e)), and these vanish at $Top(0)$. This confirms the condition required in @sec:normalisation for (N4) to admit a solution. The gluon ghost gives no $X_1$ either: the two graphs (a1), (a2) cancel, and the scalar does not couple to gluons.]

#check[*Independence and reality.* $N_1$ is independent of $MUV$ and of $m^2$ ($m^2$ has degree two). $g5$ enters only through $P_y$; a vertex and its conjugate have the same form. The combination $-N^(O_k W)_(Y L) + N^(O_k^dagger W)_(Y R)$ is the one produced by $bop$ of the hermitian monomials $Q^(f W)_R$, $Q^((W))_sigma$ (@eq:g1_cFFS-X-b); the opposite relative sign between $O_k$ and $O_k^dagger$, which an error in the relative sign of (b1) and (b2) would give, lies outside $bop("span"{Q_alpha})$ in this class. The Dirac strings, the group factors of @tab:g1_cFFS-diagrams and @tab:g1_cFFS-sm for all six vertices and all twelve ghosts, the identity @eq:g1_cFFS-extghost and the Wick signs were checked numerically with an explicit $8 times 8$ representation of the Clifford algebra and explicit $"SU"(3)$, $"SU"(2)$ matrices (@sec:supporting).]

#check[*Consistency conditions for @sec:st-solution.* In this class only $Q^(f r)_R$ and $Q^((r))_sigma$ contribute to $bop Delta S^((1))$ (@eq:btable-fermion with its Yukawa parts @eq:btable-yuk-B and @eq:btable-yuk-GW, and @eq:btable-higgs). With $Delta S^((1))_"ST" = lf sum_alpha r_alpha Q_alpha$ as in @sec:st-solution, write $r[Q_alpha]$ for the coefficient $r_alpha$ of the monomial $Q_alpha$ (it does not contain the factor $lf$). The condition $bop Delta S^((1))_"ST" = -X_1$ restricted to $gh overline(psi) psi H$ reads
$ & r[Q^(Q G)_R] = r[Q^(U G)_R] = r[Q^(D G)_R], \
  & r[Q^(Q W)_R] - r[Q^((W))_sigma] = r[Q^(L W)_R] - r[Q^((W))_sigma] = - gw^2/2 = - C_"ad"^((W))/4, \
  & -1/6 r[Q^(Q B)_R] + 2/3 r[Q^(U B)_R] - 1/2 r[Q^((B))_sigma] = 0 quad (k = u), \
  & -1/6 r[Q^(Q B)_R] - 1/3 r[Q^(D B)_R] + 1/2 r[Q^((B))_sigma] = 0 quad (k = d), \
  & 1/2 r[Q^(L B)_R] - r[Q^(E B)_R] + 1/2 r[Q^((B))_sigma] = 0 quad (k = e) . $
The last three are $sum_(phi') q_(phi') (O_k) thin r_(phi') = 0$, summed over the three fields $phi'$ of $O_k$, with the hypercharges $q_(phi') (O_k)$ of @eq:btable-yuk-B and $r_(phi')$ the coefficient $r[Q^(f B)_R]$ of the multiplet $f$ of $phi'$, or $r[Q^((B))_sigma]$ for $phi' = H, tilde(H)$. These conditions tie the source counterterms of different multiplets to each other and to the Higgs source; with (N4) all $B$ coefficients vanish and the last three hold trivially. In @sec:st-solution, $r[Q^(f r)_R] = -C_"ad"^((r))\/4$ is fixed by the class $gh gh overline(srcF) psi$ (@eq:sts-r-R) and $r[Q^((W))_sigma] = 0$ by the class $gh gh srcH^dagger H$ (@eq:sts-r-sigma); the first two conditions then hold. In @sec:g1_cFF the analogous internal-ghost graphs (c1), (c2) cancel at $q = 0$ for the same reason as (a1), (a2) here.]

#general[The vanishing of (c)--(e) rests on momentum conservation at $Dbrk_1$: at zero ghost momentum its two fermion lines carry the same momentum, and @eq:g1_cFFS-extghost vanishes for any loop momentum. This holds in every graph at every order in which the ghost of this class sits on a $Dbrk_1$ vertex (or on the part $z_R Dbrk_1$ of $bop_d delta S^((L))$, @eq:m-bYR) with two internal fermion lines, including the counterterm-inserted forests of @eq:m-Rbar. Since the graphs of the class have $dg(G) = 0$ and are evaluated at zero external momenta, $X_l$ in this class comes at every order only from graphs in which the external ghost does not sit on such a vertex: internal-ghost graphs and the insertions of the momentum-independent $gh overline(psi) psi H$ vertices of $bop_d Delta S^((L))$ ($[bop Q_R]_"Yuk"$, $[bop Q_sigma]_"Yuk"$).]
