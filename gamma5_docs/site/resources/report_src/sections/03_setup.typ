#import "../macros.typ": *
#import "../diagrams.typ": *

= Lagrangian, notation and Feynman rules <sec:setup>

This section fixes the conventions of the calculation: the Dirac algebra, the fields, the action, the BRST transformations, the Feynman rules and the loop integrals. Symbols specific to the method (@sec:method), to the operator bases (@sec:bases) or to a single class of Green's functions are defined where they first appear.

*Conventions for symbols.* @tab:alphabet states once which letters are reserved for which objects. Repeated Lorentz indices and repeated adjoint indices are summed, unless stated otherwise.

#[
#show figure: set block(breakable: true)
#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([symbol], [meaning]),
    [$mu, nu, rho, sigma, alpha, beta$], [Lorentz indices, $d$-dimensional unless marked by a bar (four-dimensional part) or a hat ($(d-4)$-dimensional part). In a formula or paragraph with a BRST source only $mu, nu, alpha, beta$ are Lorentz indices ($rho$, $sigma$ are then sources).],
    [$a, b, c, d, e, x, y, z$], [adjoint indices running over the twelve gauge fields, $a in {A, I, B}$. The letter $d$ is an adjoint index only as a superscript or subscript of a field or structure constant ($A^d_mu$, $fg^(a c d)$); otherwise $d$ is the dimension, except in the $"SU"(3)$ tensor $d^(A B C)$ (@sec:bases), the down-quark labels $d_upright(R)$, $y_d$, $O_d$ and the names of monomials ($Q_(3 d)$, $N_(d f)$).],
    [$A, B, C, D, E$ #h(0.3em) / #h(0.3em) $I, J, K$], [adjoint indices of $"SU"(3)_c$ ($1, dots, 8$) and of $"SU"(2)_L$ ($1, 2, 3$). The value $a = B$ is the hypercharge field; in a formula with $B_mu$ or $gh^B$ the $"SU"(3)$ indices are $A, C, D, E$ only.],
    [$alpha, beta$ #h(0.3em) / #h(0.3em) $i, j, k, l$], [colour indices of triplets ($1, 2, 3$; only in @tab:fields and in the Yukawa terms, which have no Lorentz index) and $"SU"(2)_L$ doublet indices ($1, 2$). A subscript or superscript $i$ is a doublet index; a factor $i$ is the imaginary unit.],
    [$f in {Q, U, D, L, E}$], [fermion multiplet (@tab:fields); $f = H$ denotes the Higgs doublet where group invariants are listed. $N^c_f$ is the colour multiplicity: $N^c_f = 3$ for $f = Q, U, D$ and $N^c_f = 1$ for $f = L, E$.],
    [$r in {G, W, B}$], [gauge factor $"SU"(3)_c$, $"SU"(2)_L$, $"U"(1)_Y$; $a in r$ means that the adjoint index $a$ belongs to the factor $r$.],
    [$c_f$; #h(0.3em) $"lh"$, $"rh"$], [chirality of the physical field of $f$ ($upright(L)$ or $upright(R)$, projectors $PL$, $PR$); the subscripts $"lh"$ and $"rh"$ mark sums over the left-handed multiplets ($Q$, $L$) and the right-handed multiplets ($U$, $D$, $E$), @eq:sumLR.],
    [$gs$, $gw$, $gy$; #h(0.3em) $g^(mu nu)$], [couplings of $"SU"(3)_c$, $"SU"(2)_L$, $"U"(1)_Y$. The metric is written with two indices only ($g^(mu nu)$, $gmb^(mu nu)$, $gmh^(mu nu)$), so $g$ without indices is always the $"SU"(2)_L$ coupling.],
    [$ep$; #h(0.3em) $epsilon$], [the regulator, $d = 4 - 2 ep$ (lunate glyph). The glyph $epsilon$ with indices is a Levi-Civita symbol: $epsilon^(mu nu rho sigma)$ (four-dimensional), $epsilon^(I J K)$ ($"SU"(2)$), $epsilon_(i j)$ (doublet).],
    [$eta$], [evanescent tag: $gmh^mu_(space mu) = -2 eta$, set to $eta = ep$ at the end (@eq:eta).],
    [$gh^a$, $agh^a$, $b^a$], [ghost, antighost and auxiliary field. The letter $omega$ is never used for anything else.],
    [$srcA^mu_a$, $srcC_a$, $srcF_f$, $srcH$; #h(0.3em) $Phi^*$], [BRST sources of $A^a_mu$, $gh^a$, $overline(Psi)_f$ and $H^dagger$ (@eq:Lext); $Phi^*$ is the source paired with the field $Phi$ (@eq:Phistar).],
    [$dg(G)$, $Top(n)$], [superficial degree of divergence of a graph $G$ and the expansion operator to degree $n$ (@sec:integrals).],
    [$lf$], [one-loop factor. Every renormalisation constant ($delta Z_x$, $delta m^2$, $delta lambda$, $delta y_k$) includes the factor $lf 1/ep$.],
  ),
  caption: [Reserved symbols.],
) <tab:alphabet>
]

== Dimensions and the BMHV algebra <sec:bmhv>

Space-time has $d = 4 - 2 ep$ dimensions and metric $g^(mu nu) = "diag"(+,-,-,-,-,dots)$. Every $d$-dimensional index is split into a four-dimensional and a $(d-4)$-dimensional ("evanescent") part, written with a bar and a hat:
$ g^(mu nu) = gmb^(mu nu) + gmh^(mu nu), quad
  k^mu = bar(k)^mu + hat(k)^mu, quad
  gamma^mu = gbar^mu + ghat^mu, quad
  sl(k) = slb(k) + slh(k) . $ <eq:split>
Here $bar(k)^mu = gmb^(mu nu) k_nu$, $hat(k)^mu = gmh^(mu nu) k_nu$, $slb(k) = gbar^mu k_mu$ and $slh(k) = ghat^mu k_mu$. The complete set of algebraic rules is the following.

*Metric.*
$ gmb^(mu nu) gmb_(nu rho) = gmb^mu_(space rho), quad
  gmh^(mu nu) gmh_(nu rho) = gmh^mu_(space rho), quad
  gmb^(mu nu) gmh_(nu rho) = 0, quad
  g^mu_(space mu) = d, quad gmb^mu_(space mu) = 4, quad gmh^mu_(space mu) = d - 4 = -2 ep . $ <eq:metric>
Consequently $k^2 = bar(k)^2 + hat(k)^2$ and $k dot p = bar(k) dot bar(p) + hat(k) dot hat(p)$.

*Clifford algebra.*
$ {gamma^mu, gamma^nu} = 2 g^(mu nu), quad
  {gbar^mu, gbar^nu} = 2 gmb^(mu nu), quad
  {ghat^mu, ghat^nu} = 2 gmh^(mu nu), quad
  {gbar^mu, ghat^nu} = 0 . $ <eq:clifford>

*The matrix $g5$ and the Levi-Civita tensor.* The Levi-Civita tensor $epsilon^(mu nu rho sigma)$ is four-dimensional (it vanishes if any index is evanescent), with $epsilon^(0123) = +1$, and
$ g5 = -i/(4!) thin epsilon_(mu nu rho sigma) gamma^mu gamma^nu gamma^rho gamma^sigma = i gamma^0 gamma^1 gamma^2 gamma^3, quad g5^2 = 1, quad g5^dagger = g5 . $ <eq:g5def>
It follows that
$ {gbar^mu, g5} = 0, quad [ghat^mu, g5] = 0, quad "hence" quad {gamma^mu, g5} = 2 ghat^mu g5 . $ <eq:g5comm>
The chiral projectors are $PL = (1 - g5)/2$ and $PR = (1 + g5)/2$, so that
$ PL^2 = PL, quad PR^2 = PR, quad PL PR = 0, quad
  gbar^mu P_(upright(L),upright(R)) = P_(upright(R),upright(L)) gbar^mu, quad
  ghat^mu P_(upright(L),upright(R)) = P_(upright(L),upright(R)) ghat^mu . $ <eq:proj>
Two consequences are used constantly:
$ PR gamma^mu PL = gbar^mu PL, quad PL gamma^mu PL = ghat^mu PL, quad "and the same with" upright(L) <-> upright(R) . $ <eq:sandwich>

*Contractions.*
$ gamma^mu gamma_mu = d, quad gbar^mu gbar_mu = 4, quad ghat^mu ghat_mu = d - 4, quad
  gbar^mu ghat_mu = 0, $ <eq:contr1>
$ gamma^mu gamma^nu gamma_mu = -(d-2) gamma^nu, quad
  gbar^mu gbar^nu gbar_mu = -2 gbar^nu, quad
  gbar^mu ghat^nu gbar_mu = -4 ghat^nu, \
  ghat^mu gbar^nu ghat_mu = -(d-4) gbar^nu, quad
  ghat^mu ghat^nu ghat_mu = -(d-6) ghat^nu . $ <eq:contr2>
All longer contractions follow from @eq:clifford by anticommuting the contracted pair together.

*Traces.* The trace is cyclic and $tr bb(1) = 4$. Traces of an odd number of $gamma$ matrices vanish, with or without $g5$. Without $g5$ the trace of a product is given by the usual recursion with the metric of each factor ($g$, $gmb$ or $gmh$ according to the two matrices paired, and $0$ for a barred with a hatted one). With one $g5$:
$ tr(g5) = tr(gamma^mu gamma^nu g5) = 0, quad
  tr(gamma^mu gamma^nu gamma^rho gamma^sigma g5) = -4 i thin epsilon^(mu nu rho sigma) . $ <eq:tr5>
Because $epsilon$ is four-dimensional only the barred parts of the four matrices contribute in @eq:tr5. Traces with $g5$ and more matrices are reduced with @eq:g5def and @eq:clifford; a hatted matrix can be moved through $g5$ freely and through barred matrices with a sign, so every trace factorises into a four-dimensional trace times a purely evanescent trace (no $g5$), divided by $tr bb(1) = 4$: $tr(overline(Gamma) hat(Gamma)) = tr(overline(Gamma)) tr(hat(Gamma)) \/ 4$ for $overline(Gamma)$ a product of barred matrices and $g5$ and $hat(Gamma)$ a product of hatted matrices.

*Contractions of two Levi-Civita tensors* are four-dimensional:
$ epsilon^(mu_1 mu_2 mu_3 mu_4) epsilon_(nu_1 nu_2 nu_3 nu_4) = - det mat(delim: "[", gmb^(mu_1)_(space nu_1), gmb^(mu_1)_(space nu_2), gmb^(mu_1)_(space nu_3), gmb^(mu_1)_(space nu_4); dots.v, , , dots.v; gmb^(mu_4)_(space nu_1), gmb^(mu_4)_(space nu_2), gmb^(mu_4)_(space nu_3), gmb^(mu_4)_(space nu_4)) , $ <eq:epseps>
i.e. minus the determinant of the $4 times 4$ matrix with entries $gmb^(mu_i)_(space nu_j)$.

*Hermiticity.* $gamma^0 (gamma^mu)^dagger gamma^0 = gamma^mu$ for barred and hatted matrices alike, and $overline(Psi) = Psi^dagger gamma^0$.

*Physical projection.* For an expression $O$ in which all index contractions have been carried out,
$ Pphys O = O |_(ghat = 0, thin gmh = 0), quad Evan = 1 - Pphys . $ <eq:Pproj>
$Pphys$ does not commute with contractions: $gmh^mu_(space mu) = d - 4$ must be evaluated before $Pphys$ is applied. An operator with $Pphys O = 0$ is called evanescent.

== Field content and multiplet notation <sec:fields>

The gauge group is $"SU"(3)_c times "SU"(2)_L times "U"(1)_Y$ with couplings $gs$, $gw$, $gy$ and gauge fields $G^A_mu$ ($A = 1,dots,8$), $W^I_mu$ ($I = 1,2,3$) and $B_mu$. We collect them in one field $A^a_mu$ whose index $a$ runs over the twelve values ${A, I, B}$. The matter content is one family of fermions and one Higgs doublet,

#figure(
  table(
    columns: 7,
    align: center,
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([multiplet], [physical field], [chirality $c_f$], [$"SU"(3)_c$], [$"SU"(2)_L$], [$Y$], [sterile partner]),
    [$Q$], [$q_upright(L)^(alpha i)$], [L], [$bold(3)$], [$bold(2)$], [$1/6$], [$chi_(Q)$, right-handed],
    [$U$], [$u_upright(R)^alpha$], [R], [$bold(3)$], [$bold(1)$], [$2/3$], [$chi_(U)$, left-handed],
    [$D$], [$d_upright(R)^alpha$], [R], [$bold(3)$], [$bold(1)$], [$-1/3$], [$chi_(D)$, left-handed],
    [$L$], [$l_upright(L)^i$], [L], [$bold(1)$], [$bold(2)$], [$-1/2$], [$chi_(L)$, right-handed],
    [$E$], [$e_upright(R)$], [R], [$bold(1)$], [$bold(1)$], [$-1$], [$chi_(E)$, left-handed],
    [$H$], [$H_i$], [scalar], [$bold(1)$], [$bold(2)$], [$1/2$], [none],
  ),
  caption: [Matter multiplets. $alpha = 1,2,3$ is a colour index and $i = 1,2$ an $"SU"(2)_L$ doublet index. The electric charge is $Q_"em" = tau^3 \/ 2 + Y$ on doublets and $Q_"em" = Y$ on singlets, with $tau^3$ the third Pauli matrix.],
) <tab:fields>

Each physical chiral fermion $psi_f$, $f in {Q,U,D,L,E}$, has a sterile partner $chi_f$ of the opposite chirality. The partner carries the same colour and $"SU"(2)_L$ indices as $psi_f$, transforms like $psi_f$ under _global_ transformations, and is a singlet under all three _gauge_ groups: it never couples to a gauge field. The pair is packaged in one Dirac spinor
$ Psi_f = psi_f + chi_f, quad psi_f = P_(c_f) Psi_f, quad chi_f = P_(overline(c)_f) Psi_f, quad
  overline(psi)_f = overline(Psi)_f P_(overline(c)_f), quad overline(chi)_f = overline(Psi)_f P_(c_f) , $ <eq:dirac-pack>
where $P_(c_f) = PL$ for $f = Q, L$ and $P_(c_f) = PR$ for $f = U, D, E$, and $P_(overline(c)_f)$ is the other projector. We write $Q, U, D, L, E$ for the five Dirac spinors $Psi_f$.

*Coupling-weighted generators.* On each multiplet the twelve generators are multiplied by their coupling,
$ Tg^A_f = gs thin T^A, quad Tg^I_f = gw thin tau^I / 2, quad Tg^B_f = gy thin Y_f , $ <eq:Tg>
where $T^A = lambda^A\/2$ acts on colour triplets and is zero on colour singlets, $tau^I$ are the Pauli matrices acting on doublets and zero on $"SU"(2)_L$ singlets, and $Y_f$ is the hypercharge of @tab:fields. Explicitly, listing the components $a = (A, I, B)$,
$ Tg^a_Q = (gs T^A, thin gw tau^I/2, thin gy/6), quad
  Tg^a_U = (gs T^A, thin 0, thin 2 gy/3), quad
  Tg^a_D = (gs T^A, thin 0, thin -gy/3), \
  Tg^a_L = (0, thin gw tau^I/2, thin -gy/2), quad
  Tg^a_E = (0, thin 0, thin -gy), quad
  Tg^a_H = (0, thin gw tau^I/2, thin gy/2) . $ <eq:Tg-explicit>
They obey
$ [Tg^a, Tg^b] = i fg^(a b c) Tg^c, quad
  fg^(A C D) = gs f^(A C D), quad fg^(I J K) = gw epsilon^(I J K), quad fg^(a b c) = 0 "otherwise", $ <eq:fg>
with $f^(A C D)$ the $"SU"(3)$ structure constants and $epsilon^(I J K)$ the Levi-Civita symbol with $epsilon^(1 2 3) = 1$. $fg^(a b c)$ is totally antisymmetric and satisfies the Jacobi identity. The covariant derivative and field strength are
$ D_mu = partial_mu - i Tg^a A^a_mu, quad
  F^a_(mu nu) = partial_mu A^a_nu - partial_nu A^a_mu + fg^(a b c) A^b_mu A^c_nu, quad
  (D_mu gh)^a = partial_mu gh^a + fg^(a b c) A^b_mu gh^c . $ <eq:covD>
On the Higgs doublet $D_mu H = (partial_mu - i gw W^I_mu tau^I\/2 - i (gy\/2) B_mu) H$, and $tilde(H)_i = epsilon_(i j) H^*_j$, with $epsilon_(i j)$ the antisymmetric symbol, $epsilon_(12) = +1$.

*Matrix-valued fields.* A fermion or Higgs line of the multiplet $f$ sees the gauge fields and ghosts only through the matrices
$ Am_mu = Tg^a_f A^a_mu = gs T^A G^A_mu + gw tau^I/2 W^I_mu + gy Y_f B_mu, quad quad wm = Tg^a_f gh^a , $ <eq:matrix-fields>
which act on the colour and $"SU"(2)_L$ indices of $f$; the label $f$ is written as $Am^f_mu$, $wm^f$ only where needed. In this notation @eq:covD and @eq:brst read $D_mu = partial_mu - i Am_mu$ and $sop Am_mu = partial_mu wm - i [Am_mu, wm]$, $sop wm = i wm wm$. Closed fermion loops produce traces $Tr_f$ over the colour and $"SU"(2)_L$ indices of $f$ of products of these matrices, summed over the multiplets; the sums over left- and right-handed multiplets are written
$ sum_"lh" = sum_(f in {Q, L}), quad sum_"rh" = sum_(f in {U, D, E}) . $ <eq:sumLR>

*Group invariants.* The one-loop results depend on the group structure only through the following numbers, defined once here. For a multiplet $f$, a gauge factor $r in {G, W, B}$ and adjoint indices $a, b in r$:
$ Tr(Tg^a_f Tg^b_f) = S_f^((r)) delta^(a b), quad
  fg^(a c e) fg^(b c e) = C_"ad"^((r)) delta^(a b), quad
  sum_(b = 1)^(12) Tg^b_f Tg^b_f = C_f bb(1) . $ <eq:invariants>
In $S_f^((r))$ the trace $Tr$ runs over all colour and $"SU"(2)_L$ components of the multiplet (it includes the colour multiplicity $N^c_f$); in $C_"ad"^((r))$ the sum runs over $c, e$; $C_f$ sums over all twelve generators. $Tr(Tg^a_f Tg^b_f) = 0$ when $a$ and $b$ belong to different factors. For $f = H$ the same definitions give $S_H^((r))$ and $C_H$.

#figure(
  table(
    columns: 7,
    align: center,
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([], [$Q$], [$U$], [$D$], [$L$], [$E$], [$H$]),
    [$S_f^((G))$], [$gs^2$], [$gs^2\/2$], [$gs^2\/2$], [$0$], [$0$], [$0$],
    [$S_f^((W))$], [$3 gw^2\/2$], [$0$], [$0$], [$gw^2\/2$], [$0$], [$gw^2\/2$],
    [$S_f^((B))$], [$gy^2\/6$], [$4 gy^2\/3$], [$gy^2\/3$], [$gy^2\/2$], [$gy^2$], [$gy^2\/2$],
    [$C_f$], [$4/3 gs^2 + 3/4 gw^2 + 1/36 gy^2$], [$4/3 gs^2 + 4/9 gy^2$], [$4/3 gs^2 + 1/9 gy^2$], [$3/4 gw^2 + 1/4 gy^2$], [$gy^2$], [$3/4 gw^2 + 1/4 gy^2$],
  ),
  caption: [Group invariants of @eq:invariants. In addition $C_"ad"^((G)) = 3 gs^2$, $C_"ad"^((W)) = 2 gw^2$, $C_"ad"^((B)) = 0$. The $H$ column gives $S_H^((r))$ and $C_H$.],
) <tab:invariants>

Sums over the left-handed and right-handed fermion multiplets appear in closed fermion loops:
$ S_"lh"^((r)) = S_Q^((r)) + S_L^((r)), quad S_"rh"^((r)) = S_U^((r)) + S_D^((r)) + S_E^((r)), quad
  S^((r)) = S_"lh"^((r)) + S_"rh"^((r)), \
  S_"lh"^((G)) = gs^2, quad S_"rh"^((G)) = gs^2, quad
  S_"lh"^((W)) = 2 gw^2, quad S_"rh"^((W)) = 0, quad
  S_"lh"^((B)) = 2/3 gy^2, quad S_"rh"^((B)) = 8/3 gy^2, \
  S^((G)) = 2 gs^2, quad S^((W)) = 2 gw^2, quad S^((B)) = 10/3 gy^2 . $ <eq:SLR>
The totally symmetric traces that control the gauge anomaly are
$ D^(a b c) = sum_"lh" Tr(Tg^a_f {Tg^b_f, Tg^c_f}) - sum_"rh" Tr(Tg^a_f {Tg^b_f, Tg^c_f}) , $ <eq:Dabc>
and they all vanish for the content of @tab:fields; the explicit sums are given where they are used, in @sec:classes.

== The regularised tree-level action <sec:action>

The tree-level action in $d$ dimensions is $S_0 = integral dif^d x thin cal(L)_0$ with
$ cal(L)_0 = cal(L)_"gauge" + cal(L)_"fermion" + cal(L)_"Higgs" + cal(L)_"Yukawa" + cal(L)_"gf" + cal(L)_"ext" . $ <eq:L0>

*Gauge fields.* All Lorentz indices are $d$-dimensional:
$ cal(L)_"gauge" = -1/4 F^a_(mu nu) F^(a mu nu)
  = -1/4 G^A_(mu nu) G^(A mu nu) - 1/4 W^I_(mu nu) W^(I mu nu) - 1/4 B_(mu nu) B^(mu nu) . $ <eq:Lgauge>

*Fermions.* For each multiplet the physical field has a four-dimensional gauge-covariant kinetic term, the sterile partner a four-dimensional free kinetic term, and the two are coupled by an evanescent kinetic term:
$ cal(L)_"fermion" = sum_f [ & i overline(psi)_f gbar^mu D_mu psi_f + i overline(chi)_f gbar^mu partial_mu chi_f \
   & + i overline(psi)_f ghat^mu partial_mu chi_f + i overline(chi)_f ghat^mu partial_mu psi_f ] . $ <eq:Lferm-weyl>
With @eq:dirac-pack and @eq:sandwich this is identical to
$ cal(L)_"fermion" = sum_f [ overline(Psi)_f thin i sl(partial) thin Psi_f + overline(Psi)_f thin gbar^mu thin Tg^a_f thin P_(c_f) Psi_f thin A^a_mu ] . $ <eq:Lferm>
The kinetic term is the full $d$-dimensional Dirac operator, so the propagator of $Psi_f$ is the ordinary $d$-dimensional one. The gauge interaction involves only $gbar^mu$ and only the physical chirality; since $gbar^mu A_mu = gbar^mu bar(A)_mu$, fermions couple only to the four-dimensional components $bar(A)^a_mu$ of the gauge fields. The last two terms of @eq:Lferm-weyl are not gauge invariant: this is the only place where the regularisation breaks the gauge symmetry.

*Higgs doublet.* With $d$-dimensional indices,
$ cal(L)_"Higgs" = (D_mu H)^dagger (D^mu H) - m^2 H^dagger H - lambda (H^dagger H)^2 . $ <eq:LHiggs>
The unbroken phase is the expansion around $H = 0$; $m^2$ is a real parameter of the Lagrangian, of either sign.

*Yukawa couplings*, with real $y_u, y_d, y_e$:
$ cal(L)_"Yukawa" = & - y_u thin overline(Q) tilde(H) PR U - y_d thin overline(Q) H PR D - y_e thin overline(L) H PR E \
   & - y_u thin overline(U) tilde(H)^dagger PL Q - y_d thin overline(D) H^dagger PL Q - y_e thin overline(E) H^dagger PL L . $ <eq:LYuk>
Colour indices are contracted between the two quark fields and the $"SU"(2)_L$ index of the doublet with that of the Higgs field: $overline(Q) tilde(H) PR U = overline(Q)^(alpha j) epsilon_(j k) H^*_k PR U_alpha$ and $overline(Q) H PR D = overline(Q)^(alpha j) H_j PR D_alpha$. Only physical components enter, because $overline(Q) PR = overline(q)_upright(L)$ and $PR U = u_upright(R)$.

*Gauge fixing and ghosts.* For each of the twelve gauge fields there is a ghost $gh^a$, an antighost $agh^a$ and an auxiliary field $b^a$. In the gauge with parameter $xi$,
$ cal(L)_"gf" = b^a partial^mu A^a_mu + xi/2 b^a b^a + partial^mu agh^a (D_mu gh)^a . $ <eq:Lgf>
Eliminating $b^a = -partial^mu A^a_mu \/ xi$ gives $-1/(2 xi) (partial^mu A^a_mu)^2$. All calculations are done in the Feynman gauge $xi = 1$.

*BRST sources.* Each non-linear BRST variation (@eq:brst below) is coupled to an external source that never propagates:
$ cal(L)_"ext" = & srcA^mu_a (D_mu gh)^a - 1/2 srcC_a fg^(a b c) gh^b gh^c \
   & + sum_f [ overline(srcF)_f thin sop Psi_f - (sop overline(Psi)_f) thin srcF_f ]
   + srcH^dagger thin sop H + (sop H^dagger) thin srcH . $ <eq:Lext>
With the variations @eq:brst-matter inserted and the ghost moved to the left of the matter field, the matter part reads explicitly
$ & sum_f [ i thin overline(srcF)_f thin gh^a Tg^a_f P_(c_f) Psi_f + i thin gh^a thin overline(Psi)_f P_(overline(c)_f) Tg^a_f thin srcF_f ] \
  & + i thin srcH^dagger gh^a Tg^a_H H - i thin gh^a H^dagger Tg^a_H srcH . $ <eq:Lext-explicit>
The source $srcA^mu_a$ is the one called $gamma^mu_a$ in the starter pack; it is renamed to avoid a clash with the Dirac matrices. In every formula that contains $srcA$ or $srcH$, the Lorentz indices are taken from $mu, nu, alpha, beta$ only. The fermion sources $srcF_f$, $overline(srcF)_f$ are distinct from the $Rop$-operation and from $Rbar$, which are always written in calligraphic type. Statistics, ghost number and mass dimension of all fields and sources are listed in @tab:quantum-numbers.

#figure(
  table(
    columns: 11,
    align: center,
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([], [$A^a_mu$], [$Psi_f$], [$H$], [$gh^a$], [$agh^a$], [$b^a$], [$srcA^mu_a$], [$srcC_a$], [$srcF_f$], [$srcH$]),
    [statistics], [B], [F], [B], [F], [F], [B], [F], [B], [B], [F],
    [ghost number], [$0$], [$0$], [$0$], [$+1$], [$-1$], [$0$], [$-1$], [$-2$], [$-1$], [$-1$],
    [dimension], [$1$], [$3/2$], [$1$], [$1$], [$1$], [$2$], [$2$], [$2$], [$3/2$], [$2$],
  ),
  caption: [Statistics (B bosonic, F fermionic, i.e. Grassmann odd), ghost number and mass dimension at $d = 4$. The sources have the statistics opposite to their field.],
) <tab:quantum-numbers>

== BRST transformations and the Slavnov--Taylor operator <sec:brst>

The BRST operator $sop$ is Grassmann odd, raises the ghost number by one, and acts as
$ sop A^a_mu = (D_mu gh)^a = partial_mu gh^a + fg^(a b c) A^b_mu gh^c, quad
  sop gh^a = -1/2 fg^(a b c) gh^b gh^c, quad sop agh^a = b^a, quad sop b^a = 0, $ <eq:brst>
$ sop Psi_f = i gh^a Tg^a_f P_(c_f) Psi_f, quad
  sop overline(Psi)_f = i overline(Psi)_f P_(overline(c)_f) Tg^a_f gh^a, quad
  sop H = i gh^a Tg^a_H H, quad sop H^dagger = -i H^dagger Tg^a_H gh^a . $ <eq:brst-matter>
Only the physical components transform: $sop psi_f = i gh^a Tg^a_f psi_f$ and $sop chi_f = 0$. The sources are BRST invariant. On products $sop$ acts as a graded derivation, $sop (F_1 F_2) = (sop F_1) F_2 + (-1)^(|F_1|) F_1 (sop F_2)$, with $|F_1| = 1$ for Grassmann-odd $F_1$ and $|F_1| = 0$ otherwise. It is nilpotent, $sop^2 = 0$, and
$ cal(L)_"gf" + cal(L)_"ext" = sop [ & agh^a (partial^mu A^a_mu + xi/2 b^a) - srcA^mu_a A^a_mu + srcC_a gh^a \
   & + sum_f (overline(srcF)_f Psi_f - overline(Psi)_f srcF_f) - srcH^dagger H + H^dagger srcH ] . $ <eq:gf-exact>

*Slavnov--Taylor operator.* For a functional $cal(F)$ of fields and sources,
$ ST(cal(F)) = integral dif^d x [ &
   (delta cal(F))/(delta srcA^mu_a) (delta cal(F))/(delta A^a_mu)
 + (delta cal(F))/(delta srcC_a) (delta cal(F))/(delta gh^a)
 + b^a (delta cal(F))/(delta agh^a) \
 & + sum_f ( (delta cal(F))/(delta overline(srcF)_f) (delta cal(F))/(delta Psi_f) - (delta cal(F))/(delta srcF_f) (delta cal(F))/(delta overline(Psi)_f) ) \
 & + (delta cal(F))/(delta srcH^dagger) (delta cal(F))/(delta H) + (delta_R cal(F))/(delta srcH) (delta cal(F))/(delta H^dagger) ] . $ <eq:ST>
Every derivative with respect to a Grassmann-odd quantity acts from the left, $cal(F)[X + delta X] = cal(F)[X] + integral delta X thin delta cal(F) \/ delta X$, except the single one marked $delta_R$, which acts from the right, $cal(F)[X + delta X] = cal(F)[X] + integral (delta_R cal(F) \/ delta X) thin delta X$. Spinor and group indices are contracted between the two factors of each term.

*The source of a field.* The source paired with the field $Phi$ in @eq:ST is written $Phi^*$:
$ (A^a_mu)^* = srcA^mu_a, quad (gh^a)^* = srcC_a, quad Psi_f^* = overline(srcF)_f, quad overline(Psi)_f^* = srcF_f, quad H^* = srcH^dagger, quad (H^dagger)^* = srcH . $ <eq:Phistar>
The antighost and $b^a$ have no source. For each of these six fields let $delta' cal(F) \/ delta Phi^*$ be the first factor of its term in @eq:ST: the left derivative $delta cal(F) \/ delta Phi^*$ for $Phi = A^a_mu, gh^a, Psi_f, H$; the left derivative with a minus sign, $- delta cal(F) \/ delta srcF_f$, for $Phi = overline(Psi)_f$; the right derivative $delta_R cal(F) \/ delta srcH$ for $Phi = H^dagger$. Then @eq:ST reads
$ ST(cal(F)) = integral dif^d x [ sum_Phi (delta' cal(F))/(delta Phi^*) (delta cal(F))/(delta Phi) + b^a (delta cal(F))/(delta agh^a) ], quad
  Phi in {A^a_mu, gh^a, Psi_f, overline(Psi)_f, H, H^dagger} . $ <eq:ST-compact>
The tree-level action depends on the sources only through $cal(L)_"ext"$ of @eq:Lext, which is linear in them, and the signs above are chosen such that $delta' S_0 \/ delta Phi^* = sop Phi$ for all six fields. With $sop agh^a = b^a$ and $sop b^a = 0$ the last term of @eq:ST-compact has the same form, so
$ ST(S_0) = integral dif^d x sum_Phi (delta' S_0)/(delta Phi^*) (delta S_0)/(delta Phi) + integral dif^d x thin b^a (delta S_0)/(delta agh^a)
  = integral dif^d x sum_(Phi') (sop Phi') (delta S_0)/(delta Phi') = sop S_0 , $ <eq:ST-S0>
where $Phi'$ runs over all fields, $agh^a$ and $b^a$ included.
The linearised operator $ST_(cal(F))$ of a Grassmann-even functional $cal(F)$ is defined by $ST(cal(F) + kappa thin cal(H)) = ST(cal(F)) + kappa thin ST_(cal(F)) cal(H) + O(kappa^2)$. Here $kappa$ is a constant: a number when $cal(H)$ is Grassmann even (ghost number $0$), a Grassmann-odd constant when $cal(H)$ is Grassmann odd (ghost number $1$), so that $kappa thin cal(H)$ is always Grassmann even. For the tree-level action we write $bop_d = ST_(S_0)$ in $d$ dimensions and $bop$ for its four-dimensional limit. Explicitly, $bop_d$ is the graded derivation (acting from the left, like $sop$) that acts on every field as $sop$ and replaces every source by plus or minus the classical equation of motion of its field:
$ bop_d Phi = sop Phi, quad
  bop_d srcA^mu_a = (delta S_0)/(delta A^a_mu), quad
  bop_d srcC_a = (delta S_0)/(delta gh^a), \
  bop_d overline(srcF)_f = (delta S_0)/(delta Psi_f), quad
  bop_d srcF_f = - (delta S_0)/(delta overline(Psi)_f), quad
  bop_d srcH^dagger = (delta S_0)/(delta H), quad
  bop_d srcH = - (delta S_0)/(delta H^dagger) , $ <eq:b-sources>
with left derivatives throughout. For a Grassmann-even local functional $Y$ (a counterterm action) this is the same as
$ bop_d Y = sop Y + integral dif^d x sum_Phi (delta'_Y Phi) (delta S_0)/(delta Phi) , $ <eq:b-action>
where $sop$ acts on the fields of $Y$ only, the sum runs over the six fields of @eq:ST-compact, and $delta'_Y Phi = delta' Y \/ delta Phi^*$ is obtained from $Y$ by the same derivative that gives $sop Phi$ from $S_0$: $delta'_Y A^a_mu = delta Y \/ delta srcA^mu_a$, $delta'_Y gh^a = delta Y \/ delta srcC_a$, $delta'_Y Psi_f = delta Y \/ delta overline(srcF)_f$, $delta'_Y overline(Psi)_f = - delta Y \/ delta srcF_f$, $delta'_Y H = delta Y \/ delta srcH^dagger$, $delta'_Y H^dagger = delta_R Y \/ delta srcH$. In words: a counterterm $Y$ changes the action by $Y$ and, through its source-dependent terms, changes the BRST transformation of $Phi$ by $delta'_Y Phi$; $bop_d Y$ is the first-order change of "transformation applied to action". In four dimensions $bop^2 = 0$.

*Linear identities.* The following identities hold exactly in the regularised theory, because they are linear in the quantum fields, and therefore hold for every counterterm:
$ (delta Gamma)/(delta b^a) = partial^mu A^a_mu + xi b^a, quad
  (delta Gamma)/(delta agh^a) + partial_mu (delta Gamma)/(delta srcA^mu_a) = 0 . $ <eq:linear-ids>
Here $Gamma$ is the effective action, the generating functional of the 1PI Green's functions. The first identity (gauge condition) says that $b^a$ never appears in loops. The second (ghost equation) says that $Gamma$ depends on $agh^a$ and $srcA^mu_a$ only through
$ tilde(srcA)^mu_a = srcA^mu_a + partial^mu agh^a . $ <eq:rhotilde>
Hence Green's functions with an external antighost follow from those with an external $srcA^mu_a$ and are not computed separately. In addition $S_0$ is invariant under _global_ $"SU"(3)_c times "SU"(2)_L times "U"(1)_Y$ transformations under which the sterile partners transform like their physical fields and ghosts, antighosts, $b^a$, $srcA^mu_a$ and $srcC_a$ transform in the adjoint representation. This rigid invariance is exact in $d$ dimensions, so every Green's function and every counterterm is rigidly invariant.

== The tree-level breaking $Dbrk_1$ <sec:delta>

Every term of $cal(L)_0$ except the evanescent fermion kinetic term is BRST invariant in $d$ dimensions. Applying @eq:brst-matter to @eq:Lferm, the gauge-vertex variation cancels the four-dimensional part of the kinetic-term variation and what remains is
$ Dbrk_1 = ST(S_0) = sop S_0 = sum_f integral dif^d x thin gh^a [ overline(Psi)_f thin Tg^a_f thin ghat^mu P_(overline(c)_f) thin partial_mu Psi_f
   + (partial_mu overline(Psi)_f) thin ghat^mu P_(c_f) thin Tg^a_f thin Psi_f ] . $ <eq:Delta1>
In terms of the chiral components, $Dbrk_1 = sum_f integral gh^a [ overline(psi)_f Tg^a_f ghat^mu partial_mu chi_f + (partial_mu overline(chi)_f) ghat^mu Tg^a_f psi_f ]$: it couples a ghost, a physical fermion and its sterile partner through an evanescent matrix. $Dbrk_1$ is an integrated local operator of ghost number $+1$ and dimension $5$, and $Pphys Dbrk_1 = 0$. The overall sign of @eq:Delta1 is fixed by @eq:brst-matter and the position of $gh^a$ to the left of both fermion fields. It agrees with eq. (3.40) of @Belusca-Maito:2020ala and with @Lazopoulos:2026two, and is opposite to eq. (1.18) of the starter pack @Thomsen:starterpack read with the ghost as the leftmost factor; we use @eq:Delta1.

== Feynman rules <sec:feynman-rules>

*Feynman rules and vertex functions.* All rules listed in this section are Feynman rules. For a local functional $Y$ the Feynman rule of a vertex is $i V_Y$, where the vertex function $V_Y$ is the coefficient of the product of the plane-wave polarisations of the legs in $Y$ (defined precisely in @sec:dictionary). The symbol $V$ always denotes the vertex function, never the Feynman rule.

Fields are Fourier transformed with $e^(-i p dot x)$ for an incoming momentum $p$. A vertex is $i$ times the corresponding term of the action, with $partial_mu -> -i p_mu$ acting on a field of incoming momentum $p$; a propagator is $i$ times the inverse of the quadratic operator. A closed fermion loop or ghost loop carries a factor $-1$. The sum of all 1PI diagrams with given external legs equals $i Gamma$ for these legs, with $Gamma$ the effective action. All momenta at a vertex are incoming unless stated otherwise.

*Propagators* (momentum $k$ in the direction of the arrow for $Psi_f$, $gh$, $H$):
$ Psi_f: quad (i sl(k))/(k^2) , quad quad
  A^a_mu A^b_nu: quad (-i g_(mu nu) delta^(a b))/(k^2), quad quad
  gh^a agh^b: quad (i delta^(a b))/(k^2), quad quad
  H_i H^dagger_j: quad (i delta_(i j))/(k^2 - m^2) . $ <eq:propagators>
In the fermion and gauge propagators $sl(k)$ and $g_(mu nu)$ are $d$-dimensional; $k^2$ is always $d$-dimensional. Colour and $"SU"(2)_L$ indices of $Psi_f$ are conserved along the line.

*Vertices of the gauge sector.*
$ overline(Psi)_f Psi_f A^a_mu: quad i thin gbar^mu P_(c_f) thin Tg^a_f , $ <eq:V-ffv>
$ A^a_mu (p_1) A^b_nu (p_2) A^c_rho (p_3): quad fg^(a b c) [ g^(mu nu) (p_1 - p_2)^rho + g^(nu rho) (p_2 - p_3)^mu + g^(rho mu) (p_3 - p_1)^nu ] , $ <eq:V-vvv>
$ A^a_mu A^b_nu A^c_rho A^e_sigma: quad
  -i [ & fg^(a b x) fg^(c e x) (g^(mu rho) g^(nu sigma) - g^(mu sigma) g^(nu rho)) \
     & + fg^(a c x) fg^(b e x) (g^(mu nu) g^(rho sigma) - g^(mu sigma) g^(nu rho)) \
     & + fg^(a e x) fg^(b c x) (g^(mu nu) g^(rho sigma) - g^(mu rho) g^(nu sigma)) ] , $ <eq:V-vvvv>
$ agh^a (p' "outgoing") thin A^b_mu thin gh^c: quad - fg^(a b c) thin p'^mu . $ <eq:V-ghost>

*Vertices of the Higgs sector* ($H_j$ incoming with momentum $p$, $H^dagger_i$ outgoing with momentum $p'$):
$ H^dagger_i H_j A^a_mu: quad i (p + p')^mu (Tg^a_H)_(i j), quad quad
  H^dagger_i H_j A^a_mu A^b_nu: quad i g^(mu nu) {Tg^a_H, Tg^b_H}_(i j), $ <eq:V-ssv>
$ H^dagger_i H_j H^dagger_k H_l: quad -2 i lambda (delta_(i j) delta_(k l) + delta_(i l) delta_(k j)) . $ <eq:V-ssss>

*Yukawa vertices* (the scalar leg is labelled by the field operator, $H$ or $H^dagger$, that stands in the term of @eq:LYuk; $j, k$ are doublet indices, $alpha, beta$ colour indices):
$ & overline(Q)^(alpha j) U_beta H^dagger_k: quad -i y_u epsilon_(j k) delta^alpha_beta PR, quad quad
  && overline(U)^alpha Q_(beta j) H_k: quad -i y_u epsilon_(j k) delta^alpha_beta PL, \
  & overline(Q)^(alpha j) D_beta H_k: quad -i y_d delta_(j k) delta^alpha_beta PR, quad quad
  && overline(D)^alpha Q_(beta j) H^dagger_k: quad -i y_d delta_(j k) delta^alpha_beta PL, \
  & overline(L)^j E H_k: quad -i y_e delta_(j k) PR, quad quad
  && overline(E) L_j H^dagger_k: quad -i y_e delta_(j k) PL . $ <eq:V-yuk>

*Vertices with BRST sources* (the source is an external leg that carries momentum but has no propagator):
$ srcA^mu_a thin gh^b (p): quad delta^(a b) p_mu, quad quad
  srcA^mu_a A^b_nu gh^c: quad i fg^(a b c) g_mu^(space nu), quad quad
  srcC_a gh^b gh^c: quad -i fg^(a b c), $ <eq:V-src-gauge>
$ overline(srcF)_f gh^a Psi_f: quad - Tg^a_f P_(c_f), quad quad
  gh^a overline(Psi)_f srcF_f: quad - P_(overline(c)_f) Tg^a_f, quad quad
  srcH^dagger gh^a H: quad - Tg^a_H, quad quad
  gh^a H^dagger srcH: quad + Tg^a_H . $ <eq:V-src-matter>
Each factor is $i$ times the coefficient of the corresponding term of @eq:Lext and @eq:Lext-explicit with the Grassmann-odd fields in the order written; the relative sign of a diagram with ghost and source legs follows from Wick's theorem applied to these terms.

*The breaking vertex.* From @eq:Delta1, with $Psi_f$ incoming with momentum $p$, $overline(Psi)_f$ outgoing with momentum $p'$ and the ghost $gh^a$ carrying the incoming momentum $p' - p$, the Feynman rule is
$ i V_(Dbrk_1)[gh^a overline(Psi)_f Psi_f] = Tg^a_f ( slh(p) thin P_(overline(c)_f) - slh(p)' thin P_(c_f) ) , quad
  V_(Dbrk_1)[gh^a overline(Psi)_f Psi_f] = -i thin Tg^a_f ( slh(p) thin P_(overline(c)_f) - slh(p)' thin P_(c_f) ) . $ <eq:V-delta>
The vertex function follows from @eq:Delta1 with $partial_mu Psi_f -> -i p_mu Psi_f$ and $partial_mu overline(Psi)_f -> +i p'_mu overline(Psi)_f$; the generator $Tg^a_f$ commutes with the Dirac matrices. Diagrams with this vertex sum to $i thin (Dbrk_1 dot Gamma)$, the insertion of $Dbrk_1$ into the effective action (@sec:m-objects).

*Consequences for fermion lines.* Every vertex except @eq:V-delta contains only physical components, $P_(c_f) Psi_f$ and $overline(Psi)_f P_(overline(c)_f)$. A fermion propagator between two such vertices is therefore sandwiched as
$ P_(c_f) (i sl(k))/(k^2) P_(overline(c)_f) = (i slb(k))/(k^2) P_(overline(c)_f) : $ <eq:phys-prop>
its numerator is four-dimensional while its denominator is $d$-dimensional. In a diagram without the breaking vertex the whole Dirac algebra is thus four-dimensional, where $g5$ anticommutes with every matrix present, and the regularisation enters only through the integration rules of @sec:integrals. The sterile components propagate only next to a breaking vertex, through $P_(c_f) (i sl(k)\/k^2) P_(c_f) = i slh(k) P_(c_f) \/ k^2$.

== One-loop integrals and the expansion operator <sec:integrals>

*Measure.* $intk = mu^(2 ep) integral dif^d k \/ (2 pi)^d$, and $overline("MS")$ means that $mu^2$ stands for $mu^2 e^(gamma_E) \/ (4 pi)$. The operator $Kop$ extracts the pole part in $ep$ of a Laurent series.

*Expansion operator.* Let $G$ be a one-loop 1PI diagram with loop momentum $k$, external momenta $p_i$ and internal masses $m_j$. Its superficial degree of divergence $dg(G)$ is the mass dimension of the integral at $d = 4$: if the numerator is homogeneous of degree $n$ in $k$, $p_i$ and $m_j$ together and there are $N$ propagators, $dg(G) = 4 + n - 2 N$. Each denominator is rewritten identically around an auxiliary mass $MUV$,
$ 1/((k + p)^2 - m^2) = 1/(k^2 - MUV^2) sum_(j >= 0) ( (-2 k dot p - p^2 + m^2 - MUV^2)/(k^2 - MUV^2) )^j , $ <eq:uv-expansion>
and $Top(dg(G)) G$ is the sum of all terms of the expanded integrand whose total degree in the small quantities is at most $dg(G)$; when the degree is known it is written as a number, $Top(2)$, $Top(1)$, $Top(0)$. The small quantities are the $p_i$ and $m_j$, wherever they stand (in the numerator of the diagram or in the numerators of @eq:uv-expansion), and the $MUV$ in the numerators of @eq:uv-expansion; each counts as degree one. Equivalently, a term is kept if and only if its vacuum integral $intk (k^2)^a \/ (k^2 - MUV^2)^b$ has $4 + 2 a - 2 b >= 0$. $Top(dg(G)) G$ is a combination of vacuum integrals with the single scale $MUV$, and it contains the complete ultraviolet behaviour of $G$.

*Vacuum integrals.*
$ intk ((k^2)^a)/((k^2 - MUV^2)^b) = (i (-1)^(a+b))/((4 pi)^(d\/2)) (Gamma(a + d\/2) Gamma(b - a - d\/2))/(Gamma(d\/2) Gamma(b)) (MUV^2)^(d\/2 + a - b) (mu^2)^ep . $ <eq:tadpole>
For the logarithmically divergent case $b = a + 2$ the pole is universal,
$ Kop intk ((k^2)^a)/((k^2 - MUV^2)^(a+2)) = ilf 1/ep , $ <eq:logpole>
and for $b = a + 1$ and $b = a$ the poles are $ilf thin (a+1) MUV^2 \/ ep$ and $ilf thin a (a+1) MUV^4 \/ (2 ep)$. All dependence on $MUV$ cancels in $Kop Top(dg(G)) G$.

*Tensor reduction.* Under $intk$ with a denominator depending on $k^2$ only, odd powers of $k$ vanish and
$ k^mu k^nu -> (k^2)/d g^(mu nu), quad
  k^mu k^nu k^rho k^sigma -> ((k^2)^2)/(d (d+2)) (g^(mu nu) g^(rho sigma) + g^(mu rho) g^(nu sigma) + g^(mu sigma) g^(nu rho)), $ <eq:tensor-red>
and for six momenta $(k^2)^3 \/ (d (d+2)(d+4))$ times the sum of the fifteen pairings. The metric is the full $d$-dimensional one. For the projected momenta this gives
$ bar(k)^mu bar(k)^nu -> (k^2)/d gmb^(mu nu), quad hat(k)^mu hat(k)^nu -> (k^2)/d gmh^(mu nu), quad bar(k)^mu hat(k)^nu -> 0, quad
  bar(k)^2 -> 4/d k^2, quad hat(k)^2 -> (d-4)/d k^2 , $ <eq:tensor-red-split>
and the analogous projections of the four- and six-momentum formulas.

*The evanescent tag $eta$.* The trace of the evanescent metric is written
$ gmh^mu_(space mu) = d - 4 = -2 eta, $ <eq:eta>
with $eta$ kept as a separate symbol while poles in $ep$ are extracted, and set to $eta = ep$ only at the very end. This separates the two roles of $d - 4$: a factor $eta$ comes from an evanescent contraction in the numerator, a pole $1\/ep$ from an ultraviolet divergence. The rule is:

- $eta$ is produced only by the complete contraction of hatted objects: $gmh^mu_(space mu) = -2 eta$, $ghat^mu ghat_mu = -2 eta$, and under the integral $hat(k)^2 -> -2 eta thin k^2 \/ d$;
- an object written with a hat is never re-expressed through unsplit objects (for instance $hat(k)^2$ is not replaced by $k^2 - bar(k)^2$);
- every other factor of $d$ (those of @eq:tadpole and @eq:tensor-red, and those from contractions of unsplit $d$-dimensional objects, $g^mu_(space mu) = d$) is written $d = 4 - 2 ep$.

The tag is a bookkeeping device for the finite terms of the form (evanescent contraction) $times$ (pole). It is used in two places only, and in both the result does not depend on how the unsplit $d$-dimensional objects are treated: the residue $N_0$ of the pole of any diagram, and the coefficient $N_1$ of diagrams with the breaking vertex @eq:V-delta after projection onto four-dimensional external momenta and indices, where the whole diagram is proportional to $eta$ (@sec:method).

*What is extracted from a one-loop diagram.* After the expansion @eq:uv-expansion, the tensor reduction @eq:tensor-red and the algebra of @sec:bmhv, every diagram takes the form
$ Top(dg(G)) G = ilf ( 1/ep N_0 + eta/ep N_1 ) + "terms that are finite and depend on" MUV , $ <eq:N0N1>
where $N_0$ and $N_1$ are polynomials in the external momenta and masses, taken at $ep = 0$. Two quantities are read off:
$ "pole part:" & quad Kop Top(dg(G)) G = ilf thin N_0 / ep, \
  "evanescent-induced finite part:" & quad ilf thin N_1 quad "at" eta = ep . $ <eq:extraction>
$N_0$ may contain evanescent structures (hatted momenta or indices of external fields); they are kept. $N_1$ is local and independent of $MUV$. A diagram without the breaking vertex has $N_1 = 0$: its Dirac algebra is four-dimensional (@eq:phys-prop) and no closed evanescent contraction occurs. The remaining finite terms are not local counterterm data and are never used.
