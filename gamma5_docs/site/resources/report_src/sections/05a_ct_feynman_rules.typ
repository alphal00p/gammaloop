#import "../macros.typ": *

== Feynman rules of the counterterms <sec:ct-rules>

This subsection turns the counterterms of @sec:Sren1 into Feynman rules: the singular counterterm $delta S^((1))$ of @eq:c-dS-inv and @eq:c-dS-ev, and the finite counterterm $Delta S^((1))_"ST"$ of @eq:c-DS. Every rule below can be entered into a model as it stands.

*Conventions.* All momenta are incoming, except that $Psi_f$ and $H_j$ enter with momentum $p$ and $overline(Psi)_f$ and $H^dagger_i$ leave with momentum $p'$, as in @sec:feynman-rules. The Feynman rule of a local functional $Y$ is $i V_Y$, with $V_Y$ the vertex function of @sec:dictionary: the coefficient of the product of the leg polarisations, with $partial_mu -> -i p_mu$ on a leg that enters with $p$ and $partial_mu -> +i p'_mu$ on a leg that leaves with $p'$, summed over all assignments of identical legs. Grassmann-odd legs are taken in the order in which they stand in the monomial, as in @eq:V-src-matter: $(overline(Psi)_f, Psi_f)$, $(agh, gh)$, $(gh, overline(Psi)_f)$.

With this normalisation the counterterm rule of $delta S^((1))$ is minus the pole of the one-loop diagrams with the same legs, $i V_(delta S^((1))) = - sum_G Kop Top(dg(G)) G$ (@eq:m-Z1). The rules of $delta S^((1))$ carry the factor $lf thin 1\/ep$ through the constants $delta Z_x$, $delta y_k$, $delta lambda$, $delta m^2$ of @tab:c-Z. The rules of $Delta S^((1))_"ST"$ carry the factor $lf$ explicitly.

=== Singular counterterm: invariant part <cr-inv>

*How the factors follow from @eq:c-bare.* $delta S^((1))_"inv"$ is the first-order term of $S_0$ written in bare quantities. A monomial of $S_0$ is multiplied by $Z_x^(1\/2)$ for each leg of a field $x$ ($Z_A^(-1\/2)$ for $agh$ and $srcA$, $Z_gh^(-1\/2)$ for $srcC$, $Z_f^(-1\/2)$ for $srcF_f$, $overline(srcF)_f$, $Z_H^(-1\/2)$ for $srcH$, $srcH^dagger$), by $Z_(c_r)$ for each power of the coupling $c_r$ in $Tg^a$ or $fg^(a b c)$, and by $1 + (delta y_k)\/y_k$ or $1 + (delta lambda)\/lambda$. The counterterm rule of a vertex is therefore its tree rule times the first-order factor
$ sum_"legs" plus.minus 1/2 delta Z_x + sum_"couplings" delta Z_(c_r) + (delta y_k)/y_k thin "or" thin (delta lambda)/lambda , $ <cr-factor>
with $+$ for fields and $-$ for $agh$ and the sources. There are three exceptions. (i) The gauge fixing is not renormalised ($b_0 partial dot A_0 = b thin partial dot A$ and $xi_0 b_0^2 = xi b^2$), so the gauge two-point counterterm is transverse. (ii) The fermion kinetic counterterm contains only the physical components: the evanescent kinetic term of @eq:Lferm-weyl is not rescaled ($chi_(f 0) = chi_f$). (iii) The Higgs two-point function also receives $delta m^2$. @cr-tab-inv lists all rules. The tree rules are those of @sec:feynman-rules with unsplit $d$-dimensional metrics and momenta; the fermion vertices contain $gbar^mu$ already.

#figure(
  {
    set par(justify: false)
    set text(size: 8.6pt)
    table(
    columns: (auto, auto, auto, 1fr),
    align: (left, left, left, left),
    inset: (x: 3.5pt, y: 4pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else if y in (5, 9, 13) { (top: 0.3pt + luma(160)) } else { none },
    table.header([legs], [rule with factor $z$], [$z$ from @eq:c-bare], [$16 pi^2 ep thin z$: general; Standard Model]),
    [$A^a_mu (p) A^b_nu (-p)$], [$-i z thin delta^(a b) (p^2 g^(mu nu) - p^mu p^nu)$], [$delta Z^((r))_A$], [$5/3 C^((r))_"ad" - 2/3 S^((r)) - 1/3 S^((r))_H$; $11/3 gs^2$, $11/6 gw^2$, $-43/18 gy^2$],
    [$agh^a (p) gh^b (p)$], [$i z thin p^2 delta^(a b)$], [$1/2 (delta Z^((r))_gh - delta Z^((r))_A)$], [$1/2 C^((r))_"ad"$; $3/2 gs^2$, $gw^2$, $0$],
    [$overline(Psi)_f (p) Psi_f (p)$], [$i z thin slb(p) P_(c_f)$], [$delta Z_f$], [$-(C_f + cal(Y)_f)$; @tab:c-Z],
    [$H^dagger_i (p) H_j (p)$], [$i [delta Z_H (p^2 - m^2) - delta m^2] delta_(i j)$], [$delta Z_H$, $delta m^2$], [$2 C_H - tilde(Y)$, $m^2 (6 lambda - 3 C_H + tilde(Y))$],
    [$A A A$], [$z times$ @eq:V-vvv], [$delta Z_(c_r) + 3/2 delta Z^((r))_A$], [$2/3 C^((r))_"ad" - 2/3 S^((r)) - 1/3 S^((r))_H$; $2/3 gs^2$, $-1/6 gw^2$],
    [$A A A A$], [$z times$ @eq:V-vvvv], [$2 delta Z_(c_r) + 2 delta Z^((r))_A$], [$-1/3 C^((r))_"ad" - 2/3 S^((r)) - 1/3 S^((r))_H$; $-7/3 gs^2$, $-13/6 gw^2$],
    [$agh^a A^b_mu gh^c$], [$z times$ @eq:V-ghost], [$delta Z_(c_r) + 1/2 delta Z^((r))_gh$], [$-1/2 C^((r))_"ad"$; $-3/2 gs^2$, $-gw^2$],
    [$overline(Psi)_f Psi_f A^a_mu$], [$z times$ @eq:V-ffv], [$delta Z_f + delta Z_(c_r) + 1/2 delta Z^((r))_A$], [$-(C_f + cal(Y)_f + C^((r))_"ad")$],
    [$H^dagger H A^a_mu$], [$z times$ @eq:V-ssv], [$delta Z_H + delta Z_(c_r) + 1/2 delta Z^((r))_A$], [$2 C_H - tilde(Y) - C^((r))_"ad"$],
    [$H^dagger H A^a_mu A^b_nu$], [$z times$ @eq:V-ssv], [$delta Z_H + sum_(t = r, s) (delta Z_(c_t) + 1/2 delta Z^((t))_A)$], [$2 C_H - tilde(Y) - C^((r))_"ad" - C^((s))_"ad"$],
    [$H^dagger H H^dagger H$], [$z times$ @eq:V-ssss], [$(delta lambda)\/lambda + 2 delta Z_H$], [$12 lambda - 2 C_H - (Y_4 - 3/16 G_4)\/lambda$],
    [Yukawa, $O_k$, $O_k^dagger$], [$z times$ @eq:V-yuk], [$(delta y_k)\/y_k + 1/2 (delta Z_f + delta Z_(f') + delta Z_H)$], [$u$: $-(16/3 gs^2 + 3/4 gw^2 + 25/36 gy^2 + y_d^2)$ \ $d$: $-(16/3 gs^2 + 3/4 gw^2 + 1/36 gy^2 + y_u^2)$ \ $e$: $-(3/4 gw^2 + 9/4 gy^2)$],
    [$srcA^mu_a gh^b$], [$z times$ @eq:V-src-gauge], [$1/2 (delta Z^((r))_gh - delta Z^((r))_A)$], [$1/2 C^((r))_"ad"$; $3/2 gs^2$, $gw^2$, $0$],
    [$srcA A gh$, $srcC gh gh$], [$z times$ @eq:V-src-gauge], [$delta Z_(c_r) + 1/2 delta Z^((r))_gh$], [$-1/2 C^((r))_"ad"$; $-3/2 gs^2$, $-gw^2$],
    [$overline(srcF)_f gh Psi_f$, $gh overline(Psi)_f srcF_f$, \ $srcH^dagger gh H$, $gh H^dagger srcH$], [$z times$ @eq:V-src-matter], [$delta Z_(c_r) + 1/2 delta Z^((r))_gh$], [$-1/2 C^((r))_"ad"$; $-3/2 gs^2$, $-gw^2$, $0$],
  )},
  caption: [Feynman rules of $delta S^((1))_"inv"$ (@eq:c-dS-inv). $r$ is the gauge factor of the gauge leg $A^a$ or of the ghost leg; in $H^dagger H A^a A^b$ the leg $A^b$ belongs to the factor $s$. In the two-point rules $agh^a$, $overline(Psi)_f$ and $H^dagger_i$ leave with the momentum $p$. Where several values are given they are $r = G, W, B$ (in the rows $A A A$, $A A A A$, $agh A gh$ and $srcA A gh$, $srcC gh gh$ only $r = G, W$ exist). The fermion and Yukawa constants are those of @tab:c-Z; $f$ and $f'$ are the two fermion multiplets of the Yukawa monomial. Every ghost and source rule of the abelian factor vanishes, because $Z_(c_B) (Z^((B))_A)^(1\/2) = 1$ and $delta Z^((B))_gh = -2 delta Z_(c_B)$.],
) <cr-tab-inv>

The momentum dependence of each row is that of the tree rule. For example, the Yukawa row gives for $overline(Q)^(alpha i) U_beta H^dagger_j$ the rule $-i y_u epsilon_(i j) delta^alpha_beta PR (lf thin 1\/ep) (-16/3 gs^2 - 3/4 gw^2 - 25/36 gy^2 - y_d^2)$, and the quartic Higgs row gives $-2 i (delta lambda + 2 lambda delta Z_H)(delta_(i j) delta_(k l) + delta_(i l) delta_(k j))$ with $16 pi^2 ep (delta lambda + 2 lambda delta Z_H) = 12 lambda^2 - 2 lambda C_H - Y_4 + 3/16 G_4$.

=== Singular counterterm: evanescent part <cr-ev>

$delta S^((1))_"ev"$ of @eq:c-dS-ev multiplies $E^((r))_A$ by $-1/3 S^((r)) lf thin 1\/ep$, $tilde(E)^((r))_F$ by $-2/3 S^((r)) lf thin 1\/ep$, $E_H$ by $-2/3 tilde(Y) lf thin 1\/ep$ and $tilde(E)_H$ by $-tilde(Y) lf thin 1\/ep$. Because $g^(mu nu) = gmb^(mu nu) + gmh^(mu nu)$, the Higgs operator is $tilde(E)_H = - integral dif^d x thin (D_mu H)^dagger gmh^(mu nu) D_nu H$. Write $Lambda^(mu nu alpha)_3$ and $Lambda^(mu nu alpha beta)_4$ for the tree Feynman rules @eq:V-vvv and @eq:V-vvvv of the gauge factor $r$ (with $alpha, beta$ in place of $rho, sigma$), and $overline(Lambda)_3$, $overline(Lambda)_4$ for the same expressions with $g -> gmb$ and every momentum $p_i -> bar(p)_i$. The rules are, with $a, b, c, d$ in the factor $r$:
$ E^((r))_A: quad & A^a_mu (p) A^b_nu (-p): quad i/3 S^((r)) lf 1/ep thin hat(p)^2 gmb^(mu nu) delta^(a b) , \
  tilde(E)^((r))_F: quad & A^a_mu (p) A^b_nu (-p): quad (2 i)/3 S^((r)) lf 1/ep thin delta^(a b) [ (bar(p)^2 gmb^(mu nu) - bar(p)^mu bar(p)^nu) - (p^2 g^(mu nu) - p^mu p^nu) ] , \
  & A^a_mu A^b_nu A^c_alpha: quad -2/3 S^((r)) lf 1/ep thin (overline(Lambda)^(mu nu alpha)_3 - Lambda^(mu nu alpha)_3) , \
  & A^a_mu A^b_nu A^c_alpha A^d_beta: quad -2/3 S^((r)) lf 1/ep thin (overline(Lambda)^(mu nu alpha beta)_4 - Lambda^(mu nu alpha beta)_4) , \
  E_H: quad & H^dagger_i (p) H_j (p): quad -(2 i)/3 tilde(Y) lf 1/ep thin hat(p)^2 delta_(i j) , \
  tilde(E)_H: quad & H^dagger_i (p) H_j (p): quad i tilde(Y) lf 1/ep thin hat(p)^2 delta_(i j) , \
  & H^dagger_i H_j A^a_mu: quad i tilde(Y) lf 1/ep thin (hat(p) + hat(p)')^mu (Tg^a_H)_(i j) , \
  & H^dagger_i H_j A^a_mu A^b_nu: quad i tilde(Y) lf 1/ep thin gmh^(mu nu) {Tg^a_H, Tg^b_H}_(i j) . $ <cr-ev-rules>
For $r = B$ only the two-point rules exist, since $fg^(a b c) = 0$. Each rule contains a hatted momentum or metric, or the difference between a barred and an unsplit structure: all of them vanish under $Pphys$. In the Higgs rules the gauge legs belong to $W$ or $B$.

*Complete rules.* Adding @cr-ev-rules to @cr-tab-inv gives @cr-tab-full. The closed fermion loops give barred structures with the coefficient $-2/3 S^((r))$ or $-tilde(Y)$; all other loops give the unsplit structure. In the gauge two-point function the unsplit transverse tensor multiplies $delta Z^((r))_A + 2/3 S^((r)) lf thin 1\/ep$, the barred one $-2/3 S^((r)) lf thin 1\/ep$, and $hat(p)^2 gmb^(mu nu)$ is the non-transverse evanescent remainder of the fermion loop.

#figure(
  {
    set par(justify: false)
    set text(size: 8.6pt)
    table(
    columns: (auto, 1fr),
    align: (left, left),
    inset: (x: 4pt, y: 4.5pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([legs], [complete rule of $delta S^((1))$, in units of $lf thin 1\/ep$]),
    [$A^a_mu (p) A^b_nu (-p)$], [$-i delta^(a b) [ (5/3 C^((r))_"ad" - 1/3 S^((r))_H) (p^2 g^(mu nu) - p^mu p^nu) - 2/3 S^((r)) (bar(p)^2 gmb^(mu nu) - bar(p)^mu bar(p)^nu) - 1/3 S^((r)) hat(p)^2 gmb^(mu nu) ]$ \
      $G$: $(5, -4/3, -2/3) gs^2$; $W$: $(19/6, -4/3, -2/3) gw^2$; $B$: $(-1/6, -20/9, -10/9) gy^2$ for the three brackets],
    [$A^a_mu A^b_nu A^c_alpha$], [$(2/3 C^((r))_"ad" - 1/3 S^((r))_H) Lambda^(mu nu alpha)_3 - 2/3 S^((r)) overline(Lambda)^(mu nu alpha)_3$; #h(0.4em) $G$: $2 gs^2 Lambda_3 - 4/3 gs^2 overline(Lambda)_3$; $W$: $7/6 gw^2 Lambda_3 - 4/3 gw^2 overline(Lambda)_3$],
    [$A^a_mu A^b_nu A^c_alpha A^d_beta$], [$-1/3 (C^((r))_"ad" + S^((r))_H) Lambda^(mu nu alpha beta)_4 - 2/3 S^((r)) overline(Lambda)^(mu nu alpha beta)_4$; #h(0.4em) $G$: $-gs^2 Lambda_4 - 4/3 gs^2 overline(Lambda)_4$; $W$: $-5/6 gw^2 Lambda_4 - 4/3 gw^2 overline(Lambda)_4$],
    [$H^dagger_i (p) H_j (p)$], [$i [ 2 C_H p^2 - tilde(Y) bar(p)^2 - 2/3 tilde(Y) hat(p)^2 - m^2 (6 lambda - C_H) ] delta_(i j)$],
    [$H^dagger_i H_j A^a_mu$], [$i [ (2 C_H - C^((r))_"ad") (p + p')^mu - tilde(Y) (bar(p) + bar(p)')^mu ] (Tg^a_H)_(i j)$],
    [$H^dagger_i H_j A^a_mu A^b_nu$], [$i [ (2 C_H - C^((r))_"ad" - C^((s))_"ad") g^(mu nu) - tilde(Y) gmb^(mu nu) ] {Tg^a_H, Tg^b_H}_(i j)$],
  )},
  caption: [Rules of the complete singular counterterm $delta S^((1)) = delta S^((1))_"inv" + delta S^((1))_"ev"$ for the vertices with an evanescent part; all other rules are those of @cr-tab-inv. $Lambda_3$, $Lambda_4$ and $overline(Lambda)_3$, $overline(Lambda)_4$ are defined before @cr-ev-rules. $a$ is in the factor $r$, $b$ in the factor $s$.],
) <cr-tab-full>

#check[Every rule of @cr-tab-inv and @cr-tab-full equals minus the pole of the class that has the same legs: @sec:g0_VV, @sec:g0_VVV, @sec:g0_VVVV, @sec:g0_FF, @sec:g0_FFV, @sec:g0_SS, @sec:g0_SSV_SSVV, @sec:g0_SSSS, @sec:g0_FFS and @sec:g0_ghost_sources, for all gauge factors, multiplets and Yukawa couplings (exact comparison of the coefficients). For example, the barred and unsplit parts of the gluon three-point rule, $2 gs^2 Lambda_3 - 4/3 gs^2 overline(Lambda)_3$, are minus the gauge-and-ghost and the quark parts of the pole of @sec:g0_VVV. The rules @cr-ev-rules were obtained independently by functional differentiation of the operators @eq:evanescent on plane waves in $d = 6$ with two explicit hatted components (@sec:supporting).]

=== Finite counterterm <cr-fin>

$Delta S^((1))_"ST" = lf sum_alpha r_alpha Q_alpha$ (@eq:c-DS, @tab:sts-r). Each monomial has the rule $i lf r_alpha V_(Q_alpha)$. All objects are four-dimensional: in a $d$-dimensional diagram every metric is $gmb$, every momentum is its barred part and every Dirac matrix is $gbar^mu$ (@sec:Sren1, continuation to $d$ dimensions). The rules of the 39 monomials with $r_alpha != 0$ are the following.

*Gauge two-point functions* ($Q^((r))_L$, $r_alpha = -1/3 S^((r))$), for $r = G, W, B$:
$ A^a_mu (p) A^b_nu (-p): quad i/3 S^((r)) lf thin delta^(a b) bar(p)^mu bar(p)^nu , quad quad 1/3 S^((r)) = 2/3 gs^2, thin 2/3 gw^2, thin 10/9 gy^2 . $ <cr-fin-2>
This is the only longitudinal counterterm. It cancels the finite longitudinal part that the closed fermion loops give to the renormalised gauge self-energy (@sec:st-checks).

*Gauge three-point functions* ($Q^((r))_(3 f)$, $r_alpha = -1/3 S^((r))$), for $r = G, W$:
$ A^a_mu (p_1) A^b_nu (p_2) A^c_alpha (p_3): quad -1/3 S^((r)) lf thin fg^(a b c) [ gmb^(mu nu) (bar(p)_1 - bar(p)_2)^alpha + gmb^(nu alpha) (bar(p)_2 - bar(p)_3)^mu + gmb^(alpha mu) (bar(p)_3 - bar(p)_1)^nu ] , $ <cr-fin-3>
i.e. $-1/3 S^((r)) lf thin overline(Lambda)_3$, with $-1/3 S^((G)) = -2/3 gs^2$ and $-1/3 S^((W)) = -2/3 gw^2$. There is no cubic rule for $B$ and no mixed cubic rule. Together with @cr-fin-2 this is the rule of the first two terms of @eq:c-DS, $sum_f Tr_f [1/6 (partial dot Am)^2 - i/3 (partial_mu Am_nu) [Am^mu, Am^nu]]$.

*Gauge four-point functions.* For legs $A^(a_i)_(mu_i) (p_i)$, $i = 1, dots, 4$, the quartic part of @eq:c-DS has the rule
$ i lf sum_(pi in S_4) sum_f Tr_f (Tg^(a_(pi(1))) Tg^(a_(pi(2))) Tg^(a_(pi(3))) Tg^(a_(pi(4)))) [ 1/3 gmb^(mu_(pi(1)) mu_(pi(2))) gmb^(mu_(pi(3)) mu_(pi(4))) - 1/4 gmb^(mu_(pi(1)) mu_(pi(3))) gmb^(mu_(pi(2)) mu_(pi(4))) ] , $ <cr-fin-4>
summed over the 24 permutations $pi$ of the four legs, with $Tg^a$ in the representation of the multiplet $f$. In the basis of @tab:Qbasis the same rule is the sum of the rules $i lf r_alpha V_(Q_alpha)$ of @cr-tab-quartic. There $sum_((i j)(k l))$ runs over the three pairings $(12)(34)$, $(13)(24)$, $(14)(23)$ of the legs, and $e$ is summed over the eight gluons.

#figure(
  {
    set par(justify: false)
    set text(size: 8.6pt)
    table(
    columns: (auto, auto, auto, 1fr),
    align: (left, left, left, left),
    inset: (x: 4pt, y: 4.5pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$Q_alpha$], [$r_alpha$], [legs], [$V_(Q_alpha)$]),
    [$Q_(G 4 a)$, $Q_(W 4 a)$], [$11/72$, $7/24$], [four $G$; four $W$], [$8 c_r^4 sum_((i j)(k l)) delta^(a_i a_j) delta^(a_k a_l) gmb^(mu_i mu_j) gmb^(mu_k mu_l)$, $c_r = gs, gw$],
    [$Q_(G 4 b)$, $Q_(W 4 b)$], [$-1/4$, $-1/4$], [four $G$; four $W$], [$4 c_r^4 sum_((i j)(k l)) gmb^(mu_i mu_j) gmb^(mu_k mu_l) (delta^(a_i a_k) delta^(a_j a_l) + delta^(a_i a_l) delta^(a_j a_k))$],
    [$Q_(G 4 c)$], [$5/12$], [four $G$], [$8 gs^4 sum_((i j)(k l)) d^(a_i a_j e) d^(a_k a_l e) gmb^(mu_i mu_j) gmb^(mu_k mu_l)$],
    [$Q_(B 4)$], [$95/648$], [four $B$], [$8 gy^4 sum_((i j)(k l)) gmb^(mu_i mu_j) gmb^(mu_k mu_l)$],
    [$Q^((r s))_a$], [$1/24$, $11/216$, $1/36$], [$1, 2 in r$; $3, 4 in s$], [$4 c_r^2 c_s^2 delta^(a_1 a_2) delta^(a_3 a_4) gmb^(mu_1 mu_2) gmb^(mu_3 mu_4)$],
    [$Q^((r s))_b$], [$1/12$, $11/108$, $1/18$], [$1, 2 in r$; $3, 4 in s$], [$2 c_r^2 c_s^2 delta^(a_1 a_2) delta^(a_3 a_4) (gmb^(mu_1 mu_3) gmb^(mu_2 mu_4) + gmb^(mu_1 mu_4) gmb^(mu_2 mu_3))$],
    [$Q_(G G G B)$], [$1/18$], [$1, 2, 3 in G$; $4 = B$], [$2 gs^3 gy thin d^(a_1 a_2 a_3) sum_((i j)(k l)) gmb^(mu_i mu_j) gmb^(mu_k mu_l)$],
  )},
  caption: [Vertex functions of the quartic gauge monomials of $Delta S^((1))_"ST"$; the rule is $i lf r_alpha V_(Q_alpha)$. For $Q^((r s))_a$ and $Q^((r s))_b$ the three values of $r_alpha$ are $r s = G W, G B, W B$. $d^(a b c)$ is the $"SU"(3)$ symmetric tensor of @sec:Qbasis.],
) <cr-tab-quartic>

*Fermion--gauge vertices* ($Q^(f r)_V$, $r_alpha = -C_f - 1/2 nu_(f r)$), for every multiplet $f$ and every gauge factor $r$ that acts on it:
$ overline(Psi)_f (p') Psi_f (p) A^a_mu: quad -i (C_f + 1/2 nu_(f r)) lf thin gbar^mu Tg^a_f P_(c_f) quad (a in r) , $ <cr-fin-ffv>
with the Standard-Model values of $C_f + 1/2 nu_(f r)$ in @cr-tab-ffv.

#figure(
  {
    set text(size: 8.6pt)
    table(
    columns: 4,
    align: (left, left, left, left),
    inset: (x: 4pt, y: 4pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$(f, r)$], [$C_f + 1/2 nu_(f r)$], [$(f, r)$], [$C_f + 1/2 nu_(f r)$]),
    [$(Q, G)$], [$4/3 gs^2 + 3/4 gw^2 + 1/36 gy^2 + 1/2 (y_u^2 + y_d^2)$], [$(Q, W)$], [$4/3 gs^2 + 3/4 gw^2 + 1/36 gy^2$],
    [$(Q, B)$], [$4/3 gs^2 + 3/4 gw^2 + 1/36 gy^2 + 2 y_u^2 - y_d^2$], [$(U, G)$], [$4/3 gs^2 + 4/9 gy^2 + y_u^2$],
    [$(U, B)$], [$4/3 gs^2 + 4/9 gy^2 + 1/4 y_u^2$], [$(D, G)$], [$4/3 gs^2 + 1/9 gy^2 + y_d^2$],
    [$(D, B)$], [$4/3 gs^2 + 1/9 gy^2 - 1/2 y_d^2$], [$(L, W)$], [$3/4 gw^2 + 1/4 gy^2$],
    [$(L, B)$], [$3/4 gw^2 + 1/4 gy^2 + y_e^2$], [$(E, B)$], [$gy^2 + 1/2 y_e^2$],
  )},
  caption: [Coefficients of the finite fermion--gauge rules @cr-fin-ffv.],
) <cr-tab-ffv>

*Source vertices* ($Q^(f r)_R$, $r_alpha = -1/4 C^((r))_"ad"$), for $r = G, W$ and every multiplet $f$ on which $r$ acts. The rule is $r_alpha lf$ times the tree rule @eq:V-src-matter:
$ overline(srcF)_f gh^a Psi_f: quad 1/4 C^((r))_"ad" lf thin Tg^a_f P_(c_f) , quad quad
  gh^a overline(Psi)_f srcF_f: quad 1/4 C^((r))_"ad" lf thin P_(overline(c)_f) Tg^a_f , quad quad 1/4 C^((r))_"ad" = 3/4 gs^2, thin 1/2 gw^2 . $ <cr-fin-src>

*Higgs--gauge vertices*, with $H_j$ entering with $p$ and $H^dagger_i$ leaving with $p'$:
$ Q^((r))_J, thin r = W, B: quad & H^dagger_i H_j A^a_mu: quad -i/3 tilde(Y) lf thin (bar(p) + bar(p)')^mu (Tg^a_H)_(i j) , \
  Q_(W W): quad & H^dagger_i H_j W^I_mu W^J_nu: quad -i/2 tilde(Y) gw^2 lf thin delta^(I J) delta_(i j) gmb^(mu nu) , \
  Q_(W B): quad & H^dagger_i H_j W^I_mu B_nu: quad -i/2 (y_u^2 + 2 y_d^2) gw gy lf thin (tau^I)_(i j) gmb^(mu nu) , \
  Q_(B B): quad & H^dagger_i H_j B_mu B_nu: quad -i/2 (13/3 y_u^2 + 7/3 y_d^2 + 3 y_e^2) gy^2 lf thin delta_(i j) gmb^(mu nu) , \
  Q_(G G): quad & H^dagger_i H_j G^A_mu G^C_nu: quad -i (y_u^2 + y_d^2) gs^2 lf thin delta^(A C) delta_(i j) gmb^(mu nu) . $ <cr-fin-higgs>
The rule of $Q^((r))_J$ is $-tilde(Y)\/3$ times the barred tree rule @eq:V-ssv, because $i [(partial_mu H)^dagger Am^mu H - H^dagger Am^mu partial_mu H]$ has the vertex function $-(bar(p) + bar(p)')^mu Tg^a_H$. The two-gauge rules follow from $V_(Q_(W W)) = {Tg^I_H, Tg^J_H} gmb^(mu nu) = 1/2 gw^2 delta^(I J) gmb^(mu nu)$, ${Tg^I_H, Tg^B_H} = 1/2 gw gy tau^I$, ${Tg^B_H, Tg^B_H} = 1/2 gy^2$ and $V_(Q_(G G)) = 2 gs^2 delta^(A C) gmb^(mu nu)$.

$Delta S^((1))_"ST"$ has no other vertex: no ghost, antighost, $srcA$, $srcC$ or $srcH$ vertex, no fermion kinetic, Yukawa, $H^4$ or $m^2$ vertex, and no vertex with $epsilon^(mu nu alpha beta)$.

#check[Each rule of @cr-fin-2 to @cr-fin-higgs and of @cr-tab-quartic was compared with the vertex function obtained by functional differentiation of the monomial on random plane waves with the plane-wave engine of @sec:bases, at generic couplings; all agree to machine precision. The trace form @cr-fin-4 and the matrix form of the first line of @eq:c-DS agree with the basis form in every combination of gauge factors.]

=== The counterterm vertices as model input <cr-model>

Each counterterm vertex enters a model as a vertex with three attributes. (i) Its loop order, $1$ here: a diagram with $L$ loops and one counterterm vertex contributes at order $L + 1$. (ii) Its type: singular, with the coefficient proportional to $lf thin 1\/ep$ (@cr-tab-inv, @cr-ev-rules), or finite, proportional to $lf$ (@cr-fin-2 to @cr-fin-higgs). The finite vertices are present only because of the non-anticommuting $g5$ and are not fixed by a pole condition. (iii) Its Lorentz structure in $d$ dimensions. The singular rules contain unsplit, barred and hatted metrics and momenta, and the model must keep them apart as in @sec:bmhv: replacing an unsplit structure by its barred part changes the evanescent operators of @cr-ev-rules, and these enter the two-loop finite parts. The finite rules are purely barred. At two loops the ghost-number-1 vertices $bop_d delta S^((1))$ and $bop_d Delta S^((1))_"ST"$ are also needed; they follow from these rules with @eq:dict. @sec:automation lists what the current tools need to represent these vertices.
