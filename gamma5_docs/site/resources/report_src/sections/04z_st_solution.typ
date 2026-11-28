#import "../macros.typ": *

== Solution of the Slavnov--Taylor conditions <sec:st-solution>

This subsection collects the eighteen classes of @sec:g0 and @sec:g1 into the one-loop counterterm action $delta S^((1)) + Delta S^((1))_"ST"$. It first assembles the singular counterterm (@sec:st-dS), then checks the pole identity of @sec:m-oneloop (@sec:st-poles), collects the finite breaking $X_1$ (@sec:st-X), solves $bop Delta S^((1))_"ST" = -X_1$ by hand (@sec:st-solve), states the result (@sec:st-result), compares it with the literature (@sec:st-lit) and tests the restored identity on renormalised Green's functions with complete finite parts (@sec:st-checks).

*Symbols used throughout.* $S^((r)) = S^((r))_"lh" + S^((r))_"rh"$ is the sum over the left-handed and the right-handed multiplets of @eq:SLR, so $S^((G)) = 2 gs^2$, $S^((W)) = 2 gw^2$, $S^((B)) = 10/3 gy^2$. $S^((r))_H$ and $C_H$ are $S^((r))_f$ and $C_f$ of @tab:invariants for $f = H$: $S^((W))_H = gw^2\/2$, $S^((B))_H = gy^2\/2$, $S^((G))_H = 0$, $C_H = 3/4 gw^2 + 1/4 gy^2$. The Yukawa sums are
$ tilde(Y) = 3 y_u^2 + 3 y_d^2 + y_e^2, quad Y_4 = 3 y_u^4 + 3 y_d^4 + y_e^4, \
  cal(Y)_f = ((y_u^2 + y_d^2)/2, thin y_u^2, thin y_d^2, thin y_e^2/2, thin y_e^2) "on" (Q, U, D, L, E), $ <eq:sts-yuk>
and $G_4 = 3 gw^4 + 2 gw^2 gy^2 + gy^4$. The Yukawa numbers $nu_(f r)$ are those of @sec:g1_cFF: $nu_(Q G) = y_u^2 + y_d^2$, $nu_(Q W) = 0$, $nu_(Q B) = 4 y_u^2 - 2 y_d^2$, $nu_(U G) = 2 y_u^2$, $nu_(U B) = y_u^2\/2$, $nu_(D G) = 2 y_d^2$, $nu_(D B) = -y_d^2$, $nu_(L W) = 0$, $nu_(L B) = 2 y_e^2$, $nu_(E B) = y_e^2$.

*Normalisation of coefficients.* The renormalisation constants $delta Z_x$, $delta y_k$, $delta lambda$, $delta m^2$ include the pole factor $lf 1\/ep$. Where a number without this factor is meant, it is written as $16 pi^2 ep thin delta Z_x$. The restoring counterterm is $Delta S^((1))_"ST" = lf sum_alpha r_alpha Q_alpha$; the coefficients $r_alpha$ do not include $lf$.

=== The singular counterterm $delta S^((1))$ <sec:st-dS>

*Assembly.* The ten classes of @sec:g0 give $delta S^((1)) = -Kop Gamma^((1))$ vertex by vertex (@eq:m-Z1). In four dimensions the result is a linear combination of the 97 monomials of @tab:Qbasis, and it lies in the 25-dimensional kernel of #bop (@eq:closed). Twenty-five monomials fix the coefficients of the 25 closed directions, and the other 72 monomial coefficients are then predictions; all 72 hold exactly. The non-zero coefficients define the renormalisation constants. With bare quantities written as
$ A_0 = Z_A^(1\/2) A, quad tilde(srcA)_0 = Z_A^(-1\/2) tilde(srcA), quad gh_0 = Z_gh^(1\/2) gh, quad srcC_0 = Z_gh^(-1\/2) srcC, quad c_(r 0) = mu^ep Z_(c_r) c_r, \
  psi_(f 0) = Z_f^(1\/2) psi_f, quad srcF_(f 0) = Z_f^(-1\/2) srcF_f, quad H_0 = Z_H^(1\/2) H, quad srcH_0 = Z_H^(-1\/2) srcH, \
  y_(k 0) = mu^ep (y_k + delta y_k), quad lambda_0 = mu^(2 ep) (lambda + delta lambda), quad m_0^2 = m^2 + delta m^2 , $ <eq:sts-bare>
with $Z = 1 + delta Z$, one $Z_A$, $Z_gh$, $Z_(c_r)$ per gauge factor and the sterile partners not renormalised, the singular counterterm is
$ delta S^((1)) = & sum_(r = G, W) [ 1/2 delta Z^((r))_A thin cal(O)^((r))_A + delta Z_(c_r) thin cal(O)^((r))_g + 1/2 delta Z^((r))_gh thin cal(O)^((r))_c ] + delta Z^((B))_A thin Q^((B))_T \
  & + sum_f 1/2 delta Z_f thin cal(O)^f_psi + 1/2 delta Z_H thin cal(O)_H + sum_k delta y_k thin Q^+_(y_k) - delta lambda thin Q_lambda - delta m^2 thin Q_m \
  & - lf 1/ep [ sum_r S^((r)) (2/3 tilde(E)^((r))_F + 1/3 E^((r))_A) + tilde(Y) (tilde(E)_H + 2/3 E_H) ] . $ <eq:sts-dS>
The closed directions of @eq:closed are continued to $d$ dimensions as stated in @sec:evanescent: unsplit indices in the bosonic monomials, the barred vertex in the fermion monomials. The graphs without a closed fermion loop give exactly such $d$-dimensional combinations. A closed fermion loop gives the barred invariants $-1/4 bar(F) bar(F)$ and $(bar(D) H)^dagger bar(D) H$, whose difference to the continued form is $tilde(E)^((r))_F$ and $tilde(E)_H$ of @eq:evanescent, with the same coefficient in the two-, three- and four-point functions. The hatted structures of the two-point functions that no invariant carries are $E^((r))_A$ and $E_H$. In the abelian factor $cal(O)^((B))_A$ and $cal(O)^((B))_c$ have zero coefficient; the $B$ field renormalisation sits in $Q^((B))_T$ and in the $B$ couplings of $cal(O)^f_psi$ and $cal(O)_H$. The constants are listed in @tab:sts-Z.

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([constant $x$], [$16 pi^2 ep thin x$, general form], [$16 pi^2 ep thin x$, Standard Model]),
    [$delta Z^((r))_A$], [$5/3 C^((r))_"ad" - 2/3 S^((r)) - 1/3 S^((r))_H$], [$11/3 gs^2$, $11/6 gw^2$, $-43/18 gy^2$],
    [$delta Z_(c_r)$], [$-11/6 C^((r))_"ad" + 1/3 S^((r)) + 1/6 S^((r))_H$], [$-29/6 gs^2$, $-35/12 gw^2$, $43/36 gy^2$],
    [$delta Z^((r))_gh$], [$8/3 C^((r))_"ad" - 2/3 S^((r)) - 1/3 S^((r))_H$], [$20/3 gs^2$, $23/6 gw^2$ ($r = G, W$)],
    [$tilde(Z)^((r))_3 - 1$], [$1/2 C^((r))_"ad"$], [$3/2 gs^2$, $gw^2$],
    [$tilde(Z)^((r))_1 - 1$], [$-1/2 C^((r))_"ad"$], [$-3/2 gs^2$, $-gw^2$],
    [$delta Z_f$], [$-(C_f + cal(Y)_f)$], [$Q$: $-(4/3 gs^2 + 3/4 gw^2 + 1/36 gy^2 + (y_u^2 + y_d^2)\/2)$],
    [], [], [$U$: $-(4/3 gs^2 + 4/9 gy^2 + y_u^2)$],
    [], [], [$D$: $-(4/3 gs^2 + 1/9 gy^2 + y_d^2)$],
    [], [], [$L$: $-(3/4 gw^2 + 1/4 gy^2 + y_e^2\/2)$],
    [], [], [$E$: $-(gy^2 + y_e^2)$],
    [$delta Z_H$], [$2 C_H - tilde(Y)$], [$3/2 gw^2 + 1/2 gy^2 - tilde(Y)$],
    [$(delta y_u) \/ y_u$], [], [$9/4 y_u^2 + 3/4 y_d^2 + 1/2 y_e^2 - 4 gs^2 - 9/8 gw^2 - 17/24 gy^2$],
    [$(delta y_d) \/ y_d$], [], [$3/4 y_u^2 + 9/4 y_d^2 + 1/2 y_e^2 - 4 gs^2 - 9/8 gw^2 - 5/24 gy^2$],
    [$(delta y_e) \/ y_e$], [], [$3/2 y_u^2 + 3/2 y_d^2 + 5/4 y_e^2 - 9/8 gw^2 - 15/8 gy^2$],
    [$delta lambda$], [], [$12 lambda^2 - 6 lambda C_H + 2 lambda tilde(Y) - Y_4 + 3/16 G_4$],
    [$delta m^2$], [], [$m^2 (6 lambda - 3 C_H + tilde(Y))$],
  ),
  caption: [Renormalisation constants of @eq:sts-dS, Feynman gauge, $overline("MS")$. Each constant $x$ includes the factor $lf 1\/ep$; the second and third columns give $16 pi^2 ep thin x$. The three entries of a row are $r = G, W, B$. $tilde(Z)^((r))_3 - 1 = 1/2 (delta Z^((r))_gh - delta Z^((r))_A)$ and $tilde(Z)^((r))_1 - 1 = delta Z_(c_r) + 1/2 delta Z^((r))_gh$ are the conventional ghost-propagator and ghost-vertex constants; $tilde(Z)^((r))_1 - 1$ is the common coefficient of all source counterterms. In the abelian factor $Z_(c_B) (Z^((B))_A)^(1\/2) = 1$; the ghost $gh^B$ is free, and $delta Z^((B))_gh = -2 delta Z_(c_B)$, i.e. $16 pi^2 ep thin delta Z^((B))_gh = -43/18 gy^2$, is a convention without effect.],
) <tab:sts-Z>

*Beta functions.* From $c_(r 0) = mu^ep Z_(c_r) c_r$ and the analogous relations of @eq:sts-bare, at one loop $16 pi^2 beta_x = 2 x thin (16 pi^2 ep thin (delta x) \/ x)$ for $x = c_r, y_k$ (with $delta c_r = c_r thin delta Z_(c_r)$) and $16 pi^2 beta_lambda = 2 (16 pi^2 ep thin delta lambda)$:
$ 16 pi^2 beta_(gs) & = -29/3 gs^3, quad 16 pi^2 beta_(gw) = -35/6 gw^3, quad 16 pi^2 beta_(gy) = 43/18 gy^3, \
  16 pi^2 beta_(y_u) & = y_u (9/2 y_u^2 + 3/2 y_d^2 + y_e^2 - 8 gs^2 - 9/4 gw^2 - 17/12 gy^2), \
  16 pi^2 beta_(y_d) & = y_d (3/2 y_u^2 + 9/2 y_d^2 + y_e^2 - 8 gs^2 - 9/4 gw^2 - 5/12 gy^2), \
  16 pi^2 beta_(y_e) & = y_e (3 y_u^2 + 3 y_d^2 + 5/2 y_e^2 - 9/4 gw^2 - 15/4 gy^2), \
  16 pi^2 beta_lambda & = 24 lambda^2 - 3 lambda (3 gw^2 + gy^2) + 4 lambda tilde(Y) - 2 Y_4 + 3/8 G_4, \
  16 pi^2 thin (dif ln m^2) / (dif ln mu) & = 12 lambda + 2 tilde(Y) - 9/2 gw^2 - 3/2 gy^2 . $ <eq:sts-beta>
These are the known one-loop beta functions of the Standard Model with one family @Machacek:1983tz: $16 pi^2 beta_(c_r) = c_r^3 (-11/3 C^((r))_"ad" + 2/3 S^((r)) + 1/3 S^((r))_H) \/ c_r^2$, i.e. $b_3 = -11 + 4/3$, $b_2 = -22/3 + 4/3 + 1/6$, $b_1 = 20/9 + 1/6$ for $g'$ without GUT normalisation. The four-dimensional part of $delta S^((1))$ therefore coincides with the standard $overline("MS")$ counterterm.

*Evanescent part.* The non-invariant evanescent counterterms of @eq:sts-dS have the Standard-Model values
$ -lf 1/ep [ 2/3 gs^2 E^((G))_A + 2/3 gw^2 E^((W))_A + 10/9 gy^2 E^((B))_A + 2/3 tilde(Y) E_H ] , $
and the $bop_d$-closed evanescent counterterms are $-lf 1/ep [4/3 gs^2 tilde(E)^((G))_F + 4/3 gw^2 tilde(E)^((W))_F + 20/9 gy^2 tilde(E)^((B))_F + tilde(Y) tilde(E)_H]$. The class subsections @sec:g0_VV to @sec:g0_VVVV quote the same coefficient of $tilde(E)^((r))_F$ in their two-, three- and four-point parts.

#check[*Cross-class conditions.* Write $[Q_alpha]$ for the coefficient of $Q_alpha$ (@tab:Qbasis) in $16 pi^2 ep thin delta S^((1))$, bosonic plus fermionic parts. The assembly requires, and the class results satisfy:
- $[cal(Q)^((r))_"YM"] = 2 [Q^((r))_(3 f)] - [Q^((r))_T]$, from @sec:g0_VV, @sec:g0_VVV, @sec:g0_VVVV ($G$: $-7/3 gs^2 = 4/3 gs^2 - 11/3 gs^2$; $W$: $-13/6 gw^2 = -1/3 gw^2 - 11/6 gw^2$);
- $[Q^((r))_(3 f)] - [Q^((r))_T] = -C^((r))_"ad"$, which links the gauge classes to @sec:g0_ghost_sources;
- $[Q^(f r)_V] - [Q^f_K] = -C^((r))_"ad"$ for all ten pairs $(f, r)$, from @sec:g0_FF and @sec:g0_FFV;
- $[Q_rho [fg_((r))]] = [Q_zeta [-1/2 fg_((r))]] = [Q^(f r)_R] = [Q^((W))_sigma] = -C^((r))_"ad"\/2$ and $[Q^((r))_(rho partial)] = C^((r))_"ad"\/2$ (@sec:g0_ghost_sources);
- the five coefficients $[Q^((W))_J]$, $[Q^((B))_J]$, $[Q_(W W)]$, $[Q_(B B)]$, $[Q_(W B)]$ of @sec:g0_SSV_SSVV follow from $delta Z_H$ of @sec:g0_SS and from $16 pi^2 ep (delta Z^((W))_A \/ 2 + delta Z_(c_W)) = -C^((W))_"ad"$;
- in the abelian factor $[Q^((B))_sigma] = [Q^((B))_(rho partial)] = [Q^(f B)_R] = 0$ and $[Q^(f B)_V] = [Q^f_K]$.
Changing any single class entry breaks the fit. The plane-wave engine of @sec:bases gives $bop (delta S^((1))|_(d = 4)) = 0$ on all 39 ghost-number-1 leg classes, and every $Z$ factor agrees with the standard Feynman-gauge values @Machacek:1983tz @Belusca-Maito:2020ala.]

=== The pole identity <sec:st-poles>

*The variations.* In @eq:sts-dS, $cal(O)^((r))_A$, $Q^((B))_T$, $cal(O)^f_psi$, $cal(O)_H$, $Q^+_(y_k)$, $Q_lambda$, $Q_m$, $tilde(E)^((r))_F$ and $tilde(E)_H$ are $bop_d$-closed (@sec:evanescent). The other terms give, with @eq:bd-Oc, @eq:b-evanescent and @eq:Nhat,
$ bop_d cal(O)^((r))_g = bop_d cal(O)^((r))_c = Dbrk^((r))_1, quad bop_d E^((r))_A = hat(N)^((r))_A + hat(N)^((r))_(A A), quad bop_d E_H = hat(N)_H . $
The coefficient of $Dbrk^((r))_1$ is $delta Z_(c_r) + 1/2 delta Z^((r))_gh = tilde(Z)^((r))_1 - 1 = -lf 1/ep thin 1/2 C^((r))_"ad"$ (@tab:sts-Z). Hence
$ bop_d delta S^((1)) = - lf 1/ep [ & sum_(r = G, W) 1/2 C^((r))_"ad" Dbrk^((r))_1 \
  & + sum_(r = G, W, B) 1/3 S^((r)) (hat(N)^((r))_A + hat(N)^((r))_(A A)) + 2/3 tilde(Y) hat(N)_H ] . $ <eq:sts-bd-dS>
It is purely evanescent: the four-dimensional parts of $bop_d delta S^((1))$ cancel in every class, which is the statement $bop (delta S^((1))|_(d=4)) = 0$ above. The identity $[(Dbrk_1 dot Gamma)^((1))]_"div" + bop_d delta S^((1)) = 0$ of @sec:m-oneloop therefore predicts the one-loop poles of the ghost-number-1 classes. Equivalently, $Dbrk^((r))_1$ is the evanescent part of $bop_d Q^(f r)_R$ (@eq:m-bYR): the source counterterm $-lf 1/ep thin 1/2 C^((r))_"ad" Q^(f r)_R$ is the only term that cancels the pole of the $gh overline(Psi) Psi$ graphs with an internal ghost.

#figure(
  {
    set par(justify: false)
    table(
    columns: (auto, auto, 1fr, auto),
    align: (left, left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([class], [legs], [pole of $(Dbrk_1 dot Gamma)^((1))$, in units of $lf 1\/ep$], [SM coefficient]),
    [@sec:g1_cV], [$gh^a A^b_nu (-q)$], [$1/3 S^((r)) hat(N)^((r))_A$, vertex function $i/3 S^((r)) hat(q)^2 bar(q)^nu delta^(a b)$], [$2/3 gs^2, 2/3 gw^2, 10/9 gy^2$],
    [@sec:g1_cVV], [$gh^a A^b_mu (p_1) A^c_nu (p_2)$], [$1/3 S^((r)) hat(N)^((r))_(A A)$, vertex function $-1/3 S^((r)) fg^(a b c) gmb^(mu nu) (hat(p)_1^2 - hat(p)_2^2)$], [$2/3 gs^2, 2/3 gw^2$],
    [@sec:g1_cVVV_cVVVV], [$gh A A A$, $gh A^4$], [$0$], [],
    [@sec:g1_cFF], [$gh^a overline(Psi)_f Psi_f$], [$1/2 C^((r))_"ad" Dbrk^((f r))_1$; no four-dimensional pole], [$3/2 gs^2$, $gw^2$, $0$],
    [@sec:g1_cFFV], [$gh overline(Psi) Psi A$], [$0$], [],
    [@sec:g1_cSS_cSSV], [$gh^a H_j (p_1) H^dagger_i (p_2)$], [$2/3 tilde(Y) hat(N)_H$, vertex function $-(2 i)/3 tilde(Y) (hat(p)_1^2 - hat(p)_2^2) (Tg^a_H)_(i j)$], [$2 y_u^2 + 2 y_d^2 + 2/3 y_e^2$],
    [@sec:g1_cSS_cSSV], [$gh H^dagger H A$, $A A$; $gh (H^dagger H)^2$], [$0$], [],
    [@sec:g1_cFFS], [$gh overline(Psi) Psi H$], [$0$], [],
    [@sec:g1_ccRF], [$gh gh overline(srcF) Psi$, $gh gh overline(Psi) srcF$], [$0$: the $-C_"ad"\/2$ of $Q^(f r)_R$ and of $Q_zeta$ cancel], [],
    [no graph], [$gh gh tilde(srcA)$, $gh gh tilde(srcA) A$, $gh gh gh srcC$, $gh gh srcH^dagger H$, $gh gh H^dagger srcH$], [$0$], [],
  )},
  caption: [The pole identity class by class. Third column: the pole predicted by @eq:sts-bd-dS, equal to the pole found in the class subsection; vertex functions as in @eq:Nhat-V, ghost first, incoming momenta. $r = G, W, B$ in the order of the last column (for $gh A A$ and $gh overline(Psi) Psi$ the $B$ entry is zero). $hat(N)$ and $Dbrk^((f r))_1$ are those of @eq:Nhat and @eq:bd-Oc.],
) <tab:sts-poles>

#check[The identity holds in all 146 leg assignments of the eight classes, evanescent poles included. The check was done in $d = 6$ with explicit components (two hatted directions with non-zero external hatted momenta), $8 times 8$ Dirac matrices, generic couplings and two random kinematic points per assignment, with $bop_d$ built from the $d$-dimensional $S_0$. $sop S_0 = Dbrk_1$, $bop_d cal(O)^((r))_g = bop_d cal(O)^((r))_c = Dbrk^((r))_1$, $bop_d$ of the invariant part $= 0$ and @eq:b-evanescent were verified in the same set-up.]

=== The finite breaking $X_1$ <sec:st-X>

@tab:sts-X collects $X_1$ of @eq:m-X1 from the eight classes of @sec:g1 in the basis of @tab:Nbasis. It has 71 non-zero components. No component lies outside $"span"{N_beta}$, and $X_1$ lies in $bop("span"{Q_alpha})$ (@sec:st-solve).

#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([class], [$16 pi^2 X_1$]),
    [@sec:g1_cV], [$-1/3 S^((r)) N^((r))_3$: $-2/3 gs^2 N^((G))_3 - 2/3 gw^2 N^((W))_3 - 10/9 gy^2 N^((B))_3$],
    [@sec:g1_cVV], [$1/3 S^((r)) (N^((r))_(f 1) - 2 N^((r))_(f 2))$, $r = G, W$; $epsilon$ part @eq:sts-anom, zero in the SM],
    [@sec:g1_cVVV_cVVVV], [$-4/9 gs^2 N^((G))_(f f) - 1/3 gw^2 N^((W))_(f f) - 1/18 N_(G 4 a) - 1/9 N_(G 4 b) - 1/6 N_(W 4 a) - 95/162 N_(B 4)$ \ $- 1/9 N_(G G G B 1) - 1/18 (N_(G G G B 2) + N_(G G G B 3))$ \ $- sum_(r s) k_(r s) (N^((r s))_a + N^((r s))_b + 2 N^((r s))_c + 2 N^((r s))_d)$, \ $k_(G W) = 1/12$, $k_(G B) = 11/108$, $k_(W B) = 1/18$; $gh A A A A$: $0$],
    [@sec:g1_cFF], [$sum_(f, r) x_(f r) N^(f r)_F$, $x_(f r) = -(C_f - 1/4 C^((r))_"ad") - 1/2 nu_(f r)$ (@tab:sts-x)],
    [@sec:g1_cFFV], [$1/4 C^((r))_"ad" thin i (N^(f r r)_(F V) - N^(f r)_(V F))$, $r = G, W$: $3/4 gs^2$ ($Q, U, D$), $1/2 gw^2$ ($Q, L$)],
    [@sec:g1_cSS_cSSV], [$-1/3 tilde(Y) (N^((W))_S + N^((B))_S) + 2/3 tilde(Y) N^((W W))_(S V) + (y_d^2 - 1/3 y_e^2)(N^((W B))_(S V) + N^((B W))_(S V))$ \ $+ (10/3 y_u^2 + 4/3 y_d^2 + 8/3 y_e^2) N^((B B))_(S V) + (y_u^2 + y_d^2) N_(S G)$],
    [@sec:g1_cFFS], [$1/2 gw^2 sum_(k = u, d, e) (-N^(O_k W)_(Y L) + N^(O_k^dagger W)_(Y R))$; gluon and $gh^B$ components $0$],
    [@sec:g1_ccRF], [$1/4 C^((r))_"ad" (N^(f r r)_(R R) + overline(N)^(f r r)_(R R))$, $r = G, W$: $3/4 gs^2$ ($Q, U, D$), $1/2 gw^2$ ($Q, L$)],
    [source classes], [$0$ in $gh gh tilde(srcA)$, $gh gh tilde(srcA) A$, $gh gh gh srcC$, $gh gh srcH^dagger H$, $gh gh H^dagger srcH$],
  ),
  caption: [The one-loop breaking $X_1$ in the basis $N_beta$ of @tab:Nbasis. All components not listed vanish, among them every component along $N^((r))_1$, $S_i [t]$, $N^((r))_(r B i)$, $N_(S 2)$, $N_(S E)$, $N_(S m)$, $N_(S C)$, $N_(S V E i)$, $N_(S lambda)$, $N_cal(J) [t]$ and the mixed $N^(f r s)_(F V)$, $N^(f r s)_(R R)$ with $r != s$.],
) <tab:sts-X>

#figure(
  table(
    columns: (auto, 1fr, auto, 1fr),
    align: (left, left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$(f, r)$], [$x_(f r)$], [$(f, r)$], [$x_(f r)$]),
    [$(Q, G)$], [$-7/12 gs^2 - 3/4 gw^2 - 1/36 gy^2 - (y_u^2 + y_d^2)\/2$], [$(Q, W)$], [$-4/3 gs^2 - 1/4 gw^2 - 1/36 gy^2$],
    [$(Q, B)$], [$-4/3 gs^2 - 3/4 gw^2 - 1/36 gy^2 - 2 y_u^2 + y_d^2$], [$(U, G)$], [$-7/12 gs^2 - 4/9 gy^2 - y_u^2$],
    [$(U, B)$], [$-4/3 gs^2 - 4/9 gy^2 - y_u^2\/4$], [$(D, G)$], [$-7/12 gs^2 - 1/9 gy^2 - y_d^2$],
    [$(D, B)$], [$-4/3 gs^2 - 1/9 gy^2 + y_d^2\/2$], [$(L, W)$], [$-1/4 gw^2 - 1/4 gy^2$],
    [$(L, B)$], [$-3/4 gw^2 - 1/4 gy^2 - y_e^2$], [$(E, B)$], [$-gy^2 - y_e^2\/2$],
  ),
  caption: [Standard-Model values of the coefficients $x_(f r)$ of $N^(f r)_F$ in $16 pi^2 X_1$.],
) <tab:sts-x>

*Wess--Zumino condition.* $bop X_1 = 0$ (@eq:m-Xl) was evaluated with the plane-wave engine of @sec:bases on all 34 ghost-number-2 leg classes ($gh gh$ with gauge, Higgs, fermion and Yukawa legs, $gh gh gh overline(srcF) psi$, $gh gh gh overline(psi) srcF$) at a coupling point not used anywhere else. In the 17 classes where single terms are non-zero ($gh gh A$, $gh gh A A$, $gh gh H^dagger H$, $gh gh overline(psi) psi$, $gh gh overline(psi) psi A$, $gh gh overline(psi) psi H$), $|bop X_1|$ is below $10^(-14)$ times the largest single term; in the other classes every term vanishes.

*Anomaly.* Each multiplet gives a parity-odd part of $X_1$. Summed over the multiplets it is proportional to $D^(a b c)$ of @eq:Dabc:
$ 16 pi^2 X_"anom" = & -1/3 integral D^(a b c) epsilon^(mu nu rho sigma) gh^a partial_mu A^b_nu partial_rho A^c_sigma \
  & + 1/12 integral D^(a e d) fg^(b c e) epsilon^(mu nu rho sigma) (partial_mu gh^a) A^b_nu A^c_rho A^d_sigma , $ <eq:sts-anom>
from @sec:g1_cVV and @sec:g1_cVVV_cVVVV; per Weyl multiplet the first term is $plus.minus 1\/(24 pi^2)$ times $epsilon^(mu nu rho sigma) Tr_f (wm partial_mu Am_nu partial_rho Am_sigma)$, the consistent anomaly. With $D^(A B C) = D_(G G G) gs^3 d^(A B C)$, $D^(B B B) = D_(B B B) gy^3$ and $D^(B a b) = D_(B r r) c_r^2 gy delta^(a b)$ for $a, b$ in the factor $r$,
$ 16 pi^2 X_"anom" = & -1/3 [ D_(G G G) N_(epsilon d) + D_(B B B) N_(epsilon B) + sum_(r = G, W) D_(B r r) (N^((r))_(epsilon B r r) + 2 N^((r))_(epsilon r B)) ] \
  & + 1/12 D_(G G G) N_(epsilon d f) + sum_(r = G, W) 1/12 D_(B r r) (N^((r))_(epsilon r B f) - N^((r))_(epsilon 4 B)) . $ <eq:sts-anom-basis>
$bop X_"anom" = 0$ holds multiplet by multiplet. The only parity-odd monomials of @tab:Qbasis are $Q^((r))_(B epsilon)$ and $Q^((r))_(4 epsilon)$; adding $-sum_r D_(B r r) (2/3 Q^((r))_(B epsilon) + 1/4 Q^((r))_(4 epsilon))$ to $Delta S$ leaves
$ -1/3 D_(G G G) N_(epsilon d) + 1/12 D_(G G G) N_(epsilon d f) - 1/3 D_(B B B) N_(epsilon B) - sum_(r = G, W) D_(B r r) (N^((r))_(epsilon B r r) + 1/3 N^((r))_(epsilon 4 B)) , $
which lies outside $bop("span"{Q_alpha})$ unless the four numbers $D$ vanish: these are the anomaly directions. For @tab:fields, with $N^c_f$ and $N^w_f$ the colour and $"SU"(2)_L$ multiplicities of the multiplet $f$, $Tr_f (Tg^B {Tg^I, Tg^J}) = gy gw^2 Y_f N^c_f delta^(I J)$ for doublets, $Tr_f (Tg^B {Tg^A, Tg^C}) = gy gs^2 Y_f N^w_f delta^(A C)$ for triplets and $Tr_f (Tg^A {Tg^B, Tg^C}) = 1/2 gs^3 N^w_f d^(A B C)$. The left-handed multiplets $Q, L$ enter with $+$, the right-handed $U, D, E$ with $-$:
$ D_(B B B) & = 2 [6 (1/6)^3 + 2 (-1/2)^3 - 3 (2/3)^3 - 3 (-1/3)^3 - (-1)^3] \
  & = 2 [1/36 - 1/4 - 8/9 + 1/9 + 1] = 0, \
  D_(B W W) & = 3 dot 1/6 - 1/2 = 0, quad D_(B G G) = 2 dot 1/6 - (2/3 - 1/3) = 0, quad D_(G G G) = 1/2 (2 - 1 - 1) = 0 . $ <eq:sts-D>
All other components of $D^(a b c)$ vanish because $T^A$ and $tau^I$ are traceless and $Tr(tau^I {tau^J, tau^K}) = 0$. Hence $X_"anom" = 0$ in the Standard Model. The plane-wave engine confirms it: on $gh A A$ and $gh A A A$ the parity-odd parts of the single multiplets are of order $10^2$ and their sum is below $10^(-13)$.

=== Solving $bop Delta S^((1))_"ST" = -X_1$ <sec:st-solve>

The restoring counterterm is written $Delta S^((1))_"ST" = lf sum_alpha r_alpha Q_alpha$ over the 64 monomials of @tab:Qbasis that (N1)--(N4) of @sec:normalisation leave. $16 pi^2 X_1$ is given in @tab:sts-X. The equation is
$ bop sum_alpha r_alpha Q_alpha = - 16 pi^2 X_1 . $ <eq:sts-eq>
It is solved class of external legs by class: in each class the coefficient of every $N_beta$ on the left must equal that on the right. The coefficient of $Q_alpha$ is written $r_alpha$ or $r[Q_alpha]$. Every #bop variation below follows from @eq:brst, @eq:brst-matter, @eq:b-sources and the dictionary @eq:dict; vertex functions follow @sec:dictionary (incoming momenta, ghost leg $gh(u, q)$ first, $partial -> -i p$). All of them agree with @eq:btable-gauge, @eq:btable-fermion and @eq:btable-higgs.

==== The variations used

*Source-free monomials.* For $Y$ without sources, $bop Y = sop Y$. In position space @eq:dict-i replaces one gauge field by $partial gh$, and @eq:dict-ii is the rotation $sop A = -i [Am, wm]$, $sop psi = i wm psi$, $sop overline(psi) = i overline(psi) wm$, $sop H = i wm H$, $sop H^dagger = -i H^dagger wm$ of every leg. For a rigidly invariant $Y$ the rotations cancel at constant $gh$, so @eq:dict-ii keeps only the terms in which a derivative of $Y$ acts on the ghost.

*Source monomials.* For a gauge factor $r$ let $Omega_r$ be the odd derivation (acting from the left) that rotates a matter field by the ghost of the factor $r$:
$ Omega_r psi_f = i wm^((r)) psi_f, quad Omega_r overline(psi)_f = i overline(psi)_f wm^((r)), quad Omega_r H = i wm^((r)) H, quad Omega_r H^dagger = -i H^dagger wm^((r)) , $ <eq:sts-theta>
zero on all other fields and sources; $Omega^((f))_r$ is its restriction to $psi_f, overline(psi)_f$ and $Omega^((H))_r$ to $H, H^dagger$. Then $Q^(f r)_R = integral [overline(srcF)_f Omega_r psi_f - (Omega_r overline(psi)_f) srcF_f]$ and $Q^((r))_sigma = integral [srcH^dagger Omega_r H + (Omega_r H^dagger) srcH]$. In @eq:b-action $delta'_Y Phi = Omega_r Phi$ for the fields of $Y$, so the second term of @eq:b-action is $Omega_r S_0$. Its part from $cal(L)_"ext"$ combines with $sop Y$ into $overline(srcF)_f [sop (Omega_r psi_f) + Omega_r (sop psi_f)]$ and its analogues. With $sop wm^((r)) = i wm^((r)) wm^((r))$ and ${wm^((r)), wm} = 2 wm^((r)) wm^((r))$ (generators of different factors commute, ghosts anticommute),
$ sop (Omega_r psi_f) + Omega_r (sop psi_f) = wm^((r)) wm^((r)) psi_f, quad quad sop (Omega_r overline(psi)_f) + Omega_r (sop overline(psi)_f) = -overline(psi)_f wm^((r)) wm^((r)) , $
and the same with $H$, $H^dagger$. $overline(srcF)_f$ is even and $srcH^dagger$ is odd, which gives the relative sign of
$ bop Q^(f r)_R = Omega^((f))_r S_0 + N^(f r r)_(R R) + overline(N)^(f r r)_(R R), quad quad
  bop Q^((r))_sigma = Omega^((H))_r S_0 - N^((r r))_sigma - overline(N)^((r r))_sigma , $ <eq:sts-bR>
where $S_0$ is its four-dimensional source-free part, $N^(f r r)_(R R) = integral overline(srcF)_f wm^((r)) wm^((r)) psi_f$, $N^((r r))_sigma = integral srcH^dagger wm^((r)) wm^((r)) H$ (@tab:Nbasis), and $N^((r r))_sigma$, $overline(N)^((r r))_sigma$ vanish for $r = B$. From $Omega_r (overline(psi) Y) = (Omega_r overline(psi)) Y - overline(psi) Omega_r Y$,
$ Omega^((f))_r integral overline(psi)_f (i sl(partial) + gamma^mu Am_mu) psi_f & = integral overline(psi)_f gamma^mu (partial_mu wm^((r))) psi_f + i integral overline(psi)_f gamma^mu [wm^((r)), Am^((r))_mu] psi_f \
  & = N^(f r)_F + i (N^(f r r)_(F V) - N^(f r)_(V F)) . $ <eq:sts-thetaF>
In the classes $gh H^dagger H$ and $gh H^dagger H A$ ($Am$ in the doublet representation, all factors)
$ Omega^((H))_r integral (D_mu H)^dagger D^mu H = & N^((r))_S - integral H^dagger {partial_mu wm^((r)), Am^mu} H \
  & - integral [(partial_mu H)^dagger [wm^((r)), Am^mu] H - H^dagger [wm^((r)), Am^mu] partial_mu H] + dots , $ <eq:sts-thetaH>
with the dots the terms with two gauge fields; $Omega^((H))_r$ annihilates $m^2 H^dagger H$ and $lambda (H^dagger H)^2$.

*Ghost-sector monomials.* $bop srcC_a$ (@eq:b-zeta) contains $partial_mu tilde(srcA)^mu_a + i (overline(srcF)_f Tg^a_f psi_f + overline(psi)_f Tg^a_f srcF_f) - i (srcH^dagger Tg^a_H H + H^dagger Tg^a_H srcH)$, and $tilde(srcA)$ is odd. With $fg^(a b c) Tg^a gh^b gh^c = -2 i wm wm$,
$ bop Q_rho [t] & supset - integral t^(a b c) tilde(srcA)^mu_a (partial_mu gh^b) gh^c, \
  bop Q_zeta [-1/2 fg_((r))] & supset integral fg^(a b c)_((r)) tilde(srcA)^mu_a (partial_mu gh^b) gh^c - sum_f (N^(f r r)_(R R) + overline(N)^(f r r)_(R R)) + N^((r r))_sigma + overline(N)^((r r))_sigma . $ <eq:sts-bghost>
The first term of $bop Q_rho$ comes from $sop A^b_mu supset partial_mu gh^b$ passing the odd $tilde(srcA)$. The other terms of $bop Q_rho [t]$ and $bop Q_zeta [t]$ lie in the classes $gh gh tilde(srcA) A$ and $gh gh gh srcC$ and in classes without sources.

==== Step 0: source classes

$X_1 = 0$ in $gh gh tilde(srcA)$, $gh gh tilde(srcA) A$, $gh gh gh srcC$, $gh gh srcH^dagger H$ and $gh gh H^dagger srcH$ (@tab:sts-X). In $gh gh tilde(srcA)$ the allowed monomials are $Q_rho [t]$ for $t = gs d, (B, W^I, W^I), (W^I, B, W^I), (B, G^A, G^A), (G^A, B, G^A)$ and $Q_zeta [-1/2 fg_((r))]$, $r = G, W$; (N3) has removed $Q_rho [fg_((r))]$ and $Q^((r))_(rho partial)$, (N4) the tensors with $gh^B$. By @eq:sts-bghost they produce seven linearly independent structures $integral t^(a b c) tilde(srcA)^mu_a (partial_mu gh^b) gh^c$, so
$ r[Q_rho [t]] = 0 quad (5 "tensors"), quad r[Q_zeta [-1/2 fg_((G))]] = r[Q_zeta [-1/2 fg_((W))]] = 0 . $ <eq:sts-r-ghost>
Without (N3) only the combination of the closed direction $cal(O)^((r))_c$ would be fixed. With @eq:sts-r-ghost, $bop Delta S$ vanishes in $gh gh tilde(srcA) A$ and $gh gh gh srcC$, as $X_1$ does. In $gh gh srcH^dagger H$ and $gh gh H^dagger srcH$ only $Q^((W))_sigma$ and $Q_zeta [-1/2 fg_((W))]$ contribute; by @eq:sts-bR and @eq:sts-bghost their sum is $(r[Q_zeta [-1/2 fg_((W))]] - r[Q^((W))_sigma]) (N^((W W))_sigma + overline(N)^((W W))_sigma)$, which must vanish, so
$ r[Q^((W))_sigma] = 0 , $ <eq:sts-r-sigma>
and $Q^((B))_sigma$ is removed by (N4). The structures $N^((W B))_sigma$, $overline(N)^((W B))_sigma$ receive no contribution.

#block(breakable: false)[
==== (a) Gauge two-point function, $gh A$

The monomials and structures of this class are (@tab:Qbasis, @tab:Nbasis)
$ Q^((r))_L = -1/2 integral (partial dot A^a)^2, quad Q^((r))_M = 1/2 integral A^a dot A^a, \
  N^((r))_3 = integral gh^a partial^2 partial dot A^a, quad N^((r))_1 = integral gh^a partial dot A^a . $
]
By @eq:dict-i, with $gh^a (u, q)$ and $A^b (e, -q)$,
$ V_(bop Q^((r))_L) = i q^2 (q dot e^a) u^a, quad V_(bop Q^((r))_M) = -i (q dot e^a) u^a, quad V_(bop Q^((r))_T) = 0, \
  V_(N^((r))_3) = -i q^2 (q dot e^a) u^a, quad V_(N^((r))_1) = i (q dot e^a) u^a . $
Hence $bop Q^((r))_L = -N^((r))_3$ and $bop Q^((r))_M = -N^((r))_1$ in this class; $Q^((r))_T$ is transverse and removed by (N2). No other allowed monomial has a $gh A$ vertex. With $16 pi^2 X_1|_(gh A) = -1/3 S^((r)) N^((r))_3$, the coefficients of $N^((r))_3$ and $N^((r))_1$ in @eq:sts-eq give $-r[Q^((r))_L] = 1/3 S^((r))$ and $-r[Q^((r))_M] = 0$:
$ r[Q^((r))_L] = -1/3 S^((r)), quad r[Q^((r))_M] = 0 quad (r = G, W, B) . $ <eq:sts-r-a>

==== (b) Gauge three-point function, $gh A A$

*$fg$ part.* The structures are $N^((r))_(f 1) = integral fg^(a b c) gh^a A^b dot partial^2 A^c$ and $N^((r))_(f 2) = integral fg^(a b c) gh^a A^b_mu partial^mu partial dot A^c$; the cubic monomial is $Q^((r))_(3 f) = -integral fg^(a b c) (partial_mu A^a_nu) A^(b mu) A^(c nu)$. With $gh^a (u, q)$, $A^b (e_1, p_1)$, $A^c (e_2, p_2)$,
$ V_(N^((r))_(f 1)) = fg^(a b c) u^a (e_1 dot e_2)(p_1^2 - p_2^2), quad V_(N^((r))_(f 2)) = fg^(a b c) u^a [(p_1 dot e_1)(p_1 dot e_2) - (p_2 dot e_1)(p_2 dot e_2)] . $
$Q^((r))_L$ gives through @eq:dict-ii $-integral (partial dot (fg^(a b c) A^b gh^c)) (partial dot A^a) = -N^((r))_(f 2)$. $Q^((r))_(3 f)$ gives through @eq:dict-i, after one integration by parts per term and with the terms symmetric in the two antisymmetric indices dropped, $-integral fg^(a b c) [(partial_mu A^a_nu)(partial^mu gh^b) A^(c nu) + (partial_mu A^a_nu) A^(b mu) partial^nu gh^c] = N^((r))_(f 1) - N^((r))_(f 2)$. Hence in this class
$ bop sum_alpha r_alpha Q_alpha = r[Q^((r))_(3 f)] N^((r))_(f 1) - (r[Q^((r))_L] + r[Q^((r))_(3 f)]) N^((r))_(f 2) + dots . $
With $16 pi^2 X_1|_(gh A A) = 1/3 S^((r)) (N^((r))_(f 1) - 2 N^((r))_(f 2))$ the coefficients of $N^((r))_(f 1)$ and $N^((r))_(f 2)$ give
$ r[Q^((r))_(3 f)] = -1/3 S^((r)), quad - r[Q^((r))_L] - r[Q^((r))_(3 f)] = 2/3 S^((r)) quad (r = G, W) . $ <eq:sts-r-b>
The second equation is a consistency condition; with @eq:sts-r-a it reads $1/3 S^((r)) + 1/3 S^((r)) = 2/3 S^((r))$.

*Symmetric part.* $X_1 = 0$ along $N_(d i)$, $N_(B i)$, $N^((r))_(B r i)$, $N^((r))_(r B i)$. The allowed monomials with these structures are $Q_(3 d)$, $Q_(3 B)$, $Q^((r))_(B 1, 2, 3)$, with the images of @eq:btable-gauge restricted to this class. Each image contains a monomial that no other contains ($N_(d 1)$, $N_(B 1)$, $N^((r))_(B r 1)$, $N^((r))_(r B 1)$, $N^((r))_(B r 4)$), so
$ r[Q_(3 d)] = r[Q_(3 B)] = r[Q^((r))_(B 1)] = r[Q^((r))_(B 2)] = r[Q^((r))_(B 3)] = 0 . $ <eq:sts-r-sym>

*$epsilon$ part.* $bop Q^((r))_(B epsilon)|_(gh A A) = N^((r))_(epsilon B r r) - N^((r))_(epsilon r B)$, while $16 pi^2 X_1|_epsilon$ is the first line of @eq:sts-anom-basis with the ratio $1 : 2$, which $bop Q^((r))_(B epsilon)$ cannot produce. For the Standard Model all $D = 0$ (@eq:sts-D), so $X_1|_epsilon = 0$ and
$ r[Q^((r))_(B epsilon)] = 0 . $ <eq:sts-r-eps>

==== (c) Gauge four- and five-point functions

*$gh A A A$.* The cubic monomials contribute through @eq:dict-ii only where the derivative acts on the ghost: $bop Q^((r))_(3 f) supset N^((r))_(f f)$, and $Q_(3 d)$, $Q^((r))_(B 2)$, $Q^((r))_(B epsilon)$ (zero coefficients) would give $N_(d f)$, $-N^((r))_(r B f)$, $-N^((r))_(epsilon r B f)$. The quartic monomials contribute through @eq:dict-i, one gauge field replaced by $partial gh$ at each of four positions:
$ bop Q_(G 4 a) = 4 N_(G 4 a), quad bop Q_(G 4 b) = 4 N_(G 4 b), quad bop Q_(W 4 a) = 4 N_(W 4 a), quad bop Q_(B 4) = 4 N_(B 4), \
  bop Q^((r s))_a = 2 (N^((r s))_a + N^((r s))_b), quad bop Q^((r s))_b = 2 (N^((r s))_c + N^((r s))_d), \
  bop Q_(G G G B) = 2 N_(G G G B 1) + N_(G G G B 2) + N_(G G G B 3) . $
The structures are $N_(G 4 a) = gs^4 integral (partial gh^A dot G^A)(G dot G)$, $N_(G 4 b) = gs^4 integral (partial gh^A dot G^B)(G^A dot G^B)$, $N_(W 4 a) = gw^4 integral (partial gh^I dot W^I)(W dot W)$ and $N^((r))_(f f) = integral fg^(e a b) fg^(e c d) (partial_mu gh^a) A^b_nu A^(c mu) A^(d nu)$. The remaining two quartic variations, $bop Q_(G 4 c) = 4 gs^4 integral T_1$ with $T_1 = d^(A B E) d^(C D E) (partial gh^A dot G^B)(G^C dot G^D)$ and $bop Q_(W 4 b) = 4 gw^4 integral (partial gh^I dot W^J)(W^I dot W^J)$, are reduced as follows. With $T_2 = d^(A C E) d^(B D E) (partial gh^A dot G^B)(G^C dot G^D)$, the cyclic identity of @sec:Qbasis gives $T_1 + 2 T_2 = 1/3 [(partial gh^A dot G^A)(G dot G) + 2 (partial gh^A dot G^B)(G^A dot G^B)]$, and $f^(A B E) f^(C D E) = 2/3 (delta^(A C) delta^(B D) - delta^(A D) delta^(B C)) + d^(A C E) d^(B D E) - d^(B C E) d^(A D E)$ gives $gs^(-2) N^((G))_(f f) = integral [2/3 (partial gh^A dot G^A)(G dot G) - 2/3 (partial gh^A dot G^B)(G^A dot G^B) + T_1 - T_2]$. Eliminating $T_2$, $integral T_1 = 2/3 gs^(-2) N^((G))_(f f) - 1/3 gs^(-4) N_(G 4 a) + 2/3 gs^(-4) N_(G 4 b)$. For $"SU"(2)$, $epsilon^(E I J) epsilon^(E K L) = delta^(I K) delta^(J L) - delta^(I L) delta^(J K)$ gives $gw^4 integral (partial gh^I dot W^J)(W^I dot W^J) = N_(W 4 a) - gw^2 N^((W))_(f f)$. Hence
$ bop Q_(G 4 c) = 8/3 gs^2 N^((G))_(f f) - 4/3 N_(G 4 a) + 8/3 N_(G 4 b), quad bop Q_(W 4 b) = 4 N_(W 4 a) - 4 gw^2 N^((W))_(f f) . $
With $r[Q^((r))_(3 f)]$ of @eq:sts-r-b, $r_(G 4 a) = r[Q_(G 4 a)]$ etc., and $16 pi^2 X_1$ of @tab:sts-X, the coefficients of $N^((G))_(f f)$, $N_(G 4 a)$, $N_(G 4 b)$, $N^((W))_(f f)$, $N_(W 4 a)$, $N_(B 4)$ give
$ -2/3 gs^2 + 8/3 gs^2 r_(G 4 c) = 4/9 gs^2, quad 4 r_(G 4 a) - 4/3 r_(G 4 c) = 1/18, quad 4 r_(G 4 b) + 8/3 r_(G 4 c) = 1/9, \
  -2/3 gw^2 - 4 gw^2 r_(W 4 b) = 1/3 gw^2, quad 4 r_(W 4 a) + 4 r_(W 4 b) = 1/6, quad 4 r_(B 4) = 95/162 , $
so that
$ r_(G 4 a) = 11/72, quad r_(G 4 b) = -1/4, quad r_(G 4 c) = 5/12, quad r_(W 4 a) = 7/24, quad r_(W 4 b) = -1/4, quad r_(B 4) = 95/648 . $ <eq:sts-r-quartic>
$N^((r s))_a$ and $N^((r s))_b$ each give $2 r[Q^((r s))_a] = k_(r s)$, $N^((r s))_c$ and $N^((r s))_d$ each give $2 r[Q^((r s))_b] = 2 k_(r s)$, and $N_(G G G B 1)$ gives $2 r[Q_(G G G B)] = 1/9$, $N_(G G G B 2, 3)$ each give $r[Q_(G G G B)] = 1/18$:
$ r[Q^((G W))_a] = 1/24, quad r[Q^((G W))_b] = 1/12, quad r[Q^((G B))_a] = 11/216, quad r[Q^((G B))_b] = 11/108, \
  r[Q^((W B))_a] = 1/36, quad r[Q^((W B))_b] = 1/18, quad r[Q_(G G G B)] = 1/18 . $ <eq:sts-r-mixed>
The duplicate equations are consistency conditions and hold: equal coefficients of $N^((r s))_a$, $N^((r s))_b$ and of $N^((r s))_c$, $N^((r s))_d$ for the three pairs, and the ratio $2 : 1 : 1$ of $N_(G G G B 1, 2, 3)$. The $epsilon$ structures of $X_1|_(gh A A A)$ are the second line of @eq:sts-anom-basis, zero in the Standard Model, and $bop Q^((r))_(4 epsilon) = 3 N^((r))_(epsilon r B f) + N^((r))_(epsilon 4 B)$ gives $r[Q^((r))_(4 epsilon)] = 0$. $X_1$ and $bop Delta S$ both vanish along $N^((r))_(B f 3)$, $N^((r))_(r B f)$, $N_(d f)$ and $N_cal(J) [t]$.

*$gh A A A A$.* No monomial has five gauge fields and the quartic monomials have no derivative, so @eq:dict-i and @eq:dict-ii give nothing: $bop Delta S|_(gh A A A A) = 0 = -X_1|_(gh A A A A)$ identically.

==== (d) Fermions

*$gh overline(psi) psi$.* The monomials are $Q^(f r)_V = integral overline(psi)_f gamma^mu Am^((r))_mu psi_f$ and $Q^(f r)_R$ of @eq:sts-bR; the structure is $N^(f r)_F = integral overline(psi)_f gamma^mu (partial_mu wm^((r))) psi_f$. $Q^(f r)_V$ gives through @eq:dict-i $-integral overline(psi)_f gamma^mu (partial_mu wm^((r))) psi_f = -N^(f r)_F$; the sign is that of $gh$ passing $overline(psi)_f$, and the rotations @eq:dict-ii cancel. $Q^(f r)_R$ gives $+N^(f r)_F$ by @eq:sts-thetaF. The vertex function, $gh^a (q)$ first, is
$ V_(N^(f r)_F)(gh^a (q), overline(Psi)_f, Psi_f) = i slb(q) thin Tg^a_f P_(c_f) . $
$Q^f_K$ (removed by (N1)) and $Q^((r))_(rho partial)$ (removed by (N3)) would add further multiples of $N^(f r)_F$. Hence $bop Delta S|_(gh overline(psi) psi) = lf sum_(f, r) (r[Q^(f r)_R] - r[Q^(f r)_V]) N^(f r)_F$. With $16 pi^2 X_1|_(gh overline(psi) psi) = sum_(f, r) x_(f r) N^(f r)_F$,
$ r[Q^(f r)_R] - r[Q^(f r)_V] = -x_(f r) . $ <eq:sts-d1>

*$gh gh overline(srcF) psi$.* By @eq:sts-bR and @eq:sts-bghost, $bop Delta S|_(gh gh overline(srcF) psi) = lf sum_(f, r) (r[Q^(f r)_R] - r[Q_zeta [-1/2 fg_((r))]]) (N^(f r r)_(R R) + overline(N)^(f r r)_(R R))$. With $16 pi^2 X_1|_(gh gh overline(srcF) psi) = sum_(f, r = G, W) 1/4 C^((r))_"ad" (N^(f r r)_(R R) + overline(N)^(f r r)_(R R))$ and @eq:sts-r-ghost,
$ r[Q^(f r)_R] = -1/4 C^((r))_"ad" quad (r = G, W), quad quad r[Q^(f B)_R] = 0 quad "(N4)" . $ <eq:sts-r-R>
The mixed structures $N^(f r s)_(R R)$, $r != s$, vanish on both sides. With @eq:sts-d1 and $x_(f r)$ of @tab:sts-X, $r[Q^(f r)_V] = r[Q^(f r)_R] + x_(f r) = -1/4 C^((r))_"ad" - (C_f - 1/4 C^((r))_"ad") - 1/2 nu_(f r)$:
$ r[Q^(f r)_V] = -C_f - 1/2 nu_(f r) quad "for all" f "and" r in f ; $ <eq:sts-r-V>
the term $C^((r))_"ad"\/4$ of $x_(f r)$ cancels against $r[Q^(f r)_R]$. For $r = B$ the two cancelling terms are zero: $C^((B))_"ad" = 0$, and $r[Q^(f B)_R] = 0$ by (N4).

*$gh overline(psi) psi A$ (consistency).* $Q^(f r)_V$ has no derivative and gives nothing. $Q^(f r)_R$ gives $i (N^(f r r)_(F V) - N^(f r)_(V F))$ by @eq:sts-thetaF, with vertex function $V(gh^a, overline(Psi)_f, Psi_f, A^b_mu) = fg^(a b c) Tg^c_f gbar^mu P_(c_f)$. The class result $16 pi^2 X_1|_(gh overline(psi) psi A) = 1/4 C^((r))_"ad" thin i (N^(f r r)_(F V) - N^(f r)_(V F))$ then requires $r[Q^(f r)_R] = -1/4 C^((r))_"ad"$, which is @eq:sts-r-R. $N^(f r s)_(F V)$ with $r != s$ and $N_cal(J) [t]$ vanish on both sides.

==== (e) Higgs

*$gh H^dagger H$.* The monomials are $Q^((r))_J = i integral [(partial_mu H)^dagger Am^((r) mu) H - H^dagger Am^((r) mu) partial_mu H]$ and $Q^((r))_D = integral H^dagger (partial dot Am^((r))) H$; the structures are $N^((r))_S = i integral [(partial_mu H)^dagger (partial^mu wm^((r))) H - H^dagger (partial^mu wm^((r))) partial_mu H]$ and $N^((r))_(S 2) = integral H^dagger (partial^2 wm^((r))) H$. By @eq:dict-i and @eq:sts-thetaH, in this class
$ bop Q^((r))_J = N^((r))_S, quad bop Q^((r))_D = N^((r))_(S 2), quad bop Q^((r))_sigma supset N^((r))_S, \
  V_(N^((r))_S)(gh^a (q), H_j (p_1), H^dagger_i (p_2)) = -i (p_1^2 - p_2^2) (Tg^a_H)_(i j) . $
$Q_(S S)$ (removed by (N1)) would add $N^((W))_S + N^((B))_S$. With $16 pi^2 X_1|_(gh H^dagger H) = -1/3 tilde(Y) (N^((W))_S + N^((B))_S)$, @eq:sts-r-sigma and (N4),
$ r[Q^((W))_J] = r[Q^((B))_J] = 1/3 tilde(Y), quad r[Q^((W))_D] = r[Q^((B))_D] = 0 . $ <eq:sts-r-J>
$X_1$ and $bop Delta S$ vanish along $N^((r))_(S E)$ and $N^((r))_(S m)$.

*$gh H^dagger H A$.* The structures are $N^((r s))_(S V) = integral H^dagger {partial_mu wm^((r)), Am^((s) mu)} H$, $N_(S G) = gs^2 integral H^dagger H thin partial gh^A dot G^A$, and $N_(S C)$, $N_(S V E i)$ of @tab:Nbasis. By @eq:dict-ii, @eq:dict-i and (for $Q^((W))_sigma$) @eq:sts-thetaH with one integration by parts, in this class
$ bop Q^((r))_J = N^((W r))_(S V) + N^((B r))_(S V), quad bop Q^((W))_D = N_(S C), quad bop Q_(W W) = N^((W W))_(S V), quad bop Q_(B B) = N^((B B))_(S V), \
  bop Q_(W B) = N^((W B))_(S V) + N^((B W))_(S V), quad bop Q_(G G) = 2 N_(S G), \
  bop Q^((W))_sigma = -N^((W B))_(S V) - i (2 N^((W W))_(S V E 1) + N^((W W))_(S V E 2) + 2 N^((W W))_(S V E 3) + N^((W W))_(S V E 4)) . $
Only $Q^((W))_sigma$ produces $N_(S V E)$ and only $Q^((W))_D$ produces $N_(S C)$; $X_1$ has neither, which confirms $r[Q^((W))_sigma] = r[Q^((W))_D] = 0$. With $16 pi^2 X_1|_(gh H^dagger H A)$ of @tab:sts-X, the coefficients of the five structures give
$ N^((W W))_(S V)&: quad r[Q^((W))_J] + r[Q_(W W)] = -2/3 tilde(Y), \
  N^((B W))_(S V)&: quad r[Q^((W))_J] + r[Q_(W B)] = -(y_d^2 - y_e^2\/3), \
  N^((W B))_(S V)&: quad r[Q^((B))_J] + r[Q_(W B)] = -(y_d^2 - y_e^2\/3), \
  N^((B B))_(S V)&: quad r[Q^((B))_J] + r[Q_(B B)] = -(10/3 y_u^2 + 4/3 y_d^2 + 8/3 y_e^2), \
  N_(S G)&: quad 2 r[Q_(G G)] = -(y_u^2 + y_d^2) . $
The second and third equations agree because $r[Q^((W))_J] = r[Q^((B))_J]$ (consistency condition). With $r[Q^((r))_J] = 1/3 tilde(Y) = y_u^2 + y_d^2 + 1/3 y_e^2$ of @eq:sts-r-J the solution is
$ r[Q_(W W)] = -tilde(Y), quad r[Q_(W B)] = -(y_u^2 + 2 y_d^2), \
  r[Q_(B B)] = -(13/3 y_u^2 + 7/3 y_d^2 + 3 y_e^2), quad r[Q_(G G)] = -1/2 (y_u^2 + y_d^2) . $ <eq:sts-r-higgs>
In $gh H^dagger H A A$ and $gh (H^dagger H)^2$ the only allowed contributions are the rotations @eq:dict-ii of derivative-free rigidly invariant monomials, which vanish, and $Q^((W))_sigma$ with $r[Q^((W))_sigma] = 0$. Hence $bop Delta S = 0 = -X_1$ there.

==== (f) Yukawa, $gh overline(psi) psi H$

$Q^plus.minus_(y_k)$ are removed by (N1) and closed. The other contributions are $Omega^((f))_r S_"Yuk"$ from $Q^(f r)_R$ and $Omega^((H))_r S_"Yuk"$ from $Q^((r))_sigma$, with $S_"Yuk" = -sum_k y_k integral (O_k + O_k^dagger)$ and $O_k$ of @eq:Ok. Each $O_k$ is invariant factor by factor, $Omega_r O_k = 0$ when $Omega_r$ acts on all three legs. For $O_u$ and $W$, $Omega^((H))_W O_u = -Omega^((Q))_W O_u$ ($U$ is an $"SU"(2)_L$ singlet); for $O_u$ and $G$, $Omega^((U))_G O_u = -Omega^((Q))_G O_u$. From $Omega_r (overline(psi) M psi) = (Omega_r overline(psi)) M psi - overline(psi) M Omega_r psi$ and @tab:Nbasis, $Omega^((Q))_W (-y_u (O_u + O_u^dagger)) = -N^(O_u W)_(Y L) + N^(O_u^dagger W)_(Y R)$. Collecting, with $f_u = f_d = Q$, $f_e = L$, $q_u = U$, $q_d = D$, and without the $gh^B$ terms (N4):
$ bop sum_alpha r_alpha Q_alpha |_(gh overline(psi) psi H) = & sum_k (r[Q^(f_k W)_R] - r[Q^((W))_sigma]) (-N^(O_k W)_(Y L) + N^(O_k^dagger W)_(Y R)) \
  & + sum_(k = u, d) (r[Q^(Q G)_R] - r[Q^(q_k G)_R]) (-N^(O_k G)_(Y L) + N^(O_k^dagger G)_(Y L)) . $
With $16 pi^2 X_1$ of @tab:sts-X the conditions are
$ r[Q^(Q W)_R] - r[Q^((W))_sigma] = r[Q^(L W)_R] - r[Q^((W))_sigma] = -1/2 gw^2, quad r[Q^(Q G)_R] = r[Q^(U G)_R] = r[Q^(D G)_R] ; $
they hold by @eq:sts-r-sigma and @eq:sts-r-R, since $C^((W))_"ad"\/4 = gw^2\/2$. The value $-1/4 C^((r))_"ad"$ of $r[Q^(f r)_R]$ is thus fixed by $gh gh overline(srcF) psi$ and confirmed independently by $gh overline(psi) psi A$ and $gh overline(psi) psi H$.

==== Normalisation and existence

(N1) removes $Q^f_K$, $Q_(S S)$, $Q_m$, $Q_lambda$ and $Q^plus.minus_(y_k)$: it is used in (d) ($Q^f_K$ would shift $r[Q^(f r)_V]$), (e) ($Q_(S S)$ would shift $r[Q^((r))_J]$) and (f). (N2) removes $Q^((r))_T$, which would shift $r[Q^((r))_(3 f)]$ in (b) and all quartic coefficients in (c). (N3) removes $Q^((r))_(rho partial)$ and $Q_rho [fg_((r))]$, used in step 0 and in (d). (N4) removes $Q^(f B)_R$, $Q^((B))_sigma$ and the ghost-sector monomials with $gh^B$, used in step 0, (d), (e) and (f). With these conditions #bop has rank 64 on the 64 remaining monomials, so the solution is unique; on all 97 monomials the rank is 72. A solution exists only if $X_1$ satisfies the eight conditions of @sec:normalisation. Five are $X_1 = 0$ in the source classes (step 0). The other three require that the component of $X_1$ along $i y_k gy integral gh^B (O_k - O_k^dagger)$ vanish for $k = u, d, e$. @sec:g1_cFFS finds no $gh^B$ component, as expected from the absence of an internal $gh^B$ line, so (N4) applies and (N4$'$) is not needed; (N4$'$) gives the identical solution. Every consistency condition above holds: 38 in total.

#check[$bop Delta S^((1))_"ST" + X_1 = 0$ was evaluated with the plane-wave engine of @sec:bases on all 39 ghost-number-1 leg classes, at a coupling point $(gs, gw, gy, y_u, y_d, y_e, lambda, m^2) = (0.97, 0.61, 1.07, 0.53, 0.67, 0.79, 0.41, 0.59)$ and kinematics not used in any fit, with $r_alpha$ from this hand solution and $X_1$ from the class results: the residual is below $2 dot 10^(-15)$ relative to the largest single term (11 classes vanish identically). A numerical solution of the full linear system at eight coupling points, fitted to exact rationals, reproduces every $r_alpha$; an independent implementation of the engine (@sec:supporting) gives the same residuals. Without the $gh gh overline(srcF) psi$ component of $X_1$ (@sec:g1_ccRF) the system has no solution.]

=== Result <sec:st-result>

#block(breakable: false)[#result[*Restoring counterterm.* With (N1)--(N4) the unique solution of $bop Delta S^((1))_"ST" = -X_1$ is, in matrix-valued notation ($Am$ and $wm$ in the representation of the multiplet at hand, all fields, indices and derivatives four-dimensional),
$ Delta S^((1))_"ST" = lf integral dif^4 x { & sum_f Tr_f [ 1/6 (partial dot Am)^2 - i/3 (partial_mu Am_nu) [Am^mu, Am^nu] + 1/3 (Am dot Am)^2 - 1/4 Am_mu Am_nu Am^mu Am^nu ] \
  & - sum_f overline(psi)_f gamma^mu [C_f Am_mu + 1/2 sum_(r in f) nu_(f r) Am^((r))_mu] psi_f
    - 1/4 sum_f [ i overline(srcF)_f wm_C psi_f - i overline(psi)_f wm_C srcF_f ] \
  & + tilde(Y)/3 thin i [(partial_mu H)^dagger Am^mu H - H^dagger Am^mu partial_mu H] - tilde(Y) H^dagger Am^((W)) dot Am^((W)) H \
  & - (y_u^2 + 2 y_d^2) H^dagger {Am^((W))_mu, Am^((B) mu)} H - (13/3 y_u^2 + 7/3 y_d^2 + 3 y_e^2) H^dagger Am^((B)) dot Am^((B)) H \
  & - 1/2 (y_u^2 + y_d^2) gs^2 H^dagger H thin G^A dot G^A } , $ <eq:sts-result>
with $wm_C = sum_(r = G, W) C^((r))_"ad" wm^((r))$ ($C^((G))_"ad" = 3 gs^2$, $C^((W))_"ad" = 2 gw^2$). In the third line $Am = Am^((W)) + Am^((B))$ on the doublet. The first line is $sum_r (r[Q^((r))_L] Q^((r))_L + r[Q^((r))_(3 f)] Q^((r))_(3 f))$ plus the quartic monomials of @eq:sts-r-quartic and @eq:sts-r-mixed, by $sum_f Tr_f (Tg^a_f Tg^b_f) = S^((r)) delta^(a b)$ and $sum_f Tr_f (Tg^a_f [Tg^b_f, Tg^c_f]) = i S^((r)) fg^(a b c)$. The only term with a ghost is the fermion-source term $-1/4 sum_f [dots]$ of the second line. There is no term with the ghost-sector sources $tilde(srcA)$ or $srcC$, no Higgs-source ($srcH$) term, no $gh^B$, no Yukawa-vertex, $H^4$ or $m^2$ term and no $epsilon$ tensor.]]

#figure(
  table(
    columns: (auto, 1fr, auto, 1fr),
    align: (left, left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$Q_alpha$], [$r_alpha$], [$Q_alpha$], [$r_alpha$]),
    [$Q^((G))_L, Q^((W))_L, Q^((B))_L$], [$-2/3 gs^2, -2/3 gw^2, -10/9 gy^2$], [$Q^((G))_(3 f), Q^((W))_(3 f)$], [$-2/3 gs^2, -2/3 gw^2$],
    [$Q_(G 4 a), Q_(G 4 b), Q_(G 4 c)$], [$11/72, -1/4, 5/12$], [$Q_(W 4 a), Q_(W 4 b), Q_(B 4)$], [$7/24, -1/4, 95/648$],
    [$Q^((G W))_a, Q^((G W))_b$], [$1/24, 1/12$], [$Q^((G B))_a, Q^((G B))_b$], [$11/216, 11/108$],
    [$Q^((W B))_a, Q^((W B))_b$], [$1/36, 1/18$], [$Q_(G G G B)$], [$1/18$],
    [$Q^(Q G)_V$], [$-C_Q - (y_u^2 + y_d^2)\/2$], [$Q^(Q W)_V$], [$-C_Q$],
    [$Q^(Q B)_V$], [$-C_Q - 2 y_u^2 + y_d^2$], [$Q^(U G)_V, Q^(U B)_V$], [$-C_U - y_u^2, thin -C_U - y_u^2\/4$],
    [$Q^(D G)_V, Q^(D B)_V$], [$-C_D - y_d^2, thin -C_D + y_d^2\/2$], [$Q^(L W)_V, Q^(L B)_V$], [$-C_L, thin -C_L - y_e^2$],
    [$Q^(E B)_V$], [$-C_E - y_e^2\/2$], [$Q^(Q G)_R, Q^(U G)_R, Q^(D G)_R$], [$-3/4 gs^2$],
    [$Q^(Q W)_R, Q^(L W)_R$], [$-1/2 gw^2$], [$Q^((W))_J, Q^((B))_J$], [$y_u^2 + y_d^2 + y_e^2\/3$],
    [$Q_(W W)$], [$-3 y_u^2 - 3 y_d^2 - y_e^2$], [$Q_(W B)$], [$-y_u^2 - 2 y_d^2$],
    [$Q_(B B)$], [$-13/3 y_u^2 - 7/3 y_d^2 - 3 y_e^2$], [$Q_(G G)$], [$-(y_u^2 + y_d^2)\/2$],
  ),
  caption: [The 39 non-zero coefficients of $Delta S^((1))_"ST" = lf sum_alpha r_alpha Q_alpha$ under (N1)--(N4), with $Q_alpha$ of @tab:Qbasis and $C_f$ of @tab:invariants; the $r_alpha$ do not include $lf$. Zero: $Q^((r))_M$, $Q_(3 d)$, $Q_(3 B)$, $Q^((r))_(B 1,2,3)$, $Q^((r))_(B epsilon)$, $Q^((r))_(4 epsilon)$, $Q^((r))_D$, $Q^((W))_sigma$, the five allowed $Q_rho [t]$ and $Q_zeta [-1/2 fg_((r))]$; the 33 monomials removed by (N1)--(N4) are zero by definition.],
) <tab:sts-r>

=== Comparison with the literature <sec:st-lit>

*Representative.* @Belusca-Maito:2020ala and @Lazopoulos:2026two fix the closed directions differently. Adding to @eq:sts-result the $bop$-closed combination
$ Delta S^((1))_"lit" - Delta S^((1))_"ST" = lf integral dif^4 x [ & 1/2 sum_f Tr_f (-1/4 bar(cal(F))_(mu nu) bar(cal(F))^(mu nu)) \
  & + sum_f C_f overline(psi)_f thin i sl(D) psi_f + tilde(Y)/3 (D_mu H)^dagger D^mu H ] , $ <eq:sts-lit-shift>
with $bar(cal(F))_(mu nu) = partial_mu Am_nu - partial_nu Am_mu - i [Am_mu, Am_nu]$, gives
$ Delta S^((1))_"lit" = lf integral dif^4 x { & sum_f Tr_f [ 1/4 Am_mu partial^2 Am^mu - 5/12 Am_mu partial^mu partial^nu Am_nu + i/6 (partial_mu Am_nu) [Am^mu, Am^nu] + 1/12 (Am dot Am)^2 ] \
  & + sum_f C_f overline(psi)_f thin i sl(partial) psi_f - 1/2 sum_f sum_(r in f) nu_(f r) overline(psi)_f gamma^mu Am^((r))_mu psi_f \
  & - 1/4 sum_f [ i overline(srcF)_f wm_C psi_f - i overline(psi)_f wm_C srcF_f ] \
  & + tilde(Y)/3 (partial_mu H)^dagger partial^mu H - 2/3 tilde(Y) H^dagger Am^((W)) dot Am^((W)) H - (y_d^2 - y_e^2\/3) H^dagger {Am^((W))_mu, Am^((B) mu)} H \
  & - (10/3 y_u^2 + 4/3 y_d^2 + 8/3 y_e^2) H^dagger Am^((B)) dot Am^((B)) H - 1/2 (y_u^2 + y_d^2) gs^2 H^dagger H thin G^A dot G^A } . $ <eq:sts-lit>
In @tab:Qbasis coefficients the shift @eq:sts-lit-shift changes $r[Q^((r))_T]$ to $1/2 S^((r))$ ($gs^2$, $gw^2$, $5/3 gy^2$), $r[Q^((r))_(3 f)]$ to $-1/3 S^((r)) + 1/2 S^((r)) = 1/6 S^((r))$, $(r_(G 4 a), r_(G 4 b), r_(G 4 c))$ to $(1/36, 0, 1/24)$, $(r_(W 4 a), r_(W 4 b))$ to $(1/24, 0)$, $r[Q^f_K]$ to $C_f$, $r[Q^(f r)_V]$ to $-1/2 nu_(f r)$, and gives $r[Q_(S S)] = tilde(Y)\/3$, $r[Q^((r))_J] = 0$, $r[Q_(W W)] = -2/3 tilde(Y)$, $r[Q_(W B)] = -(y_d^2 - y_e^2\/3)$, $r[Q_(B B)] = -(10/3 y_u^2 + 4/3 y_d^2 + 8/3 y_e^2)$. All other coefficients are those of @tab:sts-r.

*Comparison with eq. (6.46) of @Belusca-Maito:2020ala at $xi = 1$.* That paper treats a simple group with one right-handed Weyl multiplet $R$, real scalars and Majorana-form Yukawa couplings. Our Standard Model is of that form with the reducible left-handed representation $Q ⊕ L ⊕ overline(U) ⊕ overline(D) ⊕ overline(E)$, after charge conjugation of the right-handed multiplets. Under this map $g^2 S_2(R) delta^(a b) -> sum_f Tr_f (Tg^a Tg^b)$, $g^2 C_2(R) -> C_f$, $g^2 C_2(G) -> C^((r))_"ad"$ per factor, and the trace $Tr_R (T^a T^b T^c T^d)$ is reversal invariant, so left- and right-handed multiplets enter alike. Their fermion source $R$ is minus ours. Term by term:
- $g^2 S_2(R)\/6 thin (5 S_(G G) + S_(G G G) - integral G^mu partial^2 G_mu)$ is the sum of the first three terms of the first line of @eq:sts-lit, and $g^4 (T_R)^(a b c d)\/12 integral G^a dot G^b G^c dot G^d$ is the fourth;
- $g^2 C_2(R) S_(overline(psi) psi)$ is $sum_f C_f overline(psi)_f i sl(partial) psi_f$;
- the Yukawa term $-1/2 ((Y^m)^* T_(overline(R))^a Y^m) integral g overline(psi) sl(G)^a psi$ is $-1/2 sum_r nu_(f r) overline(psi)_f gamma^mu Am^((r))_mu psi_f$: with the Yukawa matrix $bold(Y)(H)$ from the right-handed to the left-handed multiplets, $sum_phi bold(Y)_phi Am_R bold(Y)_phi^dagger$ on $Q, L$ and $sum_phi bold(Y)^dagger_phi Am_L bold(Y)_phi$ on $U, D, E$ ($bold(Y)_phi = partial bold(Y) \/ partial phi$, $phi in {H_i, H_i^*}$) equal $sum_r nu_(f r) Am^((r))$;
- $-g^2 C_2(G)\/4 (S_(overline(R) c psi) + S_(R c overline(psi)))$ is $-1/4 sum_f [i overline(srcF)_f wm_C psi_f - i overline(psi)_f wm_C srcF_f]$, i.e. $r[Q^(f r)_R] = -1/4 C^((r))_"ad"$;
- $Y_2(S)\/3 thin overline(S_(Phi Phi))$ is $tilde(Y)\/3 (partial H)^dagger partial H$, since $Tr_"lh" (bold(Y)(X) bold(Y)(X)^dagger) = tilde(Y) X^dagger X$;
- the scalar--gauge term $-(cal(C)_R)^(a b)_(m n)\/3 integral g^2\/2 thin G^a dot G^b Phi^m Phi^n$ becomes
  $ -2/3 [Tr_"lh" (Am dot Am thin bold(Y) bold(Y)^dagger) + Tr_"rh" (Am dot Am thin bold(Y)^dagger bold(Y))] + 1/3 Tr_"lh" (Am^"lh"_mu bold(Y) Am^("rh" mu) bold(Y)^dagger) , $
  with $Tr_"lh"$ ($Tr_"rh"$) the trace over $Q, L$ ($U, D, E$) and $Am^"lh"$ ($Am^"rh"$) the gauge matrix on these multiplets. This equals the $H^dagger Am Am H$ terms of the last two lines of @eq:sts-lit for every $H$ and every constant gauge field.
The reference contains in addition the sterile kinetic term $g^2 C_2(R) integral overline(psi) i slb(partial) P_upright(L) psi$, which is $bop$-closed ($sop chi = 0$) and outside our basis. Every term agrees. Apart from the fermion-source term, neither result has a ghost term, and neither has a $srcC$, $srcH$, Yukawa-vertex or $H^4$ finite counterterm. The parity-odd part of $X_1$ (@eq:sts-anom) is cubic in the generators, as the triangle requires, with the magnitude $1\/(24 pi^2)$ of the consistent anomaly per Weyl multiplet; the additional factor $S_2(R)$ in front of $d^(a b c)_R$ in eq. (6.47) of @Belusca-Maito:2020ala would make it quintic and is a misprint. It vanishes for the Standard Model.

*Comparison with @Lazopoulos:2026two.* That paper treats a simple group with fermions only (no scalars), so the comparison is restricted to one simple factor with coupling inside $Am$ and one chiral multiplet. Its $S^((1))_"fct"$ at $xi = 1$ is $Delta S^((1))_"ST" + lf integral [1/3 Tr(-1/4 bar(cal(F)) bar(cal(F))) + C_f overline(psi) i sl(D) psi]$: the quadratic terms $1/6 Tr(Am partial^2 Am) - 1/3 Tr(Am partial partial Am)$, no cubic term, the quartic terms $1/6 Tr (Am dot Am)^2 - 1/12 Tr(Am_mu Am_nu Am^mu Am^nu)$ (its $T^(a b c d)$ and $f f$ terms), $C_f overline(psi) i sl(partial) psi$ (its $(5 + xi)\/6 thin C_2$ at $xi = 1$) and $-1/4 C_"ad" Q_R$. Its $tilde(S)^((1))_"fct"$ differs from $S^((1))_"fct"$ by $-1/24 C(R) integral bar(G)^a_(mu nu) bar(G)^(a mu nu)$ in its notation, i.e. by $1/6 Tr(-1/4 bar(cal(F)) bar(cal(F)))$, and is the gauge--fermion part of @eq:sts-lit. Both agree with our result modulo the closed combinations shown.

*Representative-independent checks.* (i) The longitudinal gauge two-point term: writing the bilinear part as $sum_f [a thin Tr_f (Am partial^2 Am) + b thin Tr_f (Am_mu partial^mu partial^nu Am_nu)]$, $a + b$ does not change under the addition of $Tr bar(cal(F))^2$ or of equations of motion and is fixed by $X_1|_(gh A)$. Here $a + b = -1/6$ ($0 - 1/6$ in @eq:sts-result, $1/4 - 5/12$ in @Belusca-Maito:2020ala, $1/6 - 1/3$ in @Lazopoulos:2026two), i.e. modulo invariants $lf [gs^2\/3 (partial dot G^A)^2 + gw^2\/3 (partial dot W^I)^2 + 5/9 gy^2 (partial dot B)^2]$ for one family. With three families this is $gs^2$, $gw^2$ and $5/3 gy^2$, as in @Fuentes-Martin:2025meh; the $W$ term agrees with @Kuhler:2026sm. (ii) The abelian part of the fermion counterterm: the coefficient of $overline(psi)_f gamma^mu B_mu psi_f$ without Yukawa terms is $-C_f gy Y_f$, which agrees with the general formula of @Cornella:2022hkc restricted to an abelian gauge field. The $C_"ad"$ parts are not comparable there (background-field gauge, no BRST sources).

*Diagnosis.* No genuine difference was found. All differences between representatives are the three $bop$-closed terms of @eq:sts-lit-shift and the sterile kinetic term of @Belusca-Maito:2020ala. Since $bop Delta S_"lit" = bop Delta S_"ST"$, also the breaking $X_1$ of @tab:sts-X equals that of the references class by class. The plane-wave engine confirms $bop Delta S^((1))_"lit" + X_1 = 0$ on all 39 ghost-number-1 leg classes (residual $10^(-15)$).

=== Checks with complete finite parts <sec:st-checks>

The checks above test $bop Delta S^((1))_"ST" = -X_1$ with $X_1$ taken from the class results, i.e. from the $eta$-tagged coefficients $N_1$ of @sec:m-K. The following checks test the restored Slavnov--Taylor identity directly on one-loop Green's functions of $S_"ren" = S_0 + delta S^((1)) + Delta S^((1))_"ST"$, with the complete $overline("MS")$ finite parts and without the tag $eta$. They use the Feynman rules of @sec:feynman-rules and $Delta S^((1))_"ST"$ of @eq:sts-result.

#check[*Gauge two-point function, all three factors.* For the legs $(gh^a, A^b)$ the identity states that the renormalised one-loop two-point function is transverse, $p_mu Gamma^((1) mu nu)_(A^a A^b)(p) = 0$; with the field $b^a$ kept, the tree-level two-point function is transverse. All one-loop graphs were evaluated at momentum $p$ with complete finite parts, after the subtraction $delta S^((1))$ and the limit $d -> 4$. The bosonic loops (gauge bubble, ghost loop, Higgs bubble and seagull) are exactly transverse in $d$ dimensions. The fermion loops (BMHV, both chiralities) give
$ p_mu Gamma^((1) mu nu)_(A^a A^b)(p) |_"fermion loops" = -lf S^((r))/3 thin p^2 p^nu delta^(a b) , $
at three momenta, one of them timelike. The term $lf Tr_f 1/6 (partial dot Am)^2$ of @eq:sts-result has the vertex function $lf S^((r))/3 thin p^mu p^nu delta^(a b)$ and contributes $+lf S^((r))/3 thin p^2 p^nu delta^(a b)$. The sum vanishes to working precision.

*Ward identity for the legs $(gh^B, H^dagger, H)$.* The hypercharge ghost is free, so for these legs the identity is the abelian Ward identity that relates the contracted vertex function $q_mu Gamma^mu_(B H^dagger H)$ to the difference of the Higgs self-energies at $p_1^2$ and $p_2^2$ ($H(p_1)$, $H^dagger (p_2)$ incoming, $q = -p_1 - p_2$, relative signs fixed by the tree-level identity). Only graphs with a fermion line can violate it. The fermion graphs (the Higgs self-energy bubble and the triangles with $B$ attached to either fermion line) were evaluated for each Yukawa coupling $y_k$ with exact bubble and triangle finite parts at off-shell Euclidean momenta. The violation is $-lf thin N^c_k y_k^2\/3 thin (gy\/2) (p_1^2 - p_2^2)$, with $N^c_k = 3$ for $k = u, d$ and $N^c_e = 1$; summed over $k$ it is $-lf thin tilde(Y)\/3 thin (gy\/2)(p_1^2 - p_2^2)$. This is the component $-1/3 tilde(Y) N^((B))_S$ of $16 pi^2 X_1$ in @tab:sts-X. The term $r[Q^((B))_J] Q^((B))_J$ of $Delta S^((1))_"ST"$, with $r[Q^((B))_J] = tilde(Y)\/3$, cancels it to working precision.

*Ward identity for the legs $(gh^a, overline(Psi)_f, Psi_f)$.* For these legs the identity is @eq:m-glue-wFF ($Psi_f$ incoming with $p$, $overline(Psi)_f$ outgoing with $p' = p + q$). For the hypercharge ghost $gh^B$ the functions $Gamma_(srcA gh^B)$, $Gamma_(overline(srcF) gh^B Psi)$ and $Gamma_(gh^B overline(Psi) srcF)$ have no one-loop graph ($fg^(B b c) = 0$, no ghost--gauge vertex) and no counterterm, so they keep their tree values; for the gluon ghost $gh^A$ their one-loop finite parts enter, together with the graph of $Gamma_(G overline(Psi) Psi)$ with the three-gluon vertex. All one-loop functions of @eq:m-glue-wFF were evaluated with complete $overline("MS")$ finite parts at three Euclidean off-shell points $(p, p')$, at generic couplings, with a non-zero and with a vanishing Higgs mass parameter $m^2$, for $gh^B$ with $f = Q, U, D, L, E$ and for $gh^A$ with $f = Q, U$. The poles cancel. Between the physical projectors, $Pcb [dots] Pc$, the remaining violation for $S_0 + delta S^((1))$ is
$ ST(Gamma)^((1))_(gh^a overline(Psi)_f Psi_f) = -lf thin i (C_f - 1/4 C^((r))_"ad" + 1/2 nu_(f r)) thin slb(q) thin Tg^a_f thin Pc , $
with $r$ the gauge factor of $gh^a$. This is the vertex function of the component $x_(f r) N^(f r)_F$ of $X_1$ (@tab:sts-x), to $10^(-13)$ at values of order one. The terms $r[Q^(f r)_V] Q^(f r)_V$ and $r[Q^(f r)_R] Q^(f r)_R$ of $Delta S^((1))_"ST"$ (@eq:sts-r-V, @eq:sts-r-R) cancel it to the same precision.

*Massless gauge bosons in the broken phase.* The photon, $Z$--photon and gluon two-point functions at zero momentum, computed at order $m_f^2$ from the massive fermion loops with complete finite parts, vanish once the mass-like terms of $Delta S^((1))_"ST" [H_v + phi]$ ($H_v$ the constant vacuum doublet, $phi$ the fluctuation) are included (@sec:bl-massless).]

#general[At order $l$ the same steps apply to $X_l$ of @eq:m-Xl. $X_l$ is a local functional of ghost number 1 and dimension $<= 5$, expanded in the same basis $N_beta$ (@tab:Nbasis), and $bop Delta S^((l))_"ST" = -X_l$ is solved over the same $Q_alpha$ with the same normalisation conditions (N1)--(N4), which leave a unique solution. The existence conditions change: at two loops $gh^B$ can sit on an internal line, so the component of $X_l$ along $i y_k gy integral gh^B (O_k - O_k^dagger)$ need not vanish, and (N4) must then be replaced by (N4$'$); the source classes can acquire breaking from two-loop graphs. $bop X_l = 0$ remains the consistency test, and $D^(a b c) = 0$ removes the anomaly directions at every order. Three things differ from one loop. First, $X_l$ receives contributions from the classes (a)--(d) of @sec:m-beyond, including $bop_d delta S^((L))$, $bop_d Delta S^((L))$ and the products $cal(S)_"q" (C)$ of source-dependent counterterms (@eq:m-Delta-l), so the source counterterms $Q^(f r)_R$, $Q_zeta$ enter quadratically. Second, the singular counterterm $delta S^((l))$ can contain non-evanescent non-invariant poles, so its pole identity is not purely evanescent as @eq:sts-bd-dS is. Third, the evanescent parts of $delta S^((1))$ (@eq:sts-dS, including the $d$-dimensional continuation of the closed directions) and the $d$-dimensional form of $Delta S^((1))_"ST"$ enter the two-loop graphs and change $X_2$; they must be kept exactly as stated here.]
