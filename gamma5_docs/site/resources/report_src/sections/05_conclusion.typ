#import "../macros.typ": *

// local shorthands (this file only)
#let Yt = $tilde(Y)$
#let Dl(l) = $Dbrk^([#l])$
#let STq = $cal(S)_"q"$

= Conclusion <sec:conclusion>

The one-loop renormalisation of the Standard Model with one fermion family has been carried out in the unbroken phase, in the Feynman gauge, with the BMHV $g5$ and sterile partners. Eighteen classes of 1PI functions (@sec:classes) were computed diagram by diagram. The ten classes of ghost number $0$ give the singular counterterm $delta S^((1))$. Its four-dimensional part is the standard $overline("MS")$ counterterm with the known one-loop beta functions. Its evanescent part consists of the four operator types $E^((r))_A$, $tilde(E)^((r))_F$ ($r = G, W, B$), $E_H$ and $tilde(E)_H$ of @eq:evanescent, eight operators in all. The eight classes with one insertion of $Dbrk_1$ give the finite breaking $X_1$. Their poles agree with those predicted by $delta S^((1))$. $X_1$ satisfies $bop X_1 = 0$, its anomaly part vanishes because $D^(a b c) = 0$, and $bop Delta S^((1))_"ST" = -X_1$ has a unique solution under the normalisation conditions (N1)--(N4) of @sec:normalisation. The restored identity was also checked with complete finite parts on the gauge two-point functions, on the Ward identities with Higgs legs and with fermion legs, and on the massless photon, $Z$--photon and gluon two-point functions of the broken phase (@sec:st-checks). The solution agrees with @Belusca-Maito:2020ala and @Lazopoulos:2026two up to $bop$-closed terms. No step uses an anticommuting $g5$, a reading point or a finite constant taken from another scheme.

The report delivers:
- the renormalised action at one loop in the unbroken phase, with every coefficient explicit (@sec:Sren1);
- the Feynman rules of its counterterm vertices (@sec:ct-rules);
- the ten-step recipe of @sec:m-algorithm in brief, with the place where each step was carried out and what changes at two loops (@sec:recipe);
- the same renormalised action in the broken phase, in the mass eigenbasis (@sec:Sren-broken), and the points to watch in going there (@sec:broken).

== The renormalised action at one loop <sec:Sren1>

The renormalised action to one-loop order is
#eqbox[$ S_"ren" = S_0 + delta S^((1)) + Delta S^((1))_"ST" . $ <eq:c-Sren>]
This subsection states its three parts completely. The derivations are in @sec:st-solution. All coefficients are exact.

*Conventions for coefficients.* The renormalisation constants $delta Z_x$, $delta y_k$, $delta lambda$, $delta m^2$ are counterterm coefficients: each includes the factor $lf thin 1\/ep$. Where a table lists a number without this factor, its header says that it gives $16 pi^2 ep$ times the constant. The finite counterterm is $Delta S^((1))_"ST" = lf sum_alpha r_alpha Q_alpha$, with $Q_alpha$ the monomials of @tab:Qbasis; the numbers $r_alpha$ do not include $1\/(16 pi^2)$.

*Symbols.* $C_f$, $S^((r))_f$ and $C^((r))_"ad"$ are the group invariants of @eq:invariants and @tab:invariants ($C^((G))_"ad" = 3 gs^2$, $C^((W))_"ad" = 2 gw^2$, $C^((B))_"ad" = 0$). $S^((r)) = S^((r))_"lh" + S^((r))_"rh"$ is the sum of $S^((r))_f$ over the left-handed and the right-handed multiplets (@eq:SLR): $S^((G)) = 2 gs^2$, $S^((W)) = 2 gw^2$, $S^((B)) = 10/3 gy^2$. $S^((r))_H$ and $C_H$ are the Higgs entries of @tab:invariants: $S^((G))_H = 0$, $S^((W))_H = gw^2\/2$, $S^((B))_H = gy^2\/2$, $C_H = 3/4 gw^2 + 1/4 gy^2$. The Yukawa and gauge-coupling combinations are
$ Yt = 3 y_u^2 + 3 y_d^2 + y_e^2, quad Y_4 = 3 y_u^4 + 3 y_d^4 + y_e^4, quad G_4 = 3 gw^4 + 2 gw^2 gy^2 + gy^4, \
  cal(Y)_f = ((y_u^2 + y_d^2)/2, thin y_u^2, thin y_d^2, thin y_e^2/2, thin y_e^2) quad "on" quad (Q, U, D, L, E) , $
as in @eq:sts-yuk. The Yukawa numbers $nu_(f r)$ of @sec:g1_cFF are $nu_(Q G) = y_u^2 + y_d^2$, $nu_(Q W) = nu_(L W) = 0$, $nu_(Q B) = 4 y_u^2 - 2 y_d^2$, $nu_(U G) = 2 y_u^2$, $nu_(U B) = y_u^2\/2$, $nu_(D G) = 2 y_d^2$, $nu_(D B) = -y_d^2$, $nu_(L B) = 2 y_e^2$, $nu_(E B) = y_e^2$.

=== Tree-level action

$S_0 = integral dif^d x thin cal(L)_0$ is @eq:L0: the $d$-dimensional gauge term @eq:Lgauge, the fermion term @eq:Lferm (equivalently @eq:Lferm-weyl, with the sterile partners and the evanescent kinetic term), the Higgs term @eq:LHiggs, the Yukawa term @eq:LYuk, the gauge fixing @eq:Lgf at $xi = 1$ and the BRST sources @eq:Lext. Its only non-invariance is $sop S_0 = Dbrk_1$ of @eq:Delta1.

=== The singular counterterm

$delta S^((1))$ has an invariant part and an evanescent non-invariant part, $delta S^((1)) = delta S^((1))_"inv" + delta S^((1))_"ev"$ (@eq:sts-dS).

*Invariant part.* With the closed directions $cal(O)$ of @sec:closed,
$ delta S^((1))_"inv" = & sum_(r = G, W) [ 1/2 delta Z^((r))_A thin cal(O)^((r))_A + delta Z_(c_r) thin cal(O)^((r))_g + 1/2 delta Z^((r))_gh thin cal(O)^((r))_c ] + delta Z^((B))_A thin Q^((B))_T \
  & + sum_f 1/2 delta Z_f thin cal(O)^f_psi + 1/2 delta Z_H thin cal(O)_H \
  & + sum_(k = u, d, e) delta y_k thin Q^+_(y_k) - delta lambda thin Q_lambda - delta m^2 thin Q_m . $ <eq:c-dS-inv>
Every operator is continued to $d$ dimensions as in @sec:evanescent: unsplit $d$-dimensional indices, fields and derivatives in the bosonic monomials, and the barred vertex $overline(psi)_f gbar^mu Am_mu psi_f$ in the fermion monomials ($overline(psi)_f thin i sl(partial) psi_f = overline(psi)_f thin i slb(partial) psi_f$ for a physical field). Equivalently, $delta S^((1))_"inv"$ is the first-order term of $S_0$ written in bare quantities,
$ & A_0 = Z_A^(1\/2) A, quad (srcA, agh, b)_0 = Z_A^(-1\/2) (srcA, agh, b), quad xi_0 = Z_A xi, quad gh_0 = Z_gh^(1\/2) gh, quad srcC_0 = Z_gh^(-1\/2) srcC, \
  & c_(r 0) = mu^ep Z_(c_r) c_r, quad psi_(f 0) = Z_f^(1\/2) psi_f, quad srcF_(f 0) = Z_f^(-1\/2) srcF_f, quad chi_(f 0) = chi_f, \
  & H_0 = Z_H^(1\/2) H, quad srcH_0 = Z_H^(-1\/2) srcH, \
  & y_(k 0) = mu^ep (y_k + delta y_k), quad lambda_0 = mu^(2 ep) (lambda + delta lambda), quad m_0^2 = m^2 + delta m^2 , $ <eq:c-bare>
with one $Z_A$, $Z_gh$, $Z_(c_r)$ per gauge factor and $Z_x = 1 + delta Z_x$. $c_r in {gs, gw, gy}$ is the coupling of the factor $r$. There is one exception: the evanescent kinetic term $i overline(psi)_f ghat^mu partial_mu chi_f + i overline(chi)_f ghat^mu partial_mu psi_f$ of @eq:Lferm-weyl is not rescaled, since $delta S^((1))$ contains no $chi_f$. The gauge fixing is not renormalised. The constants are listed in @tab:c-Z.

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([constant $x$], [$16 pi^2 ep thin x$, general form], [$16 pi^2 ep thin x$, Standard Model]),
    [$delta Z^((r))_A$], [$5/3 C^((r))_"ad" - 2/3 S^((r)) - 1/3 S^((r))_H$], [$11/3 gs^2$, $11/6 gw^2$, $-43/18 gy^2$],
    [$delta Z_(c_r)$], [$-11/6 C^((r))_"ad" + 1/3 S^((r)) + 1/6 S^((r))_H$], [$-29/6 gs^2$, $-35/12 gw^2$, $43/36 gy^2$],
    [$delta Z^((r))_gh$], [$8/3 C^((r))_"ad" - 2/3 S^((r)) - 1/3 S^((r))_H$], [$20/3 gs^2$, $23/6 gw^2$ ($r = G, W$)],
    [$delta Z_(c_r) + 1/2 delta Z^((r))_gh$], [$-1/2 C^((r))_"ad"$], [$-3/2 gs^2$, $-gw^2$ ($r = G, W$)],
    [$delta Z_f$], [$-(C_f + cal(Y)_f)$], [$Q$: $-(4/3 gs^2 + 3/4 gw^2 + 1/36 gy^2 + (y_u^2 + y_d^2)\/2)$],
    [], [], [$U$: $-(4/3 gs^2 + 4/9 gy^2 + y_u^2)$, #h(0.3em) $D$: $-(4/3 gs^2 + 1/9 gy^2 + y_d^2)$],
    [], [], [$L$: $-(3/4 gw^2 + 1/4 gy^2 + y_e^2\/2)$, #h(0.3em) $E$: $-(gy^2 + y_e^2)$],
    [$delta Z_H$], [$2 C_H - Yt$], [$3/2 gw^2 + 1/2 gy^2 - Yt$],
    [$(delta y_u)\/y_u$], [], [$9/4 y_u^2 + 3/4 y_d^2 + 1/2 y_e^2 - 4 gs^2 - 9/8 gw^2 - 17/24 gy^2$],
    [$(delta y_d)\/y_d$], [], [$3/4 y_u^2 + 9/4 y_d^2 + 1/2 y_e^2 - 4 gs^2 - 9/8 gw^2 - 5/24 gy^2$],
    [$(delta y_e)\/y_e$], [], [$3/2 y_u^2 + 3/2 y_d^2 + 5/4 y_e^2 - 9/8 gw^2 - 15/8 gy^2$],
    [$delta lambda$], [], [$12 lambda^2 - 6 lambda C_H + 2 lambda Yt - Y_4 + 3/16 G_4$],
    [$delta m^2$], [], [$m^2 (6 lambda - 3 C_H + Yt)$],
  ),
  caption: [Renormalisation constants of @eq:c-dS-inv and @eq:c-bare ($overline("MS")$, Feynman gauge; derived in @sec:st-dS). Each constant $x$ includes the factor $lf thin 1\/ep$; the second and third columns give $16 pi^2 ep thin x$. The three entries of a row are $r = G, W, B$. In the abelian factor $cal(O)^((B))_A$ and $cal(O)^((B))_c$ have zero coefficient, $Z_(c_B) (Z^((B))_A)^(1\/2) = 1$, and $gh^B$ is free. $delta Z_(c_r) + 1/2 delta Z^((r))_gh$ is the common coefficient of all source counterterms.],
) <tab:c-Z>

*Evanescent part.* The operators of @eq:evanescent are
$ & E^((r))_A = 1/2 integral dif^d x thin bar(A)^a_mu hat(partial)^2 bar(A)^(a mu), quad
  E_H = integral dif^d x thin (hat(partial)_mu H)^dagger hat(partial)^mu H, \
  & tilde(E)^((r))_F = -1/4 integral dif^d x (bar(F)^a_(mu nu) bar(F)^(a mu nu) - F^a_(mu nu) F^(a mu nu)), quad
  tilde(E)_H = integral dif^d x [(bar(D)_mu H)^dagger bar(D)^mu H - (D_mu H)^dagger D^mu H] , $ <eq:c-E>
where $bar(F)$, $bar(D)$ have barred indices and $bar(A)$ in place of $A$, and $a$ runs over the factor $r$. Then
#eqbox[$ delta S^((1))_"ev" = - lf 1/ep [ sum_(r = G, W, B) S^((r)) (2/3 tilde(E)^((r))_F + 1/3 E^((r))_A) + Yt (tilde(E)_H + 2/3 E_H) ] . $ <eq:c-dS-ev>]
In the Standard Model the coefficients of $(E^((G))_A, E^((W))_A, E^((B))_A, E_H)$ are $-lf 1/ep (2/3 gs^2, 2/3 gw^2, 10/9 gy^2, 2/3 Yt)$, and those of $(tilde(E)^((G))_F, tilde(E)^((W))_F, tilde(E)^((B))_F, tilde(E)_H)$ are $-lf 1/ep (4/3 gs^2, 4/3 gw^2, 20/9 gy^2, Yt)$. The coefficient of $tilde(E)^((r))_F$ is the same in the two-, three- and four-point functions. The $tilde(E)$ terms arise as follows. The graphs without a closed fermion loop give the continued directions of @eq:c-dS-inv. The closed fermion loops give $-2/3 S^((r)) integral dif^d x (-1/4 bar(F)^a_(mu nu) bar(F)^(a mu nu))$ and $-Yt integral dif^d x (bar(D)_mu H)^dagger bar(D)^mu H$ (times $lf 1\/ep$), built from barred objects only. Rewriting these two barred invariants as continued ones produces the $tilde(E)$ terms of @eq:c-dS-ev.

*Variation in $d$ dimensions.* $cal(O)^((r))_A$, $Q^((B))_T$, $cal(O)^f_psi$, $cal(O)_H$, $Q^+_(y_k)$, $Q_lambda$, $Q_m$, $tilde(E)^((r))_F$ and $tilde(E)_H$ are $bop_d$-closed. $cal(O)^((r))_g$ and $cal(O)^((r))_c$ are not: $bop_d cal(O)^((r))_g = bop_d cal(O)^((r))_c = Dbrk^((r))_1$ (@eq:bd-Oc). With @eq:b-evanescent and @eq:Nhat,
$ bop_d delta S^((1)) = - lf 1/ep [ & sum_(r = G, W) 1/2 C^((r))_"ad" Dbrk^((r))_1 \
  & + sum_(r = G, W, B) 1/3 S^((r)) (hat(N)^((r))_A + hat(N)^((r))_(A A)) + 2/3 Yt hat(N)_H ] , $ <eq:c-bd-dS>
which is purely evanescent and cancels the poles of all ghost-number-1 classes (@tab:sts-poles).

=== The restoring counterterm

Under (N1)--(N4) the unique solution of $bop Delta S^((1))_"ST" = -X_1$ is @eq:sts-result (derivation in @sec:st-solve). In matrix-valued notation, with $Am$ and $wm$ of @eq:matrix-fields in the representation of the multiplet at hand and $wm_C = sum_(r = G, W) C^((r))_"ad" wm^((r))$,
#eqbox[$ Delta S^((1))_"ST" = lf integral dif^4 x { & sum_f Tr_f [ 1/6 (partial dot Am)^2 - i/3 (partial_mu Am_nu) [Am^mu, Am^nu] + 1/3 (Am dot Am)^2 - 1/4 Am_mu Am_nu Am^mu Am^nu ] \
  & - sum_f overline(psi)_f gamma^mu [C_f Am_mu + 1/2 sum_(r in f) nu_(f r) Am^((r))_mu] psi_f
    - 1/4 sum_f [ i overline(srcF)_f wm_C psi_f - i overline(psi)_f wm_C srcF_f ] \
  & + Yt/3 thin i [(partial_mu H)^dagger Am^mu H - H^dagger Am^mu partial_mu H] - Yt H^dagger Am^((W)) dot Am^((W)) H \
  & - (y_u^2 + 2 y_d^2) H^dagger {Am^((W))_mu, Am^((B) mu)} H - (13/3 y_u^2 + 7/3 y_d^2 + 3 y_e^2) H^dagger Am^((B)) dot Am^((B)) H \
  & - 1/2 (y_u^2 + y_d^2) gs^2 H^dagger H thin G^A dot G^A } . $ <eq:c-DS>]
In the Higgs terms $Am = Am^((W)) + Am^((B))$ on the doublet. All fields, indices, derivatives and Dirac matrices are four-dimensional. The coefficients $r_alpha$ in the basis $Q_alpha$ are listed in @tab:c-summary. The only term with a ghost is the fermion-source term. There is no term with the ghost-sector sources $tilde(srcA)$ or $srcC$, no Higgs-source ($srcH$) term, no $gh^B$, no fermion kinetic, Yukawa-vertex, $H^4$ or $m^2$ term and no Levi-Civita tensor $epsilon^(mu nu alpha beta)$.

*Continuation to $d$ dimensions.* In higher orders @eq:c-DS enters as $integral dif^d x$ of the same expression with every object barred: $Am_mu -> Tg^a_f bar(A)^a_mu$, $partial_mu -> bar(partial)_mu$, $gamma^mu -> gbar^mu$, $psi_f = P_(c_f) Psi_f$, and all Lorentz contractions four-dimensional, so that $Evan Delta S^((1))_"ST" = 0$ (@sec:m-objects). This is a choice. It does not change any one-loop result, and it does change $X_2$ (@sec:m-beyond). In this form only the source term has a $bop_d$ variation with an evanescent part: $bop_d Q^(f r)_R$ contains $Dbrk^((f r))_1$ (@eq:bd-Oc, @eq:m-bYR), so that
$ Evan bop_d Delta S^((1))_"ST" = - 1/4 lf sum_(r = G, W) C^((r))_"ad" Dbrk^((r))_1 . $ <eq:c-bd-DS>

=== The scheme

- *Singular part.* $delta S^((1))$ is the pole part of the ghost-number-0 functions at $eta = ep$, with $overline("MS")$ for all parameters (@eq:m-Z1). It contains no finite term.
- *Finite part.* $Delta S^((1))_"ST"$ is purely Slavnov--Taylor restoring. Its invariant part is fixed by (N1)--(N4) of @sec:normalisation: zero coefficients for $Q^f_K$, $Q_(S S)$, $Q_m$, $Q_lambda$, $Q^plus.minus_(y_k)$ (N1), $Q^((r))_T$ (N2), $Q^((r))_(rho partial)$ and $Q_rho [fg_((r))]$ (N3), and every monomial containing $gh^B$ (N4). With these conditions $bop$ has rank 64 on the 64 remaining monomials, so the solution is unique. It exists: $X_1$ vanishes in the source classes and has no component along $i y_k gy integral gh^B (O_k - O_k^dagger)$ (@sec:g1_cFFS), so (N4$'$) is not needed.
- *Other representatives.* The representatives of @Belusca-Maito:2020ala and @Lazopoulos:2026two differ from @eq:c-DS by the $bop$-closed terms of @eq:sts-lit-shift, plus a sterile kinetic term in @Belusca-Maito:2020ala. These are finite renormalisations along the closed directions of @sec:closed, and they leave $X_1$ unchanged.

*Beta functions.* The constants of @tab:c-Z give the one-loop beta functions (derived in @sec:st-dS). They agree with the known one-loop values for one family @Machacek:1983tz:
$ 16 pi^2 beta_(gs) & = -29/3 gs^3, quad 16 pi^2 beta_(gw) = -35/6 gw^3, quad 16 pi^2 beta_(gy) = 43/18 gy^3, \
  16 pi^2 beta_(y_u) & = y_u (9/2 y_u^2 + 3/2 y_d^2 + y_e^2 - 8 gs^2 - 9/4 gw^2 - 17/12 gy^2), \
  16 pi^2 beta_(y_d) & = y_d (3/2 y_u^2 + 9/2 y_d^2 + y_e^2 - 8 gs^2 - 9/4 gw^2 - 5/12 gy^2), \
  16 pi^2 beta_(y_e) & = y_e (3 y_u^2 + 3 y_d^2 + 5/2 y_e^2 - 9/4 gw^2 - 15/4 gy^2), \
  16 pi^2 beta_lambda & = 24 lambda^2 - 3 lambda (3 gw^2 + gy^2) + 4 lambda Yt - 2 Y_4 + 3/8 G_4, \
  16 pi^2 thin (dif ln m^2) / (dif ln mu) & = 12 lambda + 2 Yt - 9/2 gw^2 - 3/2 gy^2 . $ <eq:c-beta>

=== All operators

@tab:c-summary lists every operator of $S_"ren" - S_0$ with its monomial and its two coefficients; the monomials $Q_alpha$ are those of @tab:Qbasis, written out in the second column. The divergent column gives the coefficient of the operator in $16 pi^2 ep thin delta S^((1))$, written with the group invariants and the Yukawa combinations defined above. The finite column gives $r_alpha$, the coefficient in $16 pi^2 Delta S^((1))_"ST"$. For the $Q_alpha$ the divergent column is the coefficient of the monomial of @tab:Qbasis in the four-dimensional limit. In $d$ dimensions it multiplies the continued monomial of @eq:c-dS-inv, and the four $tilde(E)$, $E$ rows complete $delta S^((1))$. In the rows $Q^((r))_T$ and $Q^f_K$ the divergent entry is $16 pi^2 ep thin delta Z^((r))_A$ and $16 pi^2 ep thin delta Z_f$; the other entries follow from @eq:c-dS-inv and @tab:c-Z.

#[
#show figure: set block(breakable: true)
#figure(
  {
    set par(justify: false)
    set text(size: 8.4pt)
    show table.cell.where(y: 0): set par(leading: 0.8em)
    show table.cell.where(x: 0): set par(leading: 0.8em)
    table(
    columns: (auto, 1.6fr, 1.05fr, 0.8fr),
    align: (left, left, left, left),
    inset: (x: 3.5pt, y: 3.4pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else if y in (13, 15, 27, 33) { (top: 0.3pt + luma(160)) } else { none },
    table.header([operator], [monomial: integrand of $integral dif^4 x$], [divergent: coefficient in $16 pi^2 ep thin delta S^((1))$], [finite: $r_alpha$, coefficient in $16 pi^2 Delta S^((1))_"ST"$]),
    [$Q^((r))_T$], [$-1/4 (partial_mu A^a_nu - partial_nu A^a_mu)^2$], [$5/3 C^((r))_"ad" - 2/3 S^((r)) - 1/3 S^((r))_H$: \ $11/3 gs^2$, $11/6 gw^2$, $-43/18 gy^2$], [$0$ (N2)],
    [$Q^((r))_L$], [$-1/2 (partial dot A^a)^2$], [$0$], [$-1/3 S^((r))$: $-2/3 gs^2$, $-2/3 gw^2$, $-10/9 gy^2$],
    [$Q^((r))_M$], [$1/2 A^a dot A^a$], [$0$], [$0$],
    [$Q^((r))_(3 f)$], [$-fg^(a b c) (partial_mu A^a_nu) A^(b mu) A^(c nu)$], [$2/3 C^((r))_"ad" - 2/3 S^((r)) - 1/3 S^((r))_H$: \ $2/3 gs^2$, $-1/6 gw^2$], [$-1/3 S^((r))$: $-2/3 gs^2$, $-2/3 gw^2$],
    [$Q_(3 d)$, $Q_(3 B)$, \ $Q^((r))_(B i)$, $Q^((r))_(B epsilon)$, \ $Q^((r))_(4 epsilon)$], [$gs^3 d^(A B C) (partial dot G^A)(G^B dot G^C)$; #h(0.2em) $gy^3 (partial dot B)(B dot B)$; #h(0.2em) $c_r^2 gy times {(partial dot B)(A^a dot A^a), thin (partial dot A^a)(A^a dot B), thin A^(a mu) A^(a nu) partial_mu B_nu}$ ($i = 1, 2, 3$); #h(0.2em) $c_r^2 gy thin epsilon^(mu nu alpha beta) A^a_mu partial_nu A^a_alpha B_beta$; #h(0.2em) $c_r^3 gy thin epsilon^(mu nu alpha beta) f^(a b c) A^a_mu A^b_nu A^c_alpha B_beta$, with $f^(I J K) = epsilon^(I J K)$ for $"SU"(2)$], [$0$], [$0$],
    [$Q_(G 4 a), Q_(G 4 b), Q_(G 4 c)$], [$gs^4 times {(G^A dot G^A)^2, thin (G^A dot G^B)^2, thin d^(A B E) d^(C D E) (G^A dot G^B)(G^C dot G^D)}$], [$7/24, -7/12, 7/8$], [$11/72, -1/4, 5/12$],
    [$Q_(W 4 a), Q_(W 4 b)$], [$gw^4 times {(W^I dot W^I)^2, thin (W^I dot W^J)^2}$], [$13/24, -13/24$], [$7/24, -1/4$],
    [$Q_(B 4)$], [$gy^4 (B dot B)^2$], [$0$], [$95/648$],
    [$Q^((G W))_a, Q^((G W))_b$], [$c_r^2 c_s^2 times {(A^a_((r)) dot A^a_((r)))(A^b_((s)) dot A^b_((s))), thin (A^a_((r)) dot A^b_((s)))^2}$, $r s = G W$], [$0$], [$1/24, 1/12$],
    [$Q^((G B))_a, Q^((G B))_b$], [the same with $r s = G B$], [$0$], [$11/216, 11/108$],
    [$Q^((W B))_a, Q^((W B))_b$], [the same with $r s = W B$], [$0$], [$1/36, 1/18$],
    [$Q_(G G G B)$], [$gs^3 gy thin d^(A C D) (G^A dot G^C)(G^D dot B)$], [$0$], [$1/18$],
    [$Q^f_K$], [$overline(psi)_f i sl(partial) psi_f$], [$-(C_f + cal(Y)_f)$], [$0$ (N1)],
    [$Q^(f r)_V$], [$overline(psi)_f gamma^mu Am^((r))_mu psi_f$], [$-(C_f + cal(Y)_f) - C^((r))_"ad"$], [$-C_f - 1/2 nu_(f r)$],
    [$Q_(S S)$], [$(partial_mu H)^dagger partial^mu H$], [$2 C_H - Yt$], [$0$ (N1)],
    [$Q_m$], [$H^dagger H$], [$m^2 (C_H - 6 lambda)$], [$0$ (N1)],
    [$Q_lambda$], [$(H^dagger H)^2$], [$2 lambda C_H - 12 lambda^2 + Y_4 - 3/16 G_4$], [$0$ (N1)],
    [$Q^+_(y_u)$], [$-(O_u + O_u^dagger)$, #h(0.2em) $O_u = overline(Q) tilde(H) PR U$], [$-y_u (3/4 gw^2 + 16/3 gs^2 + 25/36 gy^2 + y_d^2)$], [$0$ (N1)],
    [$Q^+_(y_d)$], [$-(O_d + O_d^dagger)$, #h(0.2em) $O_d = overline(Q) H PR D$], [$-y_d (3/4 gw^2 + 16/3 gs^2 + 1/36 gy^2 + y_u^2)$], [$0$ (N1)],
    [$Q^+_(y_e)$], [$-(O_e + O_e^dagger)$, #h(0.2em) $O_e = overline(L) H PR E$], [$-y_e (3/4 gw^2 + 9/4 gy^2)$], [$0$ (N1)],
    [$Q^-_(y_k)$, $Q^((r))_D$], [$-i (O_k - O_k^dagger)$; #h(0.2em) $H^dagger (partial dot Am^((r))) H$], [$0$], [$0$],
    [$Q^((W))_J$, $Q^((B))_J$], [$i [(partial_mu H)^dagger Am^((r) mu) H - H^dagger Am^((r) mu) partial_mu H]$], [$C^((W))_"ad" - 2 C_H + Yt$, #h(0.3em) $-2 C_H + Yt$], [$Yt\/3$, $Yt\/3$],
    [$Q_(W W)$], [$H^dagger Am^((W)) dot Am^((W)) H$], [$2 C_H - Yt - 2 C^((W))_"ad"$], [$-Yt$],
    [$Q_(B B)$], [$H^dagger Am^((B)) dot Am^((B)) H$], [$2 C_H - Yt$], [$-(13/3 y_u^2 + 7/3 y_d^2 + 3 y_e^2)$],
    [$Q_(W B)$], [$H^dagger {Am^((W))_mu, Am^((B) mu)} H$], [$2 C_H - Yt - C^((W))_"ad"$], [$-(y_u^2 + 2 y_d^2)$],
    [$Q_(G G)$], [$gs^2 H^dagger H thin G^A dot G^A$], [$0$], [$-1/2 (y_u^2 + y_d^2)$],
    [$Q^((r))_(rho partial)$], [$tilde(srcA)^mu_a partial_mu gh^a$], [$1/2 C^((r))_"ad"$], [$0$ (N3)],
    [$Q_rho [fg_((r))]$, \ $Q_zeta [-1/2 fg_((r))]$], [$fg^(a b c)_((r)) tilde(srcA)^mu_a A^b_mu gh^c$; #h(0.2em) $-1/2 fg^(a b c)_((r)) srcC_a gh^b gh^c$], [$-1/2 C^((r))_"ad"$], [$0$ (N3); $0$],
    [$Q^(f r)_R$ ($r = G, W$)], [$i overline(srcF)_f wm^((r)) psi_f + i gh^a overline(psi)_f Tg^a_((r)) srcF_f$], [$-1/2 C^((r))_"ad"$], [$-1/4 C^((r))_"ad"$],
    [$Q^((W))_sigma$], [$i srcH^dagger wm^((W)) H - i gh^a H^dagger Tg^a_((W)) srcH$], [$-1/2 C^((W))_"ad" = -gw^2$], [$0$],
    [other $Q_rho [t]$, $Q_zeta [t]$], [$t^(a b c) tilde(srcA)^mu_a A^b_mu gh^c$; #h(0.2em) $t^(a b c) srcC_a gh^b gh^c$ (tensors $t$ of @sec:Qbasis)], [$0$], [$0$],
    [$Q^(f B)_R$, $Q^((B))_sigma$], [$Q^(f r)_R$ and $Q^((r))_sigma$ at $r = B$], [$0$], [$0$ (N4)],
    [$E^((r))_A$; #h(0.2em) $E_H$], [@eq:c-E, $integral dif^d x$: $1/2 bar(A)^a_mu hat(partial)^2 bar(A)^(a mu)$; #h(0.2em) $(hat(partial)_mu H)^dagger hat(partial)^mu H$], [$-1/3 S^((r))$; #h(0.2em) $-2/3 Yt$], [$0$],
    [$tilde(E)^((r))_F$; #h(0.2em) $tilde(E)_H$], [@eq:c-E, $integral dif^d x$: $-1/4 (bar(F)^a_(mu nu) bar(F)^(a mu nu) - F^a_(mu nu) F^(a mu nu))$; #h(0.2em) $(bar(D)_mu H)^dagger bar(D)^mu H - (D_mu H)^dagger D^mu H$], [$-2/3 S^((r))$; #h(0.2em) $-Yt$], [$0$],
  )},
  caption: [All operators of $delta S^((1))$ and $Delta S^((1))_"ST"$. The second column gives each monomial $Q_alpha$ of @tab:Qbasis as the integrand of $integral dif^4 x$, with the notation of @sec:Qbasis: four-dimensional fields, indices, derivatives and Dirac matrices, bars omitted; the quartic and Higgs--gluon monomials contain their couplings. The last two rows are the $d$-dimensional evanescent operators of @eq:c-E. Both coefficient columns are free of the factor $1\/(16 pi^2)$, and the divergent column is also free of $1\/ep$. Where three values are given they are $r = G, W, B$; for $Q^((r))_(3 f)$ and the ghost-sector rows $r = G, W$ ($C^((r))_"ad"\/2 = 3/2 gs^2, gw^2$). "(N$i$)" marks a coefficient fixed to zero by a normalisation condition. Derivations: @eq:sts-dS, @tab:sts-Z, @tab:sts-r.],
) <tab:c-summary>
]

#include "05a_ct_feynman_rules.typ"

#include "05r_recipe.typ"

#include "05c_broken_lagrangian.typ"

#include "05b_broken_phase.typ"
