#import "../macros.typ": *

== Counterterm operators and the action of #bop <sec:bases>

This section collects the Slavnov--Taylor algebra used to solve $bop Delta S^((1))_"ST" = -X_1$ (@eq:m-X1) at one loop. The formulas are those of a construction valid at every order. Every statement below was checked numerically with exact Grassmann bookkeeping, by evaluating the local functionals on plane-wave configurations of all fields and sources. This numerical program is called the plane-wave engine in the rest of the report. A second, independent implementation in explicit components reproduces the signs of @sec:st-signs, the dictionary examples, @tab:Nbasis, @eq:btable-gauge to @eq:btable-ghost, the closed directions and the normalisation statements; the hermiticity and completeness statements were checked only with the first implementation (@sec:supporting).

=== The operator #bop on fields and sources <sec:st-signs>

*Derivatives.* For a Grassmann-odd field or source $Phi$, $cal(F)[Phi + delta Phi] = cal(F) + integral delta Phi thin delta cal(F) \/ delta Phi = cal(F) + integral (delta_R cal(F) \/ delta Phi) thin delta Phi$. For an even $cal(F)$, $delta_R cal(F) \/ delta Phi = - delta cal(F) \/ delta Phi$. #sop is a left graded derivation, $sop cal(F) = integral sum_Phi (sop Phi) thin delta cal(F) \/ delta Phi$.

*Checks of @sec:brst.* Term by term, @eq:gf-exact reproduces @eq:Lgf and @eq:Lext, @eq:ST gives @eq:ST-S0, @eq:b-sources is the linearisation of @eq:ST at $S_0$, and $bop^2 = 0$ in four dimensions. The signs of @eq:ST and @eq:gf-exact matter: with $(delta_R cal(F) \/ delta Psi_f)(delta cal(F) \/ delta overline(srcF)_f)$ in place of the fermion term, with a left derivative on the odd $srcH$ in the $H^dagger$ term, or with the opposite sign of $overline(Psi)_f srcF_f$ or $H^dagger srcH$ in @eq:gf-exact, $ST(S_0) - sop S_0$ is nonzero.

*Explicit equations of motion.* @eq:b-sources and @eq:rhotilde give $bop agh^a = b^a$, $bop b^a = 0$ and $bop tilde(srcA)^mu_a = delta tilde(S)_0 \/ delta A^a_mu$: the $b$-terms cancel. Here $tilde(S)_0$ is $S_0$ without $b^a partial^mu A^a_mu + xi b^a b^a\/2$ and with $srcA -> tilde(srcA)$. Written out in four dimensions (bars omitted, $psi_f = P_(c_f) Psi_f$), with $cal(J)^mu_a$ the gauge-field equation of motion:
$ bop tilde(srcA)^mu_a &= cal(J)^mu_a + fg^(x a e) tilde(srcA)^mu_x gh^e, \
  cal(J)^mu_a &= (D_nu F^(nu mu))^a + sum_f overline(psi)_f gamma^mu Tg^a_f psi_f + i [H^dagger Tg^a_H D^mu H - (D^mu H)^dagger Tg^a_H H], $ <eq:b-rho>
$ bop srcC_a = (D_mu tilde(srcA)^mu)_a + fg^(a b c) srcC_b gh^c + sum_f i (overline(srcF)_f Tg^a_f psi_f + overline(psi)_f Tg^a_f srcF_f) - i (srcH^dagger Tg^a_H H + H^dagger Tg^a_H srcH), $ <eq:b-zeta>
$ bop overline(srcF)_f = i (partial_mu overline(psi)_f) gamma^mu - overline(psi)_f gamma^mu Am_mu + (delta cal(L)_"Yuk")/(delta psi_f) - i overline(srcF)_f wm, quad
  bop srcF_f = - i sl(D) psi_f - (delta cal(L)_"Yuk")/(delta overline(psi)_f) + i wm srcF_f, $ <eq:b-R>
$ bop srcH^dagger = - (D_mu D^mu H)^dagger - m^2 H^dagger - 2 lambda (H^dagger H) H^dagger + (delta cal(L)_"Yuk")/(delta H) + i srcH^dagger wm, \
  bop srcH = (D_mu D^mu H) + m^2 H + 2 lambda (H^dagger H) H - (delta cal(L)_"Yuk")/(delta H^dagger) + i wm srcH . $ <eq:b-sigma>
In @eq:b-zeta, $(D_mu tilde(srcA)^mu)_a = partial_mu tilde(srcA)^mu_a + fg^(a b c) A^b_mu tilde(srcA)^(mu)_c$.

#check[In four dimensions $ST(S_0) = 0$ holds with @eq:ST on 17 classes of external legs with ghost number 1. With either of the sign changes just mentioned it fails on every class with fermion or $H^dagger$ legs. $ST(S_0 + T) = sop T$ holds for a source-free, non-invariant $T$. $bop^2 = 0$ holds on every field and source, and flipping the sign of one equation in @eq:b-sources violates it. $bop G = ST(S_0 + G) - ST(S_0) - ST(G)$ holds for an even $G$ that contains every source.]

=== Diagrammatic form of #bop <sec:dictionary>

*Vertex functions.* For a local functional $Y$ and legs $Phi_k$ with polarisations $e_k$ and incoming momenta $p_k$ ($sum_k p_k = 0$), the vertex function $V_Y (Phi_1 (e_1, p_1), dots, Phi_n (e_n, p_n))$ is the coefficient of $theta_1 dots theta_n$ in $Y[sum_k theta_k e_k e^(-i p_k dot x)]$. Here $theta_k$ is a Grassmann-odd number for an odd field and an ordinary number otherwise, and the overall factor $(2 pi)^d delta^((d))(sum_k p_k)$ is stripped. The Feynman rule of $Y$ (@sec:feynman-rules) is $i V_Y$. For $bop Y$ the ghost leg $gh(u, q)$, with adjoint polarisation $u^c$ and incoming momentum $q$, is listed first. Then
$ V_(bop Y) (gh(u, q); Phi_1, dots, Phi_n) = "(i)" + "(ii)" + "(iii)" , $ <eq:dict>
with three contributions.

(i) *Abelian part.* Add a gauge leg $A(e = -i q u, thin q)$ to the vertex function of $Y$, summed over all positions:
$ "(i)" = V_Y (A^c_mu (q) "with" e_mu = -i q_mu u^c; thin Phi_1, dots, Phi_n) . $ <eq:dict-i>
In Feynman-rule language: take the $(n+1)$-point Feynman rule of $Y$ and contract the gauge index $mu$ of $A^c$ with $-i q_mu u^c$, where $q$ is the incoming ghost momentum.

(ii) *Homogeneous part.* Rotate one leg of the $n$-point vertex function by the ghost and shift its momentum:
$ "(ii)" = sum_(i=1)^n s_i thin V_Y (Phi_1, dots, [op("rot")_u e_i] (p_i + q), dots, Phi_n), $ <eq:dict-ii>
where $op("rot")_u e$ is the polarisation rotated by the ghost polarisation $u$: $[op("rot")_u e]^a = fg^(a b c) e^b u^c$ on a gauge leg and on a ghost leg, $i u^c Tg^c_f P_(c_f) e$ on $Psi_f$, $i e P_(overline(c)_f) Tg^c_f u^c$ on $overline(Psi)_f$ with $s_i = -1$, $i u^c Tg^c_H e$ on $H$, $-i e Tg^c_H u^c$ on $H^dagger$, and $s_i = +1$ on all legs except $overline(Psi)_f$. Since $Y$ is rigidly invariant, @eq:dict-ii vanishes at $q = 0$. Only the dependence on the shift $p_i -> p_i + q$ survives; in position space the derivative acts on the ghost.

*Several external ghosts.* If the vertex function of $bop Y$ has further external ghosts besides $gh(u, q)$, each external ghost in turn is the ghost of (i) and of (ii) applied to the legs that are not ghosts, times the sign $op("sgn")$ of the permutation that moves it to the first position. A ghost leg of $Y$ is instead replaced by any two external ghosts $gh(u_i, p_i)$, $gh(u_j, p_j)$ with $i < j$: one ghost leg of polarisation $-fg^(b c e) u_i^c u_j^e$ and momentum $p_i + p_j$, listed first, times the sign $op("sgn")$ of the permutation that takes $(i, j, "others")$ to the original order. This term comes from the non-linear part $sop gh^b = -1/2 fg^(b c e) gh^c gh^e$. With a single further ghost it is the rotation of the ghost leg in @eq:dict-ii, counted once. The same rules are stated after @eq:m-b-rules.

(iii) *Source terms.* Let $Y$ depend on the source $Phi^*$ of a field $Phi$ (@eq:Phistar: $A^* = srcA$, $gh^* = srcC$, $Psi^* = overline(srcF)$, $overline(Psi)^* = srcF$, $H^* = srcH^dagger$, $(H^dagger)^* = srcH$); $bop Phi^*$ is given by @eq:b-sources. The contribution is built in three steps: remove the leg $Phi^* (k)$ from the vertex function of $Y$; multiply by the vertex function of a tree vertex of $S_0$ that has a leg $Phi(-k)$, and remove that leg as well (no propagator joins the two); keep the sign of @eq:b-sources. The remaining legs of the $S_0$ vertex become legs of $bop Y$. For classes without external sources only the ghost-independent vertices of $S_0$ enter: the kinetic two-point vertex (the equation of motion of the external leg) and the gauge, Yukawa and scalar vertices.

#[
#set par(justify: false)
*Examples* (polarisations with lower indices, $p dot e = p^mu e_mu$). All three were checked against the position-space action of #bop.
+ $Y = -1/2 integral (partial dot A^a)^2$, with $V_Y (A^a (e_1, p_1), A^b (e_2, p_2)) = delta^(a b) (p_1 dot e_1)(p_2 dot e_2)$. (i): $V_(bop Y) (gh(u, q), A(e, -q)) = i q^2 (q dot e^a) u^a$. (ii): $V_(bop Y) (gh(u, q), A^b (e_1, p_1), A^e (e_2, p_2)) = fg^(b e c) u^c [(p_2 dot e_1)(p_2 dot e_2) - (p_1 dot e_1)(p_1 dot e_2)]$.
+ $Y = integral overline(psi)_f gamma^mu Am_mu psi_f$: (i) gives $V_(bop Y) (gh^a (q), overline(Psi)_f, Psi_f) = -i sl(q) Tg^a_f P_(c_f)$. On $(gh, overline(Psi), Psi, A)$ the homogeneous part vanishes: $Y$ has no derivative.
+ $Y = integral i overline(srcF)_f wm psi_f$: only (iii) contributes. By @eq:b-R, $bop overline(srcF)_f supset i (partial_mu overline(psi)_f) gamma^mu$, so $bop Y supset - integral (partial_mu overline(psi)_f) gamma^mu wm psi_f$. Let the $overline(Psi)_f$ leg be outgoing with momentum $p'$ (incoming momentum $-p'$), so that $partial_mu -> +i p'_mu$ on $overline(Psi)_f$; moving the ghost to the first position gives one minus sign: $V_(bop Y) (gh^a, overline(Psi)_f (p'), Psi_f) = i slb(p)' thin Tg^a_f P_(c_f)$, with $slb(p)' = gbar^mu p'_mu$. This is the kinetic vertex of $S_0$ glued to the source vertex. With an additional gauge leg the term $-overline(psi)_f gamma^mu Am_mu$ of $bop overline(srcF)_f$ enters, i.e. the gauge vertex of $S_0$: $V_(bop Y) (gh(u), overline(Psi)_f, Psi_f, A(e)) = i gamma^mu (Tg^a_f e^a_mu)(Tg^b_f u^b) P_(c_f)$.
]

=== Notation for the bases

All monomials of @sec:Qbasis to @sec:btable are four-dimensional, and the bars are omitted: $A^a_mu$, $partial_mu$, $gamma^mu$ stand for $bar(A)^a_mu$, $bar(partial)_mu$, $gbar^mu$, and $psi_f = P_(c_f) Psi_f$, $overline(psi)_f = overline(Psi)_f P_(overline(c)_f)$. We write $integral equiv integral dif^4 x$. For a gauge factor $r in {G, W, B}$ with coupling $c_r in {gs, gw, gy}$, $A^a_((r))$ is its gauge field ($a$ runs over the factor), $Am^((r))_mu = Tg^a_f A^a_((r) mu)$ and $wm^((r)) = Tg^a_f gh^a_((r))$ are the restrictions of @eq:matrix-fields, $Tg^a_((r))$ the restricted generators, and $fg^(a b c)_((r))$ the restricted structure constants. $d^(A B C) = 2 Tr(T^A {T^B, T^C})$ is the $"SU"(3)$ symmetric tensor. $epsilon^(mu nu alpha beta)$ is the four-dimensional Levi-Civita symbol of @eq:g5def; its indices are written $mu, nu, alpha, beta$ in this section, because $srcA$ and $srcH$ are sources. $gh^B$ and $B_mu$ always denote the hypercharge ghost and gauge field. $A, B, C, D, E$ are $"SU"(3)$ adjoint indices (@tab:alphabet), except that in a formula that contains $gh^B$ or $B_mu$ the $"SU"(3)$ indices are written $A, C, D, E$ only. A superscript or subscript $d$ of a gauge field, ghost or structure constant is an adjoint index; subscript letters in the names of monomials ($Q_(3 d)$, $N_(d f)$, $N^((r s))_d$) are labels. The Yukawa monomials are
$ O_u = overline(Q) tilde(H) PR U, quad O_d = overline(Q) H PR D, quad O_e = overline(L) H PR E, $ <eq:Ok>
so that $cal(L)_"Yuk" = - sum_k y_k (O_k + O_k^dagger)$.

=== Basis $Q_alpha$ of the restoring counterterm <sec:Qbasis>

$Q_alpha$ runs over all monomials of ghost number 0 and dimension $<= 4$ that are four-dimensional, Lorentz invariant, rigidly invariant and hermitian; they enter with real coefficients. They are built from $A$, $psi_f$, $H$, $gh$, $tilde(srcA)$, $srcC$, $srcF_f$, $srcH$. The fields $b^a$, $agh^a$ and $chi_f$ do not appear. @eq:linear-ids removes $b$ and $agh$ outside $tilde(srcA)$. $chi_f$ couples only through $Dbrk_1$, and $Pphys$ removes its external legs. There are 97 monomials (@tab:Qbasis, last column).

#[
#show figure: set block(breakable: true)
#set par(justify: false)
#figure(
  table(
    columns: (auto, 1fr, auto, auto),
    align: (left, left, left, right),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([name], [monomial], [range], [number]),
    [$Q^((r))_T$], [$-1/4 integral (partial_mu A^a_nu - partial_nu A^a_mu)^2$], table.cell(rowspan: 3, breakable: false)[$r = G, W, B$], [3],
    [$Q^((r))_L$], [$-1/2 integral (partial dot A^a)^2$], [3],
    [$Q^((r))_M$], [$1/2 integral A^a dot A^a$], [3],
    [$Q^((r))_(3 f)$], [$- integral fg^(a b c) (partial_mu A^a_nu) A^(b mu) A^(c nu)$], [$r = G, W$], [2],
    [$Q_(3 d)$], [$gs^3 integral d^(A B C) (partial dot G^A)(G^B dot G^C)$], [], [1],
    [$Q_(3 B)$], [$gy^3 integral (partial dot B)(B dot B)$], [], [1],
    [$Q^((r))_(B 1)$], [$c_r^2 gy integral (partial dot B)(A^a dot A^a)$], table.cell(rowspan: 4, breakable: false)[$r = G, W$], [2],
    [$Q^((r))_(B 2)$], [$c_r^2 gy integral (partial dot A^a)(A^a dot B)$], [2],
    [$Q^((r))_(B 3)$], [$c_r^2 gy integral A^(a mu) A^(a nu) partial_mu B_nu$], [2],
    [$Q^((r))_(B epsilon)$], [$c_r^2 gy integral epsilon^(mu nu alpha beta) A^a_mu partial_nu A^a_alpha B_beta$], [2],
    [$Q_(G 4 a)$], [$gs^4 integral (G^A dot G^A)^2$], [], [1],
    [$Q_(G 4 b)$], [$gs^4 integral (G^A dot G^B)^2$], [], [1],
    [$Q_(G 4 c)$], [$gs^4 integral d^(A B E) d^(C D E) (G^A dot G^B)(G^C dot G^D)$], [], [1],
    [$Q_(W 4 a)$], [$gw^4 integral (W^I dot W^I)^2$], [], [1],
    [$Q_(W 4 b)$], [$gw^4 integral (W^I dot W^J)^2$], [], [1],
    [$Q_(B 4)$], [$gy^4 integral (B dot B)^2$], [], [1],
    [$Q^((r s))_a$], [$c_r^2 c_s^2 integral (A^a_((r)) dot A^a_((r)))(A^b_((s)) dot A^b_((s)))$], table.cell(rowspan: 2, breakable: false)[$r s = G W, G B, W B$], [3],
    [$Q^((r s))_b$], [$c_r^2 c_s^2 integral (A^a_((r)) dot A^b_((s)))^2$], [3],
    [$Q_(G G G B)$], [$gs^3 gy integral d^(A C D) (G^A dot G^C)(G^D dot B)$], [], [1],
    [$Q^((r))_(4 epsilon)$], [$c_r^3 gy integral epsilon^(mu nu alpha beta) f^(a b c) A^a_mu A^b_nu A^c_alpha B_beta$, with $f^(I J K) = epsilon^(I J K)$ for $"SU"(2)$], [$r = G, W$], [2],
    [$Q^f_K$], [$integral overline(psi)_f i sl(partial) psi_f$], [all $f$], [5],
    [$Q^(f r)_V$], [$integral overline(psi)_f gamma^mu Am^((r))_mu psi_f$], [all $f$; $r in f$], [10],
    [$Q_(S S)$], [$integral (partial_mu H)^dagger partial^mu H$], [], [1],
    [$Q_m$], [$integral H^dagger H$], [], [1],
    [$Q_lambda$], [$integral (H^dagger H)^2$], [], [1],
    [$Q^((r))_J$], [$i integral [(partial_mu H)^dagger Am^((r) mu) H - H^dagger Am^((r) mu) partial_mu H]$], table.cell(rowspan: 2, breakable: false)[$r = W, B$], [2],
    [$Q^((r))_D$], [$integral H^dagger (partial dot Am^((r))) H$], [2],
    [$Q_(W W)$], [$integral H^dagger Am^((W)) dot Am^((W)) H$], [], [1],
    [$Q_(B B)$], [$integral H^dagger Am^((B)) dot Am^((B)) H$], [], [1],
    [$Q_(W B)$], [$integral H^dagger {Am^((W))_mu, Am^((B) mu)} H$], [], [1],
    [$Q_(G G)$], [$gs^2 integral H^dagger H thin G^A dot G^A$], [], [1],
    [$Q^+_(y_k)$], [$- integral (O_k + O_k^dagger)$], table.cell(rowspan: 2, breakable: false)[$k = u, d, e$], [3],
    [$Q^-_(y_k)$], [$-i integral (O_k - O_k^dagger)$], [3],
    [$Q^((r))_(rho partial)$], [$integral tilde(srcA)^mu_a partial_mu gh^a$], [$r = G, W, B$], [3],
    [$Q_rho [t]$], [$integral t^(a b c) tilde(srcA)^mu_a A^b_mu gh^c$], [tensors $t$ below], [10],
    [$Q_zeta [t]$], [$integral t^(a b c) srcC_a gh^b gh^c$], [tensors $t$ below], [4],
    [$Q^(f r)_R$], [$integral [i overline(srcF)_f wm^((r)) psi_f + i gh^a overline(psi)_f Tg^a_((r)) srcF_f]$], [all $f$; $r in f$], [10],
    [$Q^((r))_sigma$], [$integral [i srcH^dagger wm^((r)) H - i gh^a H^dagger Tg^a_((r)) srcH]$], [$r = W, B$], [2],
  ),
  caption: [The basis $Q_alpha$ (97 monomials). The subscripts $T$, $L$, $M$ name the transverse, longitudinal and mass two-point structures. "$r in f$" means that $f$ carries a non-trivial representation of the factor $r$ (@tab:fields). Last column: number of monomials of the row.],
) <tab:Qbasis>
]

*Tensors of the ghost-sector monomials.* $Q_rho [t]$ is defined for $t^(a b c) = fg^(a b c)_((G))$, $fg^(a b c)_((W))$, $gs d^(a b c)$ (all in $"SU"(3)$), and for the seven Kronecker structures $(tilde(srcA), A, gh) = (W^I, W^I, B)$, $(B, W^I, W^I)$, $(W^I, B, W^I)$, $(G^A, G^A, B)$, $(B, G^A, G^A)$, $(G^A, B, G^A)$, $(B, B, B)$. $Q_zeta [t]$ is defined for $t = -1/2 fg_((G))$, $-1/2 fg_((W))$, and $(srcC^A gh^A gh^B)$, $(srcC^I gh^I gh^B)$. The tree-level source terms of $S_0$ are $Q_rho [fg_((r))] + Q^((r))_(rho partial)$, $Q_zeta [-1/2 fg_((r))]$, $Q^(f r)_R$ and $Q^((r))_sigma$.

*Completeness.* Gauge monomials: rigid invariance allows the tensors $fg$, $d^(A B C)$ and $delta$ in the combinations listed. Only $"SU"(3)$ has a $d$-tensor. $(G G W)$, $(G W B)$, $(G B B)$ and $(W B B)$ have no invariant, because the generators are traceless. Lorentz structures were counted in momentum space. Three gauge legs with one momentum span six tensors. Their Bose-symmetric parts give one structure for a totally symmetric or a totally antisymmetric group tensor, and three for $(A_((r)) A_((r)) B)$. The $epsilon$ structures vanish for totally symmetric and totally antisymmetric tensors. For four legs, two $"SU"(3)$ identities reduce the tensors with the symmetry of $(A dot A)(A dot A)$ to the three of $Q_(G 4 a)$, $Q_(G 4 b)$, $Q_(G 4 c)$: $d^(A B E)d^(C D E) + d^(A C E)d^(B D E) + d^(A D E)d^(B C E) = 1/3 (delta^(A B)delta^(C D) + delta^(A C)delta^(B D) + delta^(A D)delta^(B C))$, and the expression of $f f$ through $delta delta$ and $d d$ (used in @eq:QYM). The four-point $epsilon$ structures need a totally antisymmetric tensor: $f^(a b c)$ with one abelian leg. Matter: same-chirality bilinears allow only $overline(psi)_f gamma^mu psi_f$, hypercharge forbids bilinears that mix multiplets, and the Yukawa list is complete. This basis is also complete for the structures generated by fermion loops: the traces $Tr_f (Tg^(a_1) dots Tg^(a_n))$ of @eq:invariants decompose into exactly these tensors. The $Q_alpha$ are linearly independent. The image of #bop has dimension 72, the kernel has dimension 25 (@sec:closed), and $72 + 25 = 97$.

*Hermiticity.* With $gh^dagger = gh$, $agh^dagger = -agh$, $tilde(srcA)^dagger = -tilde(srcA)$, $srcC^dagger = -srcC$, $overline(srcF) = srcF^dagger gamma^0$, every monomial of @tab:Qbasis is hermitian. The other relative sign in $Q^(f r)_R$ or $Q^((r))_sigma$ gives anti-hermitian monomials. Hermiticity of the action excludes them.

=== The invariant directions <sec:closed>

The kernel of #bop on $"span"{Q_alpha}$ has dimension 25. It consists of the 20 directions of the renormalisation of fields and parameters, the three imaginary Yukawa couplings and two charge deformations of $"U"(1)_Y$. With $cal(Q)^((r))_"YM"$ the quartic Yang--Mills term,
$ cal(Q)^((G))_"YM" = -1/(4 gs^2) (1/2 Q_(G 4 a) - Q_(G 4 b) + 3/2 Q_(G 4 c)), quad cal(Q)^((W))_"YM" = -1/(4 gw^2) (Q_(W 4 a) - Q_(W 4 b)), $ <eq:QYM>
which uses the $"SU"(3)$ identity
$ fg^(E A B) fg^(E C D) (G^A dot G^C)(G^B dot G^D) = gs^2 [1/2 (G^A dot G^A)^2 - (G^A dot G^B)^2 \
  + 3/2 d^(A B E) d^(C D E) (G^A dot G^B)(G^C dot G^D)], $
and the Higgs pieces $h_W = -Q^((W))_J + 2 Q_(W W) + Q_(W B)$, $h_B = -Q^((B))_J + 2 Q_(B B) + Q_(W B)$, $h_G = 0$, the 25 closed operators are:
#[
#show math.equation: set block(breakable: true)
$ cal(O)^((r))_A &= 2 Q^((r))_T + 3 Q^((r))_(3 f) + 4 cal(Q)^((r))_"YM" + sum_(f in r) Q^(f r)_V + h_r - Q^((r))_(rho partial) quad (r = G, W), \
  cal(O)^((r))_c &= Q^((r))_(rho partial) + Q_rho [fg_((r))] + Q_zeta [-1/2 fg_((r))] + sum_(f in r) Q^(f r)_R + delta_(r W) Q^((W))_sigma quad (r = G, W), \
  cal(O)^((r))_g &= Q^((r))_(3 f) + 2 cal(Q)^((r))_"YM" + sum_(f in r) (Q^(f r)_V + Q^(f r)_R) + h_r + delta_(r W) Q^((W))_sigma \
  &quad + Q_rho [fg_((r))] + Q_zeta [-1/2 fg_((r))] quad (r = G, W), \
  & Q^((B))_T, quad cal(O)^((B))_A = 2 Q^((B))_T + sum_f Q^(f B)_V + h_B - Q^((B))_(rho partial), quad cal(O)^((B))_c = Q^((B))_(rho partial) + sum_f Q^(f B)_R + Q^((B))_sigma , \
  cal(O)^f_psi &= 2 (Q^f_K + sum_(r in f) Q^(f r)_V) + sum_(k: thin f in O_k) y_k Q^+_(y_k), \
  cal(O)_H &= 2 (Q_(S S) - Q^((W))_J - Q^((B))_J + Q_(W W) + Q_(B B) + Q_(W B)) - 2 m^2 Q_m - 4 lambda Q_lambda + sum_k y_k Q^+_(y_k), \
  & Q^+_(y_k), quad Q^-_(y_k), quad Q_lambda, quad Q_m, quad
  cal(O)_q = sum_f q_f / Y_f (Q^(f B)_V + Q^(f B)_R) quad "for" q = q_"baryon", thin q_"lepton" . $ <eq:closed>
]
Here $q_"baryon" = (1/3, 1/3, 1/3, 0, 0)$ and $q_"lepton" = (0, 0, 0, 1, 1)$ on $(Q, U, D, L, E)$, with $q_H = 0$. $cal(O)^((r))_A = bop integral tilde(srcA)^mu_a A^a_((r) mu)$ and $cal(O)^((r))_c = -bop integral srcC_a gh^a_((r))$ are exact. $cal(O)^f_psi = -bop integral (overline(srcF)_f psi_f - overline(psi)_f srcF_f)$ and $cal(O)_H = bop integral (srcH^dagger H - H^dagger srcH)$ are exact as well. $cal(O)^((r))_g = c_r partial S_0 \/ partial c_r$. In the abelian factor, $Q^((B))_T$ is separately closed and $cal(O)^((B))_g = cal(O)^((B))_A + cal(O)^((B))_c - 2 Q^((B))_T$. $cal(O)_q$ gauges a conserved $"U"(1)$ with the field $B$, i.e. it shifts the hypercharges by $q$; for $q = Y$ it is $cal(O)^((B))_g$. The 25 operators are closed, linearly independent and span the kernel (all three statements checked numerically).

=== Evanescent operators with one-loop poles <sec:evanescent>

Fermion loops couple only to $bar(A)$, $bar(partial)$ and to the physical $psi_f$. Gauge, ghost and scalar loops are $d$-dimensional and covariant. One-loop poles therefore contain evanescent structures only through hatted external momenta. Since $bar(A)$ has no hatted index, $hat(partial)$ can only appear as $hat(partial)^2$, i.e. in two-point functions. $overline(psi)_f ghat^mu psi_f = 0$ for one chirality. The evanescent operators that can carry one-loop poles are, with $r = G, W, B$,
$ E^((r))_A &= 1/2 integral dif^d x thin bar(A)^a_mu hat(partial)^2 bar(A)^(a mu), quad E_H = integral dif^d x thin (hat(partial)_mu H)^dagger hat(partial)^mu H, \
  tilde(E)^((r))_F &= -1/4 integral dif^d x (bar(F)^a_(mu nu) bar(F)^(a mu nu) - F^a_(mu nu) F^(a mu nu)), quad
  tilde(E)_H = integral dif^d x [(bar(D)_mu H)^dagger bar(D)^mu H - (D_mu H)^dagger D^mu H] . $ <eq:evanescent>
$bar(F)$ and $bar(D)$ have barred indices and $bar(A)$ in place of $A$. $tilde(E)^((r))_F$ and $tilde(E)_H$ are $bop_d$-closed. $E^((r))_A$ and $E_H$ are not:
$ bop_d E^((r))_A = integral dif^d x (bar(partial)^mu gh^a + fg^(a b c) bar(A)^(b mu) gh^c) hat(partial)^2 bar(A)^a_mu, \
  bop_d E_H = i integral dif^d x [ (hat(partial)_mu H)^dagger hat(partial)^mu (wm H) - (hat(partial)_mu (H^dagger wm)) hat(partial)^mu H ] . $ <eq:b-evanescent>
They are the two non-invariant evanescent singular counterterms. Fermion bilinears cannot carry an evanescent pole: $chi_f$ has no tree vertex, and $P_(overline(c)_f) ghat^mu P_(c_f) = 0$.

*Continuation to $d$ dimensions.* A closed direction of @eq:closed is continued with the field content of $S_0$: unsplit $d$-dimensional indices, fields and derivatives in the bosonic monomials, and the barred vertex $overline(psi)_f gbar^mu Am_mu psi_f$ in the fermion monomials ($overline(psi)_f thin i sl(partial) psi_f = overline(psi)_f thin i slb(partial) psi_f$ for a physical field). The poles of the graphs without a closed fermion loop are exactly such combinations. A closed fermion loop gives the same invariants built from $bar(A)$ and $bar(partial)$ only. The difference to the continued form is $tilde(E)^((r))_F$ or $tilde(E)_H$, with one coefficient for the two-, three- and four-point parts. In this form the evanescent part of $delta S^((1))$ is a combination of the four operators of @eq:evanescent (@sec:st-solution).

*Closed directions in $d$ dimensions.* Since $ST(S_0) = Dbrk_1 != 0$, $bop_d^2$ need not vanish on functionals that depend on $srcC$, $overline(srcF)_f$ or $srcF_f$: these are the sources of the fields on which $Dbrk_1$ depends. $cal(O)^((r))_A$, $cal(O)_H$, $cal(O)^f_psi$ (in its source-free form of @eq:closed), $Q^((B))_T$, $Q^plus.minus_(y_k)$, $Q_lambda$, $Q_m$, $tilde(E)^((r))_F$ and $tilde(E)_H$ remain $bop_d$-closed. The two directions with ghost sources are not:
$ bop_d cal(O)^((r))_g = bop_d cal(O)^((r))_c = Dbrk^((r))_1 = sum_f Dbrk^((f r))_1 , \
  Dbrk^((f r))_1 = integral dif^d x thin gh^a_((r)) [ overline(Psi)_f Tg^a_f ghat^mu P_(overline(c)_f) partial_mu Psi_f + (partial_mu overline(Psi)_f) ghat^mu P_(c_f) Tg^a_f Psi_f ] , $ <eq:bd-Oc>
with $gh^a_((r))$ the ghosts of the factor $r$, so that $Dbrk_1 = sum_(f, r) Dbrk^((f r))_1$. The first equality follows from $partial ST(S_0) \/ partial c_r = bop_d (partial S_0 \/ partial c_r)$ and the linearity of $Dbrk_1$ in $Tg_((r))$, the second from $cal(O)^((r))_c = - bop_d integral srcC_a gh^a_((r))$. Of the monomials of @tab:Qbasis only $Q^(f r)_R$ has a $bop_d$ variation with a $Dbrk_1$ part: on $gh overline(Psi) Psi$, $V_(bop_d Q^(f r)_R) = V_(Dbrk^((f r))_1) + V_(N^(f r)_F)$ (@eq:m-bYR). Hence $bop_d delta S^((1))$ contains, besides the variations @eq:b-evanescent, the term $Dbrk^((r))_1$ times the common coefficient of the source counterterms.

*Names.* The ghost-number-1 evanescent monomials that carry the one-loop poles of $Dbrk_1 dot Gamma$ are $Dbrk^((f r))_1$ and
$ hat(N)^((r))_A = integral dif^d x thin (bar(partial)^mu gh^a) hat(partial)^2 bar(A)^a_mu, quad
  hat(N)^((r))_(A A) = integral dif^d x thin fg^(a b c) bar(A)^(b mu) gh^c thin hat(partial)^2 bar(A)^a_mu, quad
  hat(N)_H = bop_d E_H , $ <eq:Nhat>
so that $bop_d E^((r))_A = hat(N)^((r))_A + hat(N)^((r))_(A A)$. Their vertex functions are listed with the ghost $gh^a (q)$ first. All momenta are incoming, except in the last line, where $p'$ is the outgoing momentum of $overline(Psi)_f$ as in @eq:V-delta (its incoming momentum is $-p'$, so $q = p' - p$):
$ V_(hat(N)_A) (gh^a (q), A^b_nu (-q)) &= i hat(q)^2 bar(q)^nu delta^(a b), \
  V_(hat(N)_(A A)) (gh^a (q), A^b_mu (p_1), A^c_nu (p_2)) &= -fg^(a b c) gmb^(mu nu) (hat(p)_1^2 - hat(p)_2^2), \
  V_(hat(N)_H) (gh^a (q), H_j (p_1), H^dagger_i (p_2)) &= -i (hat(p)_1^2 - hat(p)_2^2) (Tg^a_H)_(i j), \
  V_(Dbrk^((f r))_1) (gh^a (q), overline(Psi)_f (p'), Psi_f (p)) &= -i Tg^a_f (slh(p) P_(overline(c)_f) - slh(p)' P_(c_f)) , $ <eq:Nhat-V>
the last one being the Feynman rule @eq:V-delta divided by $i$.

=== Basis $N_beta$ of the breaking <sec:Nbasis>

$N_beta$ runs over four-dimensional monomials of ghost number 1 and dimension 5. The basis of @tab:Nbasis spans $bop("span"{Q_alpha})$ in every class of external legs. In the gauge classes $gh A A$ it is the complete basis with an undifferentiated ghost, and it contains the $epsilon$ structures that carry the anomaly. Since $bop X_1 = 0$, the one-loop breaking lies in $bop("span"{Q_alpha})$ plus the anomaly directions; a component outside this span signals an anomaly or an error. The table also contains monomials that no $bop Q_alpha$ produces, so that a component of $X_1$ along them would be seen. Among them are: in the class $gh A A A$ the monomial $N_(epsilon d f)$, which multiplies the $"SU"(3)^3$ anomaly coefficient; in the class $gh overline(psi) psi$ the ten $N^(f r)_E$; in the class $gh H^dagger H$ the four $N^((r))_(S E)$ and $N^((r))_(S m)$; in the class $gh H^dagger H H^dagger H$ the monomials $N^((r))_(S lambda)$, since no $Q_alpha$ contributes there, $bop("span"{Q_alpha}) = 0$. There are 244 independent $N_beta$. Class names: c = ghost, V = gauge field, F = fermion, S = Higgs; the number after the class name counts its monomials, and the last column counts the monomials of the row.

#[
#show figure: set block(breakable: true)
#set par(justify: false)
#let cl(n, body) = table.cell(rowspan: n, breakable: false, body)
#figure(
  table(
    columns: (2.1cm, auto, 1fr, auto),
    align: (left, left, left, right),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([class], [name], [monomial], [number]),
    // ---- cV
    cl(2)[cV (6)], [$N^((r))_3$], [$integral gh^a partial^2 partial dot A^a$], [3],
    [$N^((r))_1$], [$integral gh^a partial dot A^a$], [3],
    // ---- cVV
    cl(2)[cVV, $fg$ (4)], [$N^((r))_(f 1)$], [$integral fg^(a b c) gh^a A^b dot partial^2 A^c$ ($r = G, W$)], [2],
    [$N^((r))_(f 2)$], [$integral fg^(a b c) gh^a A^b_mu partial^mu partial dot A^c$ ($r = G, W$)], [2],
    cl(4)[cVV, symmetric (20)], table.cell(colspan: 2)[$S_i [t] = integral t^(a b c) gh^a {A^b dot partial^2 A^c, thin A^b_mu partial^mu partial dot A^c, thin partial_mu A^b_nu partial^mu A^(c nu), thin partial_mu A^b_nu partial^nu A^(c mu), thin (partial dot A^b)(partial dot A^c)}$, $i = 1, dots, 5$, for the tensors:], [],
    [$N_(d i)$], [$t = gs^3 d^(a b c)$ (all in $"SU"(3)$)], [5],
    [$N_(B i)$], [$t = gy^3$ on $(gh^B, B, B)$], [5],
    [$N^((r))_(B r i)$], [$t = c_r^2 gy delta^(b c)$ on $(gh^B, A_((r)), A_((r)))$, $r = G, W$], [10],
    cl(1)[cVV, mixed (14)], [$N^((r))_(r B i)$], [$c_r^2 gy integral gh^a {B dot partial^2 A^a, thin B_mu partial^mu partial dot A^a, thin A^a dot partial^2 B, thin A^a_mu partial^mu partial dot B, thin partial_mu A^a_nu partial^mu B^nu, thin partial_mu A^a_nu partial^nu B^mu, thin (partial dot A^a)(partial dot B)}$, $i = 1, dots, 7$, $r = G, W$], [14],
    cl(4)[cVV, $epsilon$ (6)], [$N^((r))_(epsilon B r r)$], [$c_r^2 gy integral epsilon^(mu nu alpha beta) gh^B partial_mu A^a_nu partial_alpha A^a_beta$ ($r = G, W$)], [2],
    [$N^((r))_(epsilon r B)$], [$c_r^2 gy integral epsilon^(mu nu alpha beta) gh^a partial_mu A^a_nu partial_alpha B_beta$ ($r = G, W$)], [2],
    [$N_(epsilon d)$], [$gs^3 integral d^(A B C) epsilon^(mu nu alpha beta) gh^A partial_mu G^B_nu partial_alpha G^C_beta$], [1],
    [$N_(epsilon B)$], [$gy^3 integral epsilon^(mu nu alpha beta) gh^B partial_mu B_nu partial_alpha B_beta$], [1],
    // ---- cVVV
    cl(13)[cVVV (34)], [$N^((r))_(f f)$], [$integral fg^(e a b) fg^(e c y)(partial_mu gh^a) A^b_nu A^(c mu) A^(y nu)$ ($r = G, W$)], [2],
    [$N^((r))_(B f 3)$], [$integral gh^B fg^(a b c) A^(a mu) A^(b nu) partial_mu A^c_nu$ ($r = G, W$)], [2],
    [$N_(d f)$], [$gs^3 integral d^(A B C) fg^(A D E) (partial gh^E dot G^D)(G^B dot G^C)$], [1],
    [$N^((r))_(r B f)$], [$c_r^2 gy integral fg^(a b c)(partial gh^a dot A^b)(A^c dot B)$ ($r = G, W$)], [2],
    [$N^((r))_(epsilon r B f)$], [$c_r^2 gy integral epsilon^(mu nu alpha beta) fg^(a b c)(partial_mu gh^a) A^b_nu A^c_alpha B_beta$ ($r = G, W$)], [2],
    [$N^((r))_(epsilon 4 B)$], [$c_r^3 gy integral epsilon^(mu nu alpha beta) f^(a b c) A^a_mu A^b_nu A^c_alpha partial_beta gh^B$ ($r = G, W$)], [2],
    [$N_(G 4 a), N_(G 4 b)$], [$gs^4 integral (partial gh^A dot G^A)(G^B dot G^B)$, #h(0.4em) $gs^4 integral (partial gh^A dot G^B)(G^A dot G^B)$], [2],
    [$N_(W 4 a)$], [$gw^4 integral (partial gh^I dot W^I)(W^J dot W^J)$], [1],
    [$N_(B 4)$], [$gy^4 integral (partial gh^B dot B)(B dot B)$], [1],
    [$N^((r s))_(a, b, c, d)$], [$c_r^2 c_s^2 integral {(partial gh^a_((r)) dot A^a_((r)))(A^b_((s)) dot A^b_((s))), thin (A^a_((r)) dot A^a_((r)))(partial gh^b_((s)) dot A^b_((s))), thin (partial gh^a_((r)) dot A^b_((s)))(A^a_((r)) dot A^b_((s))), thin (A^a_((r)) dot partial gh^b_((s)))(A^a_((r)) dot A^b_((s)))}$, $r s = G W, G B, W B$], [12],
    [$N_(G G G B 1, 2, 3)$], [$gs^3 gy integral d^(A C D) {(partial gh^A dot G^C)(G^D dot B), thin (G^A dot G^C)(G^D dot partial gh^B), thin (G^A dot G^C)(partial gh^D dot B)}$], [3],
    [$N_cal(J) [t]$], [$t = gs d$, $(W^I, B, W^I)$, $(G^A, B, G^A)$], [3],
    [$N_(epsilon d f)$], [$gs^4 integral epsilon^(mu nu alpha beta) d^(A E D) f^(C C' E) (partial_mu gh^A) G^C_nu G^(C')_alpha G^D_beta$; no $bop Q_alpha$ contains it, it enters only the anomaly (@sec:st-solution)], [1],
    // ---- cVVVV
    cl(2)[cVVVV (5)], [$N^((r))_(B Y M)$], [$integral gh^B fg^(e a b) fg^(e c y) A^a_mu A^b_nu A^(c mu) A^(y nu)$ ($r = G, W$)], [2],
    [$N_cal(J) [t]$], [$t = gs d$, $(W^I, B, W^I)$, $(G^A, B, G^A)$], [3],
    // ---- cFF
    cl(2)[cFF (20)], [$N^(f r)_F$], [$integral overline(psi)_f gamma^mu (partial_mu wm^((r))) psi_f$ (all $f$; $r in f$)], [10],
    [$N^(f r)_E$], [$i integral (partial_mu overline(psi)_f) gamma^mu wm^((r)) psi_f$ (all $f$; $r in f$)], [10],
    // ---- cFFV
    cl(3)[cFFV (35)], [$N^(f r s)_(F V)$], [$integral overline(psi)_f gamma^mu wm^((r)) Am^((s))_mu psi_f$ (all $f$; all $r, s in f$)], [22],
    [$N^(f r)_(V F)$], [$integral overline(psi)_f gamma^mu Am^((r))_mu wm^((r)) psi_f$ ($r = G, W$; $r in f$)], [5],
    [$N_cal(J) [t]$], [$t = gs d$ ($f = Q, U, D$), $(B, W^I, W^I)$ ($f = U, D, E$), $(B, G^A, G^A)$ ($f = L, E$)], [8],
    // ---- cFFS
    cl(2)[cFFS (16)], [$N^(O r)_(Y L)$, $N^(O r)_(Y R)$], [$i y_k integral O[overline(psi) -> overline(psi) wm^((r))]$, #h(0.4em) $i y_k integral O[psi -> wm^((r)) psi]$, for $O = O_k, O_k^dagger$ and $r = G, W$: per $O$ one gluon structure ($N^(O G)_(Y R) = N^(O G)_(Y L)$; none for $O_e$, $O_e^dagger$) and one $W$ structure (the $W$ ghost on the doublet leg)], [10],
    [$N^O_(Y B)$], [$i y_k gy integral gh^B O$ for $O = O_k, O_k^dagger$], [6],
    // ---- cSS
    cl(4)[cSS (8)], [$N^((r))_S$], [$i integral [(partial_mu H)^dagger (partial^mu wm^((r))) H - H^dagger (partial^mu wm^((r))) partial_mu H]$ ($r = W, B$)], [2],
    [$N^((r))_(S 2)$], [$integral H^dagger (partial^2 wm^((r))) H$ ($r = W, B$)], [2],
    [$N^((r))_(S E)$], [$i integral (partial^2 H)^dagger wm^((r)) H$ ($r = W, B$)], [2],
    [$N^((r))_(S m)$], [$i m^2 integral H^dagger wm^((r)) H$ ($r = W, B$)], [2],
    // ---- cSSV
    cl(8)[cSSV (17)], [$N^((r s))_(S V)$], [$integral H^dagger {partial_mu wm^((r)), Am^((s) mu)} H$ ($r, s = W, B$)], [4],
    [$N^((s r))_(S V E 1)$], [$i integral (partial_mu H)^dagger Am^((s) mu) wm^((r)) H$ ($r, s = W, B$)], [4],
    [$N^((s r))_(S V E 2)$], [$i integral H^dagger (partial dot Am^((s))) wm^((r)) H$ ($r, s = W, B$)], [4],
    [$N^((W W))_(S V E 3)$], [$i integral H^dagger wm^((W)) Am^((W))_mu partial^mu H$], [1],
    [$N^((W W))_(S V E 4)$], [$i integral H^dagger wm^((W)) (partial dot Am^((W))) H$], [1],
    [$N_(S C)$], [$-i integral H^dagger [Am^((W))_mu, partial^mu wm^((W))] H$], [1],
    [$N_(S G)$], [$gs^2 integral H^dagger H thin partial gh^A dot G^A$], [1],
    [$N_cal(J) [t]$], [$t = (B, G^A, G^A)$], [1],
    // ---- cSSVV
    cl(3)[cSSVV (9)], [$N^((s t r))_(S V V 1)$], [$i integral H^dagger Am^((s)) dot Am^((t)) wm^((r)) H$ ($s t = W W, W B, B B$; $r = W, B$)], [6],
    [$N^((W B W))_(S V V 2)$], [$i integral H^dagger wm^((W)) Am^((W)) dot Am^((B)) H$], [1],
    [$N_cal(J) [t]$], [$t = (B, W^I, W^I)$, $(B, G^A, G^A)$], [2],
    // ---- cSSSS
    cl(1)[cSSSS (2)], [$N^((r))_(S lambda)$], [$i lambda integral (H^dagger H)(H^dagger wm^((r)) H)$ ($r = W, B$)], [2],
    // ---- ccRF
    cl(2)[ccRF (22)], [$N^(f r s)_(R R)$], [$integral overline(srcF)_f wm^((r)) wm^((s)) psi_f$], [11],
    [$overline(N)^(f r s)_(R R)$], [$integral overline(psi)_f wm^((r)) wm^((s)) srcF_f$; in both, $r, s in f$, $r <= s$ in the order $G, W, B$, $(r, s) != (B, B)$], [11],
    // ---- cc sigma H
    cl(2)[cc$srcH$H (4)], [$N^((r s))_sigma$], [$integral srcH^dagger wm^((r)) wm^((s)) H$ ($r s = W W, W B$)], [2],
    [$overline(N)^((r s))_sigma$], [$integral H^dagger wm^((r)) wm^((s)) srcH$ ($r s = W W, W B$)], [2],
    // ---- cc rho, ccc zeta
    cl(4)[cc$srcA$, ccc$srcC$ (22)], [$N_(rho partial) [t]$], [$integral t^(a b c) tilde(srcA)^mu_a gh^b partial_mu gh^c$ for the 10 tensors $t$ of $Q_rho [t]$], [10],
    [$N_(rho A) [T]$], [$integral T^(a b c e) tilde(srcA)^mu_a A^b_mu gh^c gh^e$ with $T^(a b c e) = t^(a x e) fg^(x b c)$, for $t = fg_((G))$, $fg_((W))$, $gs d$, $(W^I, W^I, B)$, $(B, W^I, W^I)$, $(G^A, G^A, B)$, $(B, G^A, G^A)$ (for the other three $t$ it vanishes)], [7],
    [$N_(rho A) [T']$], [the same with $T'^(a b c e) = 1/2 t^(a b x) fg^(x c e)$, for $t = gs d$, $(W^I, B, W^I)$, $(G^A, B, G^A)$], [3],
    [$N^((r))_zeta$], [$integral fg^(a b c)_((r)) srcC_a gh^b gh^c gh^B$ ($r = G, W$)], [2],
  ),
  caption: [The basis $N_beta$ (244 monomials). $N_cal(J) [t]$ denotes the part of $integral t^(a b c) gh^c A^b_mu cal(J)^mu_a$ (@eq:b-rho) in the class of its row; the same name in different classes denotes different monomials. The 20 monomials $N_(rho partial) [t]$, $N_(rho A) [T]$, $N_(rho A) [T']$ span all structures of these two forms. In the tensor labels $(X, Y, Z)$ the three entries are the fields carrying the indices $a$, $b$, $c$.],
) <tab:Nbasis>
]

=== The action of #bop on the basis <sec:btable>

The coefficients $B_(alpha beta)$ in $bop Q_alpha = sum_beta B_(alpha beta) N_beta$ are given below. They were obtained by fitting the vertex functions of $bop Q_alpha$ (computed with @eq:dict) on random kinematics in every class, at generic and at unit couplings. The fits reproduce every $bop Q_alpha$ exactly.

#block(sticky: true)[*Gauge sector* ($r = G, W, B$; terms with $fg_((B))$ are absent):]
$ bop Q^((r))_M &= - N^((r))_1, quad bop Q^((r))_L = - N^((r))_3 - N^((r))_(f 2), quad bop Q^((r))_T = - N^((r))_(f 1) + N^((r))_(f 2), \
  bop Q^((r))_(3 f) &= N^((r))_(f 1) - N^((r))_(f 2) + N^((r))_(f f), \
  bop Q_(3 d) &= 2(N_(d 1) - N_(d 2) + N_(d 3) - N_(d 5)) + N_(d f), quad bop Q_(3 B) = 2 (N_(B 1) - N_(B 2) + N_(B 3) - N_(B 5)), \
  bop Q^((r))_(B 1) &= 2 N^((r))_(B r 1) + 2 N^((r))_(B r 3) - 2 N^((r))_(r B 4) - 2 N^((r))_(r B 7), \
  bop Q^((r))_(B 2) &= - N^((r))_(B r 2) - N^((r))_(B r 5) + N^((r))_(r B 1) - N^((r))_(r B 2) + N^((r))_(r B 3) + 2 N^((r))_(r B 5) - N^((r))_(r B 7) - N^((r))_(r B f), \
  bop Q^((r))_(B 3) &= 2 N^((r))_(B r 2) + N^((r))_(B r 4) + N^((r))_(B r 5) - N^((r))_(r B 3) - N^((r))_(r B 4) - N^((r))_(r B 5) - N^((r))_(r B 6), \
  bop Q^((r))_(B epsilon) &= N^((r))_(epsilon B r r) - N^((r))_(epsilon r B) - N^((r))_(epsilon r B f), \
  bop Q_(G 4 a) &= 4 N_(G 4 a), quad bop Q_(G 4 b) = 4 N_(G 4 b), quad bop Q_(G 4 c) = 8/3 gs^2 N^((G))_(f f) - 4/3 N_(G 4 a) + 8/3 N_(G 4 b), \
  bop Q_(W 4 a) &= 4 N_(W 4 a), quad bop Q_(W 4 b) = - 4 gw^2 N^((W))_(f f) + 4 N_(W 4 a), quad bop Q_(B 4) = 4 N_(B 4), \
  bop Q^((r s))_a &= 2 (N^((r s))_a + N^((r s))_b), quad bop Q^((r s))_b = 2 (N^((r s))_c + N^((r s))_d), \
  bop Q_(G G G B) &= 2 N_(G G G B 1) + N_(G G G B 2) + N_(G G G B 3), quad bop Q^((r))_(4 epsilon) = 3 N^((r))_(epsilon r B f) + N^((r))_(epsilon 4 B) . $ <eq:btable-gauge>
The four-point monomials have no derivative, so only @eq:dict-i contributes; with @eq:dict-ii the three-point monomials pick up $N_(f f)$, $N_(d f)$, $N_(r B f)$, $N_(epsilon r B f)$.

#block(sticky: true)[*Fermion sector* ($r in f$):]
$ bop Q^f_K &= sum_(r in f) N^(f r)_F, quad bop Q^(f r)_V = - N^(f r)_F, \
  bop Q^(f r)_R &= N^(f r)_F + i (N^(f r r)_(F V) - N^(f r)_(V F)) + N^(f r r)_(R R) + overline(N)^(f r r)_(R R) + [bop Q^(f r)_R]_"Yuk" quad (r = G, W), \
  bop Q^(f B)_R &= N^(f B)_F + [bop Q^(f B)_R]_"Yuk" . $ <eq:btable-fermion>
The Yukawa parts are the ghost rotations of the $f$ legs of the Yukawa monomials. For the abelian factor, let $q_f (O_k)$ and $q_H (O_k)$ be the hypercharges with which the fields of $f$ and of the Higgs doublet enter $O_k$: $q = -Y_f$ for $overline(psi)_f$, $q = Y_(f')$ for $psi_(f')$, $q = Y_H$ for $H$ and $q = -Y_H$ for $tilde(H)$ (zero if the field does not occur in $O_k$). Then
$ [bop Q^(f B)_R]_"Yuk" = - sum_(k) q_f (O_k) (N^(O_k)_(Y B) - N^(O_k^dagger)_(Y B)), quad
  [bop Q^((B))_sigma]_"Yuk" = - sum_k q_H (O_k) (N^(O_k)_(Y B) - N^(O_k^dagger)_(Y B)) . $ <eq:btable-yuk-B>
#block(sticky: true)[For gluons and $W$:]
$ [bop Q^(Q G)_R]_"Yuk" &= -N^(O_u G)_(Y L) + N^(O_u^dagger G)_(Y L) - N^(O_d G)_(Y L) + N^(O_d^dagger G)_(Y L), \
  [bop Q^(U G)_R]_"Yuk" &= N^(O_u G)_(Y L) - N^(O_u^dagger G)_(Y L), quad [bop Q^(D G)_R]_"Yuk" = N^(O_d G)_(Y L) - N^(O_d^dagger G)_(Y L), \
  [bop Q^(Q W)_R]_"Yuk" &= -N^(O_u W)_(Y L) + N^(O_u^dagger W)_(Y R) - N^(O_d W)_(Y L) + N^(O_d^dagger W)_(Y R), \
  [bop Q^(L W)_R]_"Yuk" &= -N^(O_e W)_(Y L) + N^(O_e^dagger W)_(Y R) . $ <eq:btable-yuk-GW>
For gluons the rotations of the two fermion legs coincide by rigid invariance ($H$ is a colour singlet), $N^(O G)_(Y R) = N^(O G)_(Y L)$.

#block(sticky: true)[*Higgs sector:*]
$ bop Q_(S S) &= N^((W))_S + N^((B))_S, quad bop Q^((r))_J = N^((r))_S + N^((W r))_(S V) + N^((B r))_(S V), quad bop Q^((W))_D = N^((W))_(S 2) + N_(S C), quad bop Q^((B))_D = N^((B))_(S 2), \
  bop Q_(W W) &= N^((W W))_(S V), quad bop Q_(B B) = N^((B B))_(S V), quad bop Q_(W B) = N^((W B))_(S V) + N^((B W))_(S V), quad bop Q_(G G) = 2 N_(S G), \
  bop Q_m &= bop Q_lambda = bop Q^plus.minus_(y_k) = 0, \
  bop Q^((W))_sigma &= N^((W))_S - i (2 N^((W W))_(S V E 1) + N^((W W))_(S V E 2) + 2 N^((W W))_(S V E 3) + N^((W W))_(S V E 4)) - N^((W B))_(S V) \
  &quad + 2 (N^((W B W))_(S V V 1) - N^((W B W))_(S V V 2)) - N^((W W))_sigma - overline(N)^((W W))_sigma \
  &quad + N^(O_u W)_(Y L) - N^(O_u^dagger W)_(Y R) + N^(O_d W)_(Y L) - N^(O_d^dagger W)_(Y R) + N^(O_e W)_(Y L) - N^(O_e^dagger W)_(Y R), \
  bop Q^((B))_sigma &= N^((B))_S - N^((B W))_(S V) - N^((B B))_(S V) + [bop Q^((B))_sigma]_"Yuk" . $ <eq:btable-higgs>

*Ghost sector.* For every tensor $t$ (and $t^(a b c) -> delta^(a b) partial$ for $Q^((r))_(rho partial)$):
$ bop Q_rho [t] &= integral t^(a b c) [ (cal(J)^mu_a + fg^(x a e) tilde(srcA)^mu_x gh^e) A^b_mu gh^c - tilde(srcA)^mu_a (D_mu gh)^b gh^c + 1/2 tilde(srcA)^mu_a A^b_mu fg^(c x y) gh^x gh^y ], \
  bop Q^((r))_(rho partial) &= integral [ (cal(J)^mu_a + fg^(x a e) tilde(srcA)^mu_x gh^e) partial_mu gh^a_((r)) + 1/2 tilde(srcA)^mu_a partial_mu (fg^(a b c) gh^b gh^c) ], \
  bop Q_zeta [t] &= integral t^(a b c) [ (bop srcC_a) gh^b gh^c - 1/2 srcC_a fg^(b x y) gh^x gh^y gh^c + 1/2 srcC_a gh^b fg^(c x y) gh^x gh^y ], $ <eq:btable-ghost>
with $bop srcC_a$ from @eq:b-zeta. The parts without sources give, for example, $bop Q^((G))_(rho partial) = N^((G))_(f 1) - N^((G))_(f 2) - N^((G))_(f f) - sum_(f = Q,U,D) N^(f G)_F + ("source terms")$. The source terms lie in the classes $gh gh tilde(srcA)$, $gh gh tilde(srcA) A$, $gh gh gh srcC$, $gh gh srcH^dagger H$, $gh gh H^dagger srcH$ and $gh gh overline(srcF) psi$. The complete numerical table, including these classes, is part of the supporting material (@sec:supporting).

=== Normalisation conditions <sec:normalisation>

The restoring counterterm is written $Delta S^((1))_"ST" = lf sum_alpha r_alpha Q_alpha$; the coefficients $r_alpha$ do not contain the factor $1\/(16 pi^2)$. $Delta S^((1))_"ST"$ is fixed by $bop Delta S^((1))_"ST" = -X_1$ only up to the 25 closed directions of @eq:closed. A normalisation condition removes a set of monomials. It makes $Delta S^((1))_"ST"$ unique if no closed direction lies in the span of the remaining ones. It still allows a solution if the image of the remaining ones contains $X_1$. The conditions are:
- (N1) no monomial that survives when all gauge fields, ghosts and sources vanish: $Q^f_K$, $Q_(S S)$, $Q_m$, $Q_lambda$, $Q^plus.minus_(y_k)$ (14 monomials). This removes the five $cal(O)^f_psi$, $cal(O)_H$, $Q^plus.minus_(y_k)$, $Q_lambda$, $Q_m$: 14 closed directions.
- (N2) no $Q^((r))_T$ (3 monomials). This removes $cal(O)^((G))_A$, $cal(O)^((W))_A$, $Q^((B))_T$. The coupling directions $cal(O)^((r))_g$ contain no two-point monomial and survive (N2).
- (N3) no $Q^((r))_(rho partial)$ ($r = G, W, B$) and no $Q_rho [fg_((r))]$ ($r = G, W$) (5 monomials). This removes $cal(O)^((G,W))_c$, $cal(O)^((G,W))_g$, $cal(O)^((B))_c$.
- (N4) no monomial containing $gh^B$: $Q_rho [(W,W,B)]$, $Q_rho [(G,G,B)]$, $Q_rho [(B,B,B)]$, $Q_zeta [(srcC^A gh^A gh^B)]$, $Q_zeta [(srcC^I gh^I gh^B)]$, $Q^(f B)_R$, $Q^((B))_sigma$ (11 monomials). After (N1)--(N3) three closed directions remain: the charge deformations $cal(O)_q$ of @eq:closed for $q = Y$ (this is $cal(O)^((B))_g$), $q_"baryon"$ and $q_"lepton"$. All three contain $Q^(f B)_R$, and (N4) removes them.

#result[*Uniqueness.* With (N1)--(N4) the kernel of #bop on the remaining 64 monomials is zero: $Delta S^((1))_"ST"$ is unique. Of the 25 closed directions, the three $Q^-_(y_k)$ are fixed by (N1), and $cal(O)_"baryon"$, $cal(O)_"lepton"$ by (N4). For each non-abelian factor, (N2) fixes the field direction and (N3) the ghost and coupling directions. For $B$, (N2), (N3) and (N4) fix one direction each.]

*Existence.* (N4) removes 11 monomials but only three closed directions. The image of the 64 allowed monomials has dimension 64, while $bop("span"{Q_alpha})$ has dimension 72. A solution therefore exists only if $X_1$ satisfies eight linear conditions. They were located numerically, and all eight hold:
+ Five conditions live in the source classes $gh gh tilde(srcA)$, $gh gh tilde(srcA) A$, $gh gh gh srcC$, $gh gh srcH^dagger H$ and $gh gh H^dagger srcH$. At one loop $X_1$ vanishes in these classes. A one-loop diagram of $Dbrk_1 dot Gamma$ with an external $tilde(srcA)$, $srcC$ or $srcH$ would contain the source vertex attached to the fermion loop of $Dbrk_1$ through a single line, or it needs two loops. With $X_1 = 0$ there, five of the eight conditions hold.
+ The remaining three live only in the Yukawa classes $gh overline(psi) psi H$. In the coordinates of @tab:Nbasis, the component of $X_1$ along $N^(O_k)_(Y B) - N^(O_k^dagger)_(Y B)$, i.e. along $i y_k gy integral gh^B (O_k - O_k^dagger)$, must vanish for $k = u, d, e$. Only $Q^(f B)_R$ and $Q^((B))_sigma$ produce these structures (@eq:btable-yuk-B). At one loop $gh^B$ has no vertex with two internal lines, so in $X_1$ it sits on the vertex $Dbrk_1$. In the class $gh overline(psi) psi H$ both fermion legs of that vertex are then internal, and these graphs vanish under $Top(0)$. @sec:g1_cFFS establishes that $X_1$ has no component along $gh^B$ in this class, so the three conditions hold. For comparison, the non-vanishing contributions to $Dbrk dot Gamma_(gh overline(psi) psi H)$ in @Belusca-Maito:2020ala (Eq.~(6.26)) come from graphs with an internal ghost line and are proportional to the adjoint Casimir; such graphs do not exist for $gh^B$.
Hence (N1)--(N4) admit a solution, and @sec:st-solution constructs it. The three Yukawa conditions are the only place where the consistency of (N4) depends on the value of $X_1$, not on topology alone.

#general[The abelian ghost equation $delta Gamma \/ delta gh^B = delta S_0 \/ delta gh^B$ holds exactly in the regularised theory as long as the counterterms contain no $gh^B$. It is consistent with the Slavnov--Taylor identity, since ${delta \/ delta gh^B, bop}$ is the local hypercharge Ward operator. It is not by itself a proof that a $gh^B$-free solution exists at every order, because $X_l$ depends on $gh^B$ through $Dbrk_1$. If at some order $l$ the breaking $X_l$ (@eq:m-Xl) has a component along $integral gh^B (O_k - O_k^dagger)$, (N4) must be relaxed to (N4$'$): drop $Q^((B))_sigma$, $Q^(Q B)_R$, $Q^(L B)_R$ and the $gh^B$ ghost-sector monomials, and keep $Q^(U B)_R$, $Q^(D B)_R$, $Q^(E B)_R$. (N1)--(N3) with (N4$'$) give uniqueness, and a solution exists whenever $X_l$ vanishes in the source classes. A non-zero coefficient of $Q^(U B)_R$, $Q^(D B)_R$ or $Q^(E B)_R$ is a finite shift of the hypercharges in #sop that is not proportional to $Y$.]
