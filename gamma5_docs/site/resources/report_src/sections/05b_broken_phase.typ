#import "../macros.typ": *

// local shorthands (this file only)
#let litw(body) = [_(literature: #body)_]
#let der = [_(derived here)_]
#let expd = [_(expected, not derived)_]
#let Hv = $H_v$
#let Hb = $H_"bg"$
#let qb = $q_"bg"$
#let Ag = $A^gamma$
#let Yt = $tilde(Y)$
#let Mk = $M_"kin"$

== Subtleties in going to the broken phase <sec:broken>

@sec:Sren-broken gives the complete one-loop renormalised action of the broken phase in the linear gauge @eq:Lgf. It needs no loop input beyond @sec:classes. This subsection treats what that action leaves open: fermion propagators with sterile partners, the option of another gauge fixing, the tadpole convention, and the passage to the complex-mass scheme. Statements are marked #der, #litw[source] or #expd. The fields $Hv$, $phi$, the mixing angle, the mass matrix @eq:bp-masses and the fermion propagator @eq:bp-prop are those of @sec:bl-fields and @sec:bl-tree. Algebraic statements marked "checked" were verified with explicit matrices (@sec:supporting).

=== Why the shift needs no further loop input <sec:bp-shift>

@eq:bl-def rests on the three facts listed in @sec:bl-statement. Three points behind them are spelled out here #der.
+ *The breaking does not involve $H$.* $Dbrk_1$ (@eq:Delta1) contains only $gh$, $Psi$, $overline(Psi)$. The shift leaves @eq:ST, @eq:b-sources and @eq:linear-ids in the same form in $(phi, srcH)$ as in $(H, srcH)$, and $ST(cal(F))[Hv + phi] = ST(cal(F)[Hv + dot])[phi]$ for every functional $cal(F)$. Hence $Delta S^((1))_"ST" [Hv + phi]$ restores the identity in the broken phase if $Delta S^((1))_"ST"$ restores it at $H = 0$.
+ *Pole parts and $X_1$ are polynomials in $v$.* $Kop Top(dg(G))$ counts every mass and every vertex proportional to $v$ as a small quantity of @eq:uv-expansion. The pole part of a broken-phase graph is therefore the unbroken-phase pole part with external $H$ legs set to $Hv$ at zero momentum. The same holds for $X_1$, a residue of one-loop integrals (@eq:m-PG). This is why the coefficients of @tab:sts-r, computed at $H = 0$ with massless propagators, are the complete one-loop input.
+ *No $m^2$ in the one-loop $r_alpha$.* A one-loop graph of $Dbrk_1 dot Gamma$ is a single cycle containing $Dbrk_1$ and its internal sterile line. If the ghost of $Dbrk_1$ is external, the cycle is a loop of massless fermions. Higgs propagators occur only in classes with external fermions, whose local part has dimensionless coefficients. Hence $X_1$ does not depend on $m^2$, and $r^((r))_M = 0$: $N^((r))_1 = integral gh^a partial dot A^a$, the only image of $Q^((r))_M$, would need a coefficient of dimension two. The $H$-dependent coefficients are purely Yukawa (@tab:sts-r):
$ & r^((W))_J = r^((B))_J = Yt/3, quad r_(W W) = -Yt, quad r_(B B) = -(13/3 y_u^2 + 7/3 y_d^2 + 3 y_e^2), quad r_(W B) = -(y_u^2 + 2 y_d^2), \
  & r_(G G) = -1/2 (y_u^2 + y_d^2), quad r^((W))_D = r^((B))_D = r^((W))_sigma = 0 . $ <eq:bp-r1>
Beyond one loop $r^((r))_D = 0$ is expected from CP invariance of the regularised one-family theory with real $y_k$ #expd.

=== Fermion masses with sterile partners <sec:bp-fermion>

*Propagator* #der. After the shift @eq:LYuk gives $-m_u (overline(Psi)^1_Q PR Psi_U + overline(Psi)_U PL Psi^1_Q)$. It couples only $u_upright(L) = PL Psi^1_Q$ and $u_upright(R) = PR Psi_U$. With the physical and sterile Dirac fields
$ u = PL Psi^1_Q + PR Psi_U, quad quad u' = PR Psi^1_Q + PL Psi_U = chi^1_Q + chi_U , $ <eq:bp-uu>
@eq:Lferm-weyl gives the free Lagrangian of @eq:bl-L2F, $overline(u) (i slb(partial) - m) u + overline(u)' i slb(partial) u' + overline(u) i slh(partial) u' + overline(u)' i slh(partial) u$ ($m = m_u$; the same for $d$ and $e$). Its kinetic matrix $Mk (p)$ and the propagator $i Mk^(-1)(p)$ are @eq:bp-prop.

The physical block $S_"phys" (p) = (Mk^(-1)(p))_(11)$ is not the standard massive propagator. Its denominator $(p^2)^2 - m^2 bar(p)^2$ depends on $bar(p)^2$ and $p^2$ separately and is not $d$-dimensionally Lorentz covariant. For $hat(p) -> 0$ it vanishes at $bar(p)^2 = m^2$ and at $bar(p)^2 = 0$; the second zero is the sterile mode, with zero residue in the physical block. Graphs of ghost number $0$ contain only the physical block, because every vertex except $Dbrk_1$ contains only physical components. External fermions carry four-dimensional momenta, where the physical block is standard.

*(a) The $Rop$-operation is unchanged* #der. The mass term can be treated as a two-point vertex, $-i m PR$ for $overline(Psi)^1_Q Psi_U$ and $-i m PL$ for $overline(Psi)_U Psi^1_Q$, on the massless $d$-dimensional propagators $i sl(p) \/ p^2$. The geometric series of insertions reproduces $Mk^(-1)$, and $Top(dg(G))$ keeps finitely many insertions. Pole parts and $N_1$, hence $delta S^((1))$ and $X_1$, are polynomials in $m$ and equal the shifted unbroken-phase ones. This is point 2 of @sec:bp-shift for the fermion lines.

*(b) Finite parts are non-standard.* Amplitudes need integrals with the denominator of @eq:bp-prop. At one loop they can be traded for standard ones #der: the difference @eq:bp-AminusS between $S_"phys"$ and the standard propagator $S(p) = (slb(p) + m)\/(p^2 - m^2)$ is proportional to $hat(p)^2$, so a one-loop graph equals the same graph with $S$ plus the local term $ilf N_1$ of the difference, @eq:bp-trade (@sec:bl-massless). Beyond one loop @eq:bp-trade is not derived: the $hat(k)^2$ of one loop can multiply a divergent subloop.

*(c) Alternatives to sterile partners.* _Sterile Yukawa term._ Add $S_(y') = - sum_k y'_k integral dif^d x (O'_k + O'^dagger_k)$, where $O'_k$ is @eq:Ok with $PR -> PL$, e.g. $O'_u = overline(Q) tilde(H) PL U = overline(chi)_Q tilde(H) chi_U$. For $y'_k = y_k$ the propagator is the standard $i (sl(p) + m)\/(p^2 - m^2)$ in the $(Psi^1_Q, Psi_U)$ basis (checked). Since the sterile fields do not transform but $H$ does, $sop S_(y')$ adds to the breaking a term of ghost number $1$ and dimension $5$ without any hatted object #der. It is not evanescent: the four-dimensional limit is not BRST invariant at tree level, and the sterile fields become four-dimensional fermions coupled to $H$ (a $Psi$ loop in the Higgs self-energy gives the unprojected trace, twice the physical one). _Dirac pairing._ $u_upright(L)$ and $u_upright(R)$ are placed in one Dirac field with a $d$-dimensional kinetic term, as in @Kuhler:2026sm (where only $nu_upright(L)$ keeps a sterile partner) and @Cornella:2022hkc @OlgosoRuiz:2024dzq. The massive propagator is standard @Fuentes-Martin:2025meh, but the evanescent term $overline(u)_upright(R) slh(partial) u_upright(L) + "h.c."$ couples different $"SU"(2)_L times "U"(1)_Y$ representations and breaks the rigid symmetry; hypercharge-violating Green's functions must then be computed @Kuhler:2026sm.

#figure(
  {
    set par(justify: false)
    set text(size: 9pt)
    table(
    columns: (1.1fr, 0.8fr, 0.8fr, 1.45fr, 1.9fr),
    align: (left, left, left, left, left),
    inset: (x: 4pt, y: 3.5pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { (bottom: 0.2pt + luma(200)) },
    table.header([option], [rigid symmetry], [standard propagator], [tree breaking $ST(S_0)$], [cost]),
    [sterile partners (this report)], [exact], [no], [$Dbrk_1$: evanescent, dimension 5], [none for $delta S^((1))$, $X_1$; finite parts need @eq:bp-prop, or @eq:bp-trade at one loop],
    [sterile Yukawa], [exact], [yes], [$Dbrk_1 + sop S_(y')$: not evanescent], [not admissible: the limit $d -> 4$ is not the Standard Model],
    [Dirac pairing], [broken], [yes], [evanescent terms $overline(f)_upright(R) slh(partial) f_upright(L) + "h.c."$ with left and right generators, dimension 5], [non-invariant bases, hypercharge-violating classes; unbroken phase redone],
  )},
  caption: [Three ways to give the fermions a mass. First two rows derived here; third row from @Kuhler:2026sm.],
) <tab:bp-options>

The sterile-Yukawa option is excluded. Between the other two the choice is between a non-standard propagator in finite parts and a larger basis without rigid invariance. This report does not decide it.

*Status.* Derived and checked: @eq:bp-prop and its limits, @eq:bp-AminusS, the unchanged $Rop$-operation, the standard propagator of the sterile-Yukawa option. @eq:bp-trade is derived and was applied to the photon, $Z$--photon and gluon two-point functions at zero momentum, where an independent evaluation through the Slavnov--Taylor identity gives the same result (@sec:bl-massless). Not done: broken-phase finite parts at non-zero external momenta with the propagator @eq:bp-prop.

=== Choosing another gauge <sec:bp-gauge>

*Linear gauge.* @eq:Lgf does not depend on $H$. In the broken phase the tree-level gauge--Goldstone mixing @eq:bl-L2mix therefore stays in the quadratic action. It is part of the propagators @eq:bl-prop: mixed gauge--Goldstone propagators and double poles at $k^2 = 0$ in the longitudinal and Goldstone parts. This is a consistent gauge, and @sec:Sren-broken is complete in it.

*$R_xi$ gauge (optional).* The $R_xi$ gauge removes the mixing. It is a change of gauge, not a requirement of the broken phase. In the notation of @eq:gf-exact it is the BRST-exact term
$ cal(L)^(R_xi)_"gf" = sop [agh^a (F^a + xi/2 b^a)], quad F^a = partial^mu A^a_mu + i xi (H^dagger Tg^a_H Hv - Hv^dagger Tg^a_H H) . $ <eq:bp-Rxi>
With @eq:bl-rot, $F^Z = partial dot Z - xi M_Z G^0$, $F^plus.minus = partial dot W^plus.minus minus.plus i xi M_W G^plus.minus$, $F^gamma = partial dot Ag$, $F^A = partial dot G^A$. Eliminating $b^a = -F^a \/ xi$ gives, up to a total derivative, #der (checked)
$ - 1/(2 xi) sum_a F^a F^a = - 1/(2 xi) sum_a (partial dot A^a)^2 - cal(L)^(V G)_2 - xi (M_W^2 G^+ G^- + 1/2 M_Z^2 (G^0)^2) , $ <eq:bp-tree>
so the cross term cancels @eq:bl-L2mix and the Goldstone fields acquire the masses $xi M_W^2$, $xi M_Z^2$. The ghost part of @eq:bp-Rxi is $partial^mu agh^a (D_mu gh)^a - xi thin agh^a gh^b (H^dagger Tg^b_H Tg^a_H Hv + Hv^dagger Tg^a_H Tg^b_H H)$. At $H = Hv$ the second term is the ghost mass matrix $xi cal(M)^2_(a b)$ of @eq:bp-masses (checked); its $phi$-linear part is a ghost--ghost--scalar vertex proportional to $xi v$. This is a different action from @eq:bl-def. Using it requires input that this report does not provide:
+ *Modified ghost equation* #der. The antighost couples to $sop H = delta Gamma \/ delta srcH^dagger$, and the ghost equation of @eq:linear-ids becomes
  $ (delta Gamma)/(delta agh^a) + partial_mu (delta Gamma)/(delta srcA^a_mu) + i xi [ (delta_R Gamma)/(delta srcH) Tg^a_H Hv - Hv^dagger Tg^a_H (delta Gamma)/(delta srcH^dagger) ] = 0 . $ <eq:bp-ghost>
  $Gamma$ depends on $agh$ through $tilde(srcA)$ and $tilde(srcH) = srcH + i xi Tg^a_H Hv agh^a$, and $delta Gamma \/ delta b^a = F^a + xi b^a$.
+ *Rigid invariance.* Since $Hv$ is fixed, @eq:bp-Rxi breaks rigid $"SU"(2)_L times "U"(1)_Y$ to $"U"(1)_"em"$, while the bases of @sec:bases assume rigid invariance. The remedy (algebraic renormalisation of the electroweak Standard Model @Kraus:1997bi and its background-field formulation @Grassi:1999nb) replaces $Hv$ by an external, non-propagating doublet $Hb$ with the quantum numbers of $H$ and a Grassmann-odd partner $qb$ of ghost number $1$, $sop Hb = qb$, $sop qb = 0$. Rigid invariance is then exact, and $sop[agh^a (F^a + xi b^a \/ 2)]$ contains three new vertices #der:
  $ & "(a)" quad i xi b^a (H^dagger Tg^a_H Hb - Hb^dagger Tg^a_H H), quad quad
    "(c)" quad - i xi thin agh^a (H^dagger Tg^a_H qb - qb^dagger Tg^a_H H), \
    & "(b)" quad - xi thin agh^a gh^b (H^dagger Tg^b_H Tg^a_H Hb + Hb^dagger Tg^a_H Tg^b_H H) . $ <eq:bp-bgvertices>
  $ST$ acquires $integral (qb^dagger delta \/ delta Hb^dagger + delta \/ delta Hb thin qb)$, and the ghost equation is @eq:bp-ghost with $Hv -> Hb$ and $-i xi (H^dagger Tg^a_H qb - qb^dagger Tg^a_H H)$ on the right-hand side. $Hb = Hv$, $qb = 0$, $H = Hv + phi$ gives the $R_xi$ gauge. At $Hb = qb = 0$ all three vertices vanish, so every Green's function of this report is the $Hb = qb = 0$ part of this formulation #der.
+ *Additional Green's functions* (power counting with $dim Hb = dim qb = 1$). Ghost number $0$: $Hb H$, $Hb Hb$, $Hb H A$, $Hb Hb A$, $Hb H A A$, $Hb Hb A A$, the four-scalar classes with one to four $Hb$, $Hb overline(psi) psi$, $srcH gh Hb$, $qb overline(srcF) psi$, $srcC qb qb$. Ghost number $1$: the classes of @sec:Nbasis with $H -> Hb$, and the classes with one $qb$: $qb H$, $qb H A$, $qb H A A$, $qb H H H$, $qb overline(psi) psi$, $qb overline(psi) psi A$, $gh qb overline(srcF) psi$. The bases, the closed directions and the normalisation conditions must be extended accordingly.
+ *One-loop classes with new diagrams* #der. A one-loop 1PI graph of $Dbrk_1 dot Gamma$ is a cycle with $Dbrk_1$ on it and its sterile leg internal. If the ghost of $Dbrk_1$ is external, the cycle is a closed fermion loop with gauge and Yukawa vertices only. @eq:bp-bgvertices contains no fermion, so the bosonic classes ($gh A^n$, $gh H H A^n$, $gh H^4$) receive no new diagrams. If the ghost is internal, the cycle has the open fermion line of $Dbrk_1$ and a ghost segment from $Dbrk_1$ to an $agh$ leg. These are the classes with an internal ghost line, $gh overline(Psi) Psi$, $gh overline(Psi) Psi A$, $gh overline(Psi) Psi H$, $gh gh overline(srcF) Psi$. The segment can end at vertex (b) (external $gh$ and $Hb$, internal $H$ to a Yukawa vertex) or at vertex (c) (external $qb$). This gives new diagrams in $gh overline(psi) psi' Hb$, $qb overline(psi) psi'$, $qb overline(psi) psi' A$, $qb overline(psi) psi' H$ and $gh qb overline(srcF) psi'$, with $psi'$ the Yukawa partner of $psi$. Vertex (b) inside a ghost segment leaves $H$ and $Hb$ external, with dimension $>= 6$, and does not contribute to $X_1$. At $Hb = Hv$ these are the $xi v$-dependent parts of $gh overline(psi) psi$ and its relatives. They can require further finite restoring terms.
None of these Green's functions is computed in this report.

=== The vacuum value and tadpoles <sec:bp-tadpole>

*Convention of @sec:Sren-broken.* $v$ is the minimum of the tree-level potential in the $overline("MS")$ parameters, @eq:bl-vev, and tadpole graphs are kept in every Green's function, the two-point functions of @sec:bp-cms included. The tadpole counterterm is the term $-v thin delta T thin h$ of @eq:bl-dS-scalar, with $delta T$ of @eq:bl-dT; it cancels the pole of the one-loop tadpole. The finite tadpole $Gamma^((1))_h$ is not zero. $Delta S^((1))_"ST"$ gives no tadpole and no scalar two-point term (@tab:bl-ST). The one-loop tadpole itself, and hence the finite part of every tadpole insertion, is not computed here.

*Other conventions* _(literature: #cite(<Fleischer:1980ub>), #cite(<Denner:2016etu>))_. Keeping $v$ at the tree-level minimum in terms of the Lagrangian parameters, with explicit tadpoles, is the Fleischer--Jegerlehner treatment. The relations between $overline("MS")$ parameters and physical quantities are then gauge independent, at the price of large tadpole contributions (enhanced by $m_t^4 \/ M_h^2$ in the three-family Standard Model). If instead the renormalised tadpole is set to zero, $v$ is the minimum of the gauge-dependent effective potential, and the $overline("MS")$ masses depend on $xi$. The two conventions differ by a finite reparametrisation of $v$. With a non-anticommuting $g5$ the gauge independence of the first convention also needs the restored identity and a control of the $xi$-dependence #expd. This report works at $xi = 1$ only.

=== From $overline("MS")$ to the complex-mass scheme <sec:bp-cms>

*Finite renormalisation along closed directions* #der. Let $theta = (gs, gw, gy, lambda, m^2, y_k)$ be the $overline("MS")$ parameters and $theta_c$ the complex-mass-scheme ones, $theta = theta_c + delta theta(theta_c)$ with $delta theta = O(ħ)$ finite and complex @Denner:2005fg @Denner:2006ic. The regularised action is unchanged and re-expanded:
$ S_"ren" (theta) = S_"ren" (theta_c) + sum_i delta theta_i (partial S_0)/(partial theta_i) |_(theta_c) + delta S_Z + O(ħ^2) , $ <eq:bp-cms>
with $delta S_Z$ finite field and source renormalisations. In four dimensions $c_r partial S_0 \/ partial c_r = cal(O)^((r))_g$, $partial S_0 \/ partial m^2 = -Q_m$, $partial S_0 \/ partial lambda = -Q_lambda$, $partial S_0 \/ partial y_k = Q^+_(y_k)$, and $delta S_Z$ lies in the exact directions $cal(O)^f_psi$, $cal(O)_H$, $cal(O)^((r))_A$, $cal(O)^((r))_c$: all are in @eq:closed. So $bop$ of the change vanishes, and $X_1$ and $Delta S^((1))_"ST"$ are unchanged. The statement holds at every order: the bare theory is unchanged, $ST$ does not act on parameters, and $ST$ is invariant under $Phi -> z^(1\/2) Phi$, $Phi^* -> z^(-1\/2) Phi^*$ for every field $Phi$ with source $Phi^*$ (@sec:brst). The higher-order $r_alpha$ are the $overline("MS")$ ones expressed through $theta_c$, provided @eq:bp-cms uses the $d$-dimensional $partial S_0 \/ partial theta$ (its evanescent part matters at two loops). The algebra of @sec:bases is polynomial in the couplings, so the $r_alpha$ continue analytically to complex parameters.

*What must be computed.* The complex poles $mu_V^2$ ($V = W, Z$), $mu_h^2$ and $mu_f$ follow from the one-loop two-point functions of the broken phase: $overline("MS")$ finite parts with tadpoles and with $Delta S^((1))_"ST" [Hv + phi]$. With $Gamma^(mu nu)_(V V) = -[p^2 - M_V^2 + Sigma^T_(V V) (p^2)] (g^(mu nu) - p^mu p^nu \/ p^2) + dots$, one has $mu_V^2 = M_V^2 - Sigma^T_(V V) (mu_V^2)$ at one loop. The couplings follow from the masses ($c_W^2 = mu_W^2 \/ mu_Z^2$) and one vertex-defined coupling such as $e$. Reducing the renormalisation of $e$ to two-point functions uses the Ward identity, which holds only after $Delta S_"ST"$ is included #expd.

*Which two-point functions depend on $Delta S_"ST"$* #der. Under (N1)--(N3) the only two-point monomials of @tab:Qbasis with non-zero coefficient at $H = 0$ are $Q^((r))_L$ ($r^((r))_M = 0$). In the broken phase the shifted $H$-dependent monomials add the mass-like terms @eq:bl-M and the mixing @eq:bl-mix. @tab:bp-twopoint lists the consequences.

#figure(
  {
    set par(justify: false)
    set text(size: 9pt)
    table(
    columns: (auto, 1fr, 1.15fr),
    align: (left, left, left),
    inset: (x: 4pt, y: 3.5pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { (bottom: 0.2pt + luma(200)) },
    table.header([two-point function], [terms of $Delta S^((1))_"ST" [Hv + phi]$], [consequence]),
    [$overline(f) f$; $h h$, $G G$, $h G$], [none (N1)], [insensitive],
    [$V V'$ transverse], [mass-like constants @eq:bl-M], [$W$, $Z$ poles sensitive; $gamma Z$, $gamma gamma$, $g g$ at $p^2 = 0$: the constants cancel the fermion-loop mass terms (@eq:bl-massless); $partial Sigma^T \/ partial p^2$ insensitive],
    [$V V'$ longitudinal], [@eq:bl-L and the constants @eq:bl-M], [sensitive (unphysical)],
    [$V G$], [@eq:bl-mix], [sensitive (unphysical)],
    [$V h$; $srcH gh$], [none ($r^((r))_D = r^((W))_sigma = 0$)], [no finite restoring term at one loop],
  )},
  caption: [Two-point functions of the broken phase (linear gauge, $xi = 1$) and the terms of the shifted restoring counterterm that enter them. $V, V' in {W^plus.minus, Z, Ag, G^A}$ label gauge bosons, $G$ the Goldstone fields.],
) <tab:bp-twopoint>

With $Gamma^(mu nu)_(V V')$ normalised as above for each pair, a term $1/2 lf Delta M^2_(V V') A^V dot A^(V')$ of the action contributes $-lf Delta M^2_(V V')$ to $Sigma^T_(V V')$. With @eq:bl-M and @eq:bl-dMVV #der:
$ & Sigma^T_(W W) supset lf gw^2/2 (3 m_u^2 + 3 m_d^2 + m_e^2), quad Sigma^T_(g g) supset lf gs^2 (m_u^2 + m_d^2), \
  & Sigma^T_(V V') supset 2 lf sum_(f = u, d, e) N^c_f m_f^2 (v^V_f v^(V')_f + 3 a^V_f a^(V')_f) quad (V, V' in {Z, gamma}), \
  & Sigma^T_(gamma gamma) supset 2 e^2 lf sum_f N^c_f q_f^2 m_f^2 . $ <eq:bp-dM-1loop>
The $W$ and $Z$ poles therefore depend on $r_(W W)$, $r_(B B)$, $r_(W B)$. The photon and gluon terms of @eq:bp-dM-1loop cancel the finite mass terms of the massive fermion loops in the regularised theory: at one loop and zero momentum the fermion loops plus these terms give no $gamma gamma$, $gamma Z$ or $g g$ mass term (@eq:bl-massless, @sec:bl-massless). These are checks of the $r_alpha$.

*LSZ with $overline("MS")$ fields.* Residues are finite but not one, and complex for unstable particles. Near the pole the propagator is $i R \/ (p^2 - mu^2)$ with $R = [1 + partial Sigma^T \/ partial p^2 |_(mu^2)]^(-1)$. An $S$-matrix element is the amputated Green's function at $p_i^2 = mu_i^2$ times $R_i^(1\/2)$ per leg, with separate factors for the two fermion chiralities. For $Ag$ and $Z$ the residue is a matrix: at one loop the external photon receives the $Z$-leg amplitude with a coefficient proportional to $Sigma^T_(gamma Z) (0) \/ M_Z^2$, which no $overline("MS")$ condition removes (its fermion-loop mass part cancels against $Delta S^((1))_"ST"$, @eq:bl-massless; the bosonic loops remain). Tadpole graphs enter $mu^2$ and $R$. Poles and residues come from the complex self-energies at $mu^2$, or from a consistent expansion around it. No on-shell scheme and no real-part prescription is used. External fermions have four-dimensional momenta, where the physical block of @eq:bp-prop is standard, and by (N1) their residues come from loops and $delta S^((1))$ only #der.

=== What is needed beyond this report <sec:bp-needed>

The renormalised action of the broken phase in the linear gauge is not on this list: it is @sec:Sren-broken. Its photon, $Z$--photon and gluon two-point functions at zero momentum are checked in @sec:bl-massless. Still needed are:
+ The one-loop broken-phase two-point functions at complex momenta, with tadpoles and $Delta S^((1))_"ST" [Hv + phi]$, and the one-loop tadpole (only its $lambda^2$ pole is checked in @sec:bl-dS).
+ The local terms $N_1$ of @eq:bp-trade for each amplitude at non-zero momenta, or integrals with the denominator of @eq:bp-prop.
+ The vertex function that defines $e$, with $Delta S^((1))_"ST"$.
+ The LSZ residue matrices of the mixing systems.
+ A control of the $xi$-dependence (the report is at $xi = 1$) that establishes gauge independence of the $overline("MS")$ relations with a non-anticommuting $g5$.
+ A proof that $r^((r))_D = 0$ by CP beyond one loop, or its value.
+ Only if the $R_xi$ gauge is chosen: $Hb$, $qb$, the extended bases, closed directions and normalisation conditions, and the Green's functions listed in @sec:bp-gauge.
+ For two loops: $tilde(E)_H [Hv + phi]$, the $d$-dimensional form of @eq:bp-cms, and @eq:bp-prop in subloops.
