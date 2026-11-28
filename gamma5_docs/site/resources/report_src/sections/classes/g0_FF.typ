#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Fermion two-point functions $overline(Psi) Psi$ <sec:g0_FF>

External legs $overline(Psi)_f (-p)$ and $Psi_f (p)$ for $f in {Q, U, D, L, E}$, all momenta incoming: $Psi_f$ enters with momentum $p$, $overline(Psi)_f$ leaves with momentum $p$ (incoming momentum $-p$). Colour and $"SU"(2)_L$ indices of the two legs are $alpha' i'$ and $alpha i$. The sum of the diagrams below is $i thin Gamma^((1))_(overline(Psi)_f Psi_f)(p)$, a matrix in spinor space whose tree-level value is $Gamma^((0))_(overline(Psi)_f Psi_f)(p) = sl(p)$. The superficial degree of divergence (@sec:integrals) is $dg(G) = 4 + 1 - 2 dot 2 = 1$ for both topologies.

*Diagrams.* Every vertex that can be attached to the fermion line, @eq:V-ffv and @eq:V-yuk, contains $P_(c_f) Psi_f$ and $overline(Psi)_f P_(overline(c)_f)$ only. There is no $overline(Psi) Psi A A$, $overline(Psi) Psi H H$ or four-fermion vertex, so the one-loop 1PI graphs are the two bubbles of @fig:g0_FF-diagrams.

#figure(
  diagrow(
    oneloop(
      ((legs: (("F", [$Psi_f (p)$#h(6mm)], 1),)), (legs: (("F", $overline(Psi)_f$, -1),))),
      (("V", $A^a$, 0), ("F", $Psi_f$, -1)), caption: [(a) gauge exchange]),
    oneloop(
      ((legs: (("F", [$Psi_f (p)$#h(6mm)], 1),)), (legs: (("F", $overline(Psi)_f$, -1),))),
      (("S", $H$, 0), ("F", $Psi_(f')$, -1)), caption: [(b) Higgs exchange]),
  ),
  caption: [Topologies of $Gamma^((1))_(overline(Psi)_f Psi_f)$. In (a) the internal fermion is $Psi_f$ itself, in (b) it is a Yukawa partner $Psi_(f')$; the scalar line is $H$ or $H^dagger$ according to the vertices.],
) <fig:g0_FF-diagrams>

#figure(
  table(
    columns: 5,
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [symmetry factor], [sign], [group factor]),
    [(a)], [$Psi_f$, $A^a$ of every factor that acts on $f$: $G, W, B$ ($Q$); $G, B$ ($U$, $D$); $W, B$ ($L$); $B$ ($E$)], [$1$], [$+1$], [$sum_a Tg^a_f Tg^a_f = C_f bb(1)$],
    [(b), $f = Q$], [$Psi_U$, $H$; vertices $overline(U) Q H$, $overline(Q) U H^dagger$], [$1$], [$+1$], [$y_u^2 sum_j epsilon_(i' j) epsilon_(i j) = y_u^2 delta_(i' i)$],
    [(b), $f = Q$], [$Psi_D$, $H$; vertices $overline(D) Q H^dagger$, $overline(Q) D H$], [$1$], [$+1$], [$y_d^2 sum_j delta_(i' j) delta_(i j) = y_d^2 delta_(i' i)$],
    [(b), $f = U$], [$Psi_Q$, $H$; vertices $overline(Q) U H^dagger$, $overline(U) Q H$], [$1$], [$+1$], [$y_u^2 sum_(i j) epsilon_(i j) epsilon_(i j) = 2 y_u^2$],
    [(b), $f = D$], [$Psi_Q$, $H$; vertices $overline(Q) D H$, $overline(D) Q H^dagger$], [$1$], [$+1$], [$y_d^2 sum_(i j) delta_(i j) delta_(i j) = 2 y_d^2$],
    [(b), $f = L$], [$Psi_E$, $H$; vertices $overline(E) L H^dagger$, $overline(L) E H$], [$1$], [$+1$], [$y_e^2 delta_(i' i)$],
    [(b), $f = E$], [$Psi_L$, $H$; vertices $overline(L) E H$, $overline(E) L H^dagger$], [$1$], [$+1$], [$2 y_e^2$],
  ),
  caption: [Field assignments of @fig:g0_FF-diagrams: $10$ gauge-exchange and $6$ Higgs-exchange diagrams at multiplet level, counting each gauge factor separately ($16$ in all). There is no closed loop, so no sign; colour is conserved along the fermion line and gives $delta^(alpha' alpha)$ in (b). The first vertex listed absorbs the incoming $Psi_f$.],
) <tab:g0_FF-diagrams>

*Integrands.* The loop momentum $k$ flows along the internal fermion line, in the direction of the arrow; the boson line carries $p - k$. Writing the vertex at the outgoing end on the left,
$ G_"(a)" = intk (i gbar^mu P_(c_f) Tg^a_f) (i sl(k))/(k^2) (i gbar^nu P_(c_f) Tg^a_f) (-i g_(mu nu))/((k - p)^2) , $ <eq:g0_FF-int-a>
summed over $a$. For (b) with $f = Q$ and partner $D$ the vertices are $overline(Q)^(alpha' i') D_beta H_j$ ($-i y_d delta_(i' j) PR$) and $overline(D)^beta Q_(alpha i) H^dagger_(j')$ ($-i y_d delta_(i j') PL$), with the propagator $i delta_(j j') \/ ((k-p)^2 - m^2)$:
$ G_"(b)" = intk (-i y_d PR) (i sl(k))/(k^2) (-i y_d PL) (i)/((k - p)^2 - m^2) thin delta_(i' i) delta^(alpha' alpha) . $ <eq:g0_FF-int-b>
The other assignments of @tab:g0_FF-diagrams differ only in the coupling, in the $"SU"(2)_L$ factor and in $PL <-> PR$; in every case the projector to the right of $sl(k)$ is $P_(c_f)$ and the one to its left is $P_(overline(c)_f)$.

*$Rop$-operation.* At one loop $Rbar G = G$. _Dirac algebra._ By @eq:phys-prop, $P_(c_f) sl(k) P_(overline(c)_f) = slb(k) P_(overline(c)_f)$, and $gbar^mu g_(mu nu) gbar^nu = gbar^mu gbar_mu$. With @eq:sandwich and @eq:contr2,
$ gbar^mu P_(c_f) sl(k) gbar_mu P_(c_f) = gbar^mu slb(k) gbar_mu P_(c_f) = -2 slb(k) P_(c_f), quad quad
  PR sl(k) PL = slb(k) PL . $ <eq:g0_FF-dirac>
The prefactors are $i dot i dot i dot (-i) = -1$ in (a) and $(-i)(i)(-i)(i) = +1$ in (b), so
$ G_"(a)" = 2 C_f intk (slb(k) P_(c_f))/(k^2 (k - p)^2), quad quad
  G_"(b)" = y_d^2 intk (slb(k) P_(c_f))/(k^2 ((k - p)^2 - m^2)) thin delta_(i' i) delta^(alpha' alpha) . $ <eq:g0_FF-num>
The numerator contains only $gbar$ and $bar(k)$: no hatted object appears in it, and the Dirac algebra is four-dimensional.

_Expansion._ With @eq:uv-expansion, $1\/k^2 = 1\/(k^2 - MUV^2) + O(MUV^2)$ and $1\/((k-p)^2 - m^2) = 1\/(k^2 - MUV^2) + 2 k dot p \/ (k^2 - MUV^2)^2 + O(2)$, where $O(2)$ collects the terms of degree two in $p$, $m$ and $MUV$. Since $dg(G) = 1$, $Top(1)$ keeps the terms of degree $0$ and $1$. The degree-$0$ term $slb(k) \/ (k^2 - MUV^2)^2$ is odd in $k$ and vanishes; $MUV^2$ and $m^2$ are of degree two and do not enter. For a numerator $slb(k)$,
$ Top(1) intk slb(k)/(k^2 ((k - p)^2 - m^2)) = intk (2 slb(k) thin k dot p)/((k^2 - MUV^2)^3) = 2/d thin slb(p) intk (k^2)/((k^2 - MUV^2)^3) = ilf 1/ep (1/2 slb(p)) + O(ep^0) , $ <eq:g0_FF-T1>
by @eq:tensor-red-split ($bar(k)^mu k^nu -> k^2 gmb^(mu nu) \/ d$, so $slb(k) thin k dot p -> k^2 slb(p) \/ d$; the component $hat(k) dot hat(p)$ of $k dot p$ is paired with $bar(k)$ and gives zero) and @eq:logpole. The factor $2\/d$ comes from an unsplit #box[$d$-dimensional] pairing and is $1\/2 + O(ep)$; no trace $gmh^mu_(space mu)$ occurs. The same expansion holds for (a) with $m = 0$. Hence
$ N_0^"(a)" = C_f thin slb(p) P_(c_f), quad quad N_0^"(b)" = cal(Y)_f thin slb(p) P_(c_f), quad quad N_1^"(a)" = N_1^"(b)" = 0 , $ <eq:g0_FF-N>
where $cal(Y)_f$ is one half of the Yukawa group factor of @tab:g0_FF-diagrams, summed over the partners. These are the Yukawa sums $cal(Y)_f$ of @sec:st-solution:
$ cal(Y)_Q = (y_u^2 + y_d^2)/2, quad cal(Y)_U = y_u^2, quad cal(Y)_D = y_d^2, quad cal(Y)_L = y_e^2/2, quad cal(Y)_E = y_e^2 . $ <eq:g0_FF-Y>
In matrix form, $cal(Y)_f bb(1) = 1/2 sum_(f') sum_j bold(Y)^j_(f f') (bold(Y)^j_(f f'))^dagger$. Here the Feynman rule of the Yukawa vertex $overline(Psi)_f Psi_(f') phi_j$ of @eq:V-yuk is written $-i thin bold(Y)^j_(f f') P_(overline(c)_f)$: $bold(Y)^j_(f f')$ is the coupling matrix from the colour and $"SU"(2)_L$ indices of $f'$ to those of $f$, and $phi_j in {H_j, H^dagger_j}$ is the scalar field operator of that vertex, with $"SU"(2)_L$ index $j$. The same notation is used in @sec:g0_FFS. Example: $(bold(Y)^j_(Q D))^(alpha i)_(space beta) = y_d delta_(i j) delta^alpha_beta$ for $overline(Q) D H_j$, so that the $y_d$ part of $cal(Y)_Q$ is $1/2 y_d^2 sum_j delta_(i j) delta_(i' j) = 1/2 y_d^2 delta_(i i')$.

_No evanescent pole, no sterile leg._ The complete one-loop result has the form $P_(overline(c)_f) [dots] P_(c_f)$, with only barred matrices inside. A structure linear in $p$ with this sandwich can only be $slb(p) P_(c_f)$: the alternative $slh(p)$ vanishes between the projectors, $P_(overline(c)_f) ghat^mu P_(c_f) = 0$ (@eq:proj). Hence there is no evanescent pole. For the same reason every component with an external sterile partner vanishes: $Gamma^((1))_(overline(chi)_f chi_f) = Gamma^((1))_(overline(chi)_f psi_f) = Gamma^((1))_(overline(psi)_f chi_f) = 0$. In diagrams, an external $chi_f = P_(overline(c)_f) Psi_f$ or $overline(chi)_f = overline(Psi)_f P_(c_f)$ meets a vertex that contains $P_(c_f) Psi_f$ or $overline(Psi)_f P_(overline(c)_f)$, and $P_(c_f) P_(overline(c)_f) = 0$: there is no one-loop diagram. The evanescent part $slh(p)$ of the tree-level kinetic term, which couples $psi_f$ and $chi_f$, receives no one-loop correction.

*Result.*
#result[
  $ Gamma^((1))_(overline(Psi)_f Psi_f)(p) |_"pole" = lf 1/ep N_0, quad quad N_0 = N_0^"(a)" + N_0^"(b)" = (C_f + cal(Y)_f) thin slb(p) P_(c_f) , $ <eq:g0_FF-pole>
  with $C_f bb(1) = sum_a Tg^a_f Tg^a_f$ (@eq:invariants) and $cal(Y)_f$ of @eq:g0_FF-Y. The counterterm vertex function is $-lf (C_f + cal(Y)_f) slb(p) P_(c_f) \/ ep$. Since $integral dif^d x thin overline(Psi)_f thin i gbar^mu partial_mu P_(c_f) Psi_f = integral dif^d x thin overline(psi)_f thin i gbar^mu partial_mu psi_f$ has the vertex function $slb(p) P_(c_f)$,
  $ delta S^((1)) supset - lf 1/ep (C_f + cal(Y)_f) integral dif^d x thin overline(psi)_f thin i gbar^mu partial_mu psi_f = - lf 1/ep (C_f + cal(Y)_f) thin Q^f_K , $ <eq:g0_FF-ct>
  where $Q^f_K$ of @tab:Qbasis is continued to $d$ dimensions with $gbar^mu$ and $bar(partial)_mu$. No evanescent operator appears. The coefficient of $Q^f_K$ is the fermion field-renormalisation constant $delta Z_f$ of @sec:st-solution ($psi_(f 0) = Z_f^(1\/2) psi_f$, $Z_f = 1 + delta Z_f$; in the assembled counterterm the term $1/2 delta Z_f thin cal(O)^f_psi$ contains $delta Z_f thin Q^f_K$). In Feynman gauge
  $ delta Z_f = - lf 1/ep (C_f + cal(Y)_f) . $ <eq:g0_FF-dZ>
]

#figure(
  table(
    columns: 4,
    align: (center, left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$f$], [$C_f$], [$cal(Y)_f$], [$16 pi^2 ep thin delta Z_f$]),
    [$Q$], [$4/3 gs^2 + 3/4 gw^2 + 1/36 gy^2$], [$(y_u^2 + y_d^2)\/2$], [$-4/3 gs^2 - 3/4 gw^2 - 1/36 gy^2 - 1/2 y_u^2 - 1/2 y_d^2$],
    [$U$], [$4/3 gs^2 + 4/9 gy^2$], [$y_u^2$], [$-4/3 gs^2 - 4/9 gy^2 - y_u^2$],
    [$D$], [$4/3 gs^2 + 1/9 gy^2$], [$y_d^2$], [$-4/3 gs^2 - 1/9 gy^2 - y_d^2$],
    [$L$], [$3/4 gw^2 + 1/4 gy^2$], [$y_e^2 \/ 2$], [$-3/4 gw^2 - 1/4 gy^2 - 1/2 y_e^2$],
    [$E$], [$gy^2$], [$y_e^2$], [$-gy^2 - y_e^2$],
  ),
  caption: [Standard-Model values of @eq:g0_FF-pole and @eq:g0_FF-dZ. $C_f$ from @tab:invariants: per gauge factor $gs^2 C_F$ with $C_F = 4\/3$ for quarks, $gw^2 dot 3\/4$ for doublets, $gy^2 Y_f^2$.],
) <tab:g0_FF-values>

#check[*Auxiliary mass and Higgs mass.* $MUV$ and $m^2$ enter $Top(1)$ only at degree two, where they multiply an odd integrand or exceed $dg(G) = 1$; $N_0$ is independent of both, as required by dimension (a mass term $overline(psi)_f psi_f$ is forbidden by chirality).]

#check[*Anticommuting limit and literature.* The numerators @eq:g0_FF-num contain no hatted object, so @eq:g0_FF-pole coincides with the result for an anticommuting $g5$ and with the BMHV result of @Belusca-Maito:2020ala and @Lazopoulos:2026two: the gauge part $C_f$ is the Feynman-gauge Casimir term of a chiral fermion, and $cal(Y)_f$ is one half of the Yukawa matrix product $sum_(f') sum_j bold(Y)^j_(f f') (bold(Y)^j_(f f'))^dagger$. For $gs$ alone and a Dirac quark ($PL + PR$), $delta Z_f = - lf 1/ep thin 4/3 thin gs^2$ is the QCD quark field-renormalisation constant in Feynman gauge.]

#check[*Relation to other classes.* The poles of $(Dbrk_1 dot Gamma)^((1))$ are evanescent (@sec:m-oneloop), so the four-dimensional part of $delta S^((1))$ is $bop$-closed. In the class $gh^a overline(psi)_f psi_f$ the structure $N^(f r)_F$ is produced by $bop Q^f_K = sum_(r in f) N^(f r)_F$, $bop Q^(f r)_V = -N^(f r)_F$ (@eq:btable-fermion), $Q^(f r)_R$ and $Q^((r))_(rho partial)$. For $r = B$ the last two have no one-loop pole (the ghost $gh^B$ has no vertex with two internal lines, @sec:g0_ghost_sources), so the pole coefficient of $Q^(f B)_V$ found in @sec:g0_FFV must equal that of $Q^f_K$ in @eq:g0_FF-ct, Yukawa part included. For $r = G, W$ the same relation holds after adding the poles of $Q^(f r)_R$ and $Q^((r))_(rho partial)$ from @sec:g0_ghost_sources. The Yukawa part $cal(Y)_f$ also enters @sec:g0_FFS through the combination $y_k Q^+_(y_k)$ of $cal(O)^f_psi$ (@eq:closed). These relations are among the cross-class conditions of @sec:st-solution.]

#general[The absence of an evanescent pole and of sterile legs holds at every order for the ghost-number-$0$ function $Gamma_(overline(Psi)_f Psi_f)$. The external legs always attach to vertices containing $P_(c_f) Psi_f$, $overline(Psi)_f P_(overline(c)_f)$ (the counterterms $delta S^((L))$ and $Delta S^((L))$ contain no fermion bilinear with $chi_f$), so the function is $P_(overline(c)_f) [dots] P_(c_f)$. Its pole, of dimension one and linear in $p$, is proportional to $slb(p) P_(c_f)$ because $P_(overline(c)_f) slh(p) P_(c_f) = 0$. From two loops on, $Rbar$ subtracts the one-loop counterterms of @sec:g0_FF and @sec:g0_FFV, and the insertions of $Delta S^((1))$ and of the evanescent $E^((r))_A$, $E_H$ in the internal lines give $eta$-dependent finite parts that are kept as described in @sec:m-beyond. From three loops on, in bookkeeping B, the ghost-number-0 subgraphs with $eta$-dependent poles are subtracted with $delta S^((L))_eta$ and the remainder $E^((L))$ is inserted as a vertex (@eq:m-LPl). These change the coefficient, not the structure.]
