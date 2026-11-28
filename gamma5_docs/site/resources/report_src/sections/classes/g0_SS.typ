#import "../../macros.typ": *
#import "../../diagrams.typ": *

// local shorthands (this file only)
#let DM(n) = $(k^2 - MUV^2)^(#n)$
#let Yt = $tilde(Y)$

=== Higgs two-point function $H^dagger H$ <sec:g0_SS>

External legs $H_j (p)$ and $H^dagger_i (-p)$, all momenta incoming: the doublet line enters with momentum $p$ through $H_j$ and leaves through $H^dagger_i$. The sum of the diagrams below is $i thin Gamma^((1))_(H^dagger_i H_j) (p)$, with tree level $Gamma^((0))_(H^dagger_i H_j) = (p^2 - m^2) delta_(i j)$. Every diagram has superficial degree $dg(G) = 2$. The mass $m^2$ is kept in the Higgs propagator; the gauge factors $W$ and $B$ and the three Yukawa couplings contribute, the gluon does not ($Tg^A_H = 0$).

*Diagrams.* @fig:g0_SS-diagrams shows the four topologies, @tab:g0_SS-diagrams the field assignments.
#figure(
  diagrow(
    oneloop(
      ((legs: (("S", $H_j$, 1),)), (legs: (("S", $H^dagger_i$, -1),))),
      (("F", $psi'$, 1), ("F", $psi$, 1)), caption: [(a) fermion loop]),
    oneloop(
      ((legs: (("S", $H_j$, 1),)), (legs: (("S", $H^dagger_i$, -1),))),
      (("S", $H$, 1), ("V", $A^a$, 0)), caption: [(b) gauge--Higgs bubble]),
    oneloop(
      ((legs: (("S", $H_j$, 1), ("S", $H^dagger_i$, -1))),),
      (("V", $A^a$, 0),), caption: [(c) gauge tadpole]),
    oneloop(
      ((legs: (("S", $H_j$, 1), ("S", $H^dagger_i$, -1))),),
      (("S", $H$, 1),), caption: [(d) quartic tadpole]),
  ),
  caption: [Topologies of $Gamma^((1))_(H^dagger H)$. In (a) $psi$ enters the $H_j$ vertex and $psi'$ leaves it: $(psi, psi') = (D, Q), (E, L), (Q, U)$ for the Yukawa couplings $y_d, y_e, y_u$.],
) <fig:g0_SS-diagrams>

#figure(
  table(
    columns: 5,
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [symmetry factor], [sign], [group factor]),
    [$"(a)"_d$], [$Q$, $D$], [$1$], [$-1$ (fermion loop)], [$y_d^2 thin delta^alpha_alpha thin delta_(l i) delta_(l j) = 3 y_d^2 thin delta_(i j)$],
    [$"(a)"_e$], [$L$, $E$], [$1$], [$-1$ (fermion loop)], [$y_e^2 thin delta_(l i) delta_(l j) = y_e^2 thin delta_(i j)$],
    [$"(a)"_u$], [$Q$, $U$], [$1$], [$-1$ (fermion loop)], [$y_u^2 thin delta^alpha_alpha thin epsilon_(l i) epsilon_(l j) = 3 y_u^2 thin delta_(i j)$],
    [(b)], [$H$, $A^a$ ($a = I, B$)], [$1$], [$+1$], [$(Tg^a_H Tg^a_H)_(i j) = C_H delta_(i j)$],
    [(c)], [$A^a$ ($a = I, B$)], [$1\/2$], [$+1$], [${Tg^a_H, Tg^a_H}_(i j) = 2 C_H delta_(i j)$],
    [(d)], [$H$], [$1$], [$+1$], [$(delta_(i j) delta_(n l) + delta_(i l) delta_(n j)) delta_(l n) = 3 delta_(i j)$],
  ),
  caption: [Field assignments, symmetry factors, signs and group factors. $l$ and $n$ are internal $"SU"(2)_L$ indices, $alpha$ the colour index of the quark loop. The vertex @eq:V-ssss is the complete fourfold derivative of $-lambda (H^dagger H)^2$, so (d) has symmetry factor $1$ (an oriented charged loop).],
) <tab:g0_SS-diagrams>

The Yukawa factors are collected in
$ Yt = 3 y_u^2 + 3 y_d^2 + y_e^2 , $ <eq:g0_SS-Y>
the sum over the three fermion loops of the colour multiplicity times the squared coupling ($tilde(Y)$ of @sec:st-solution). For the Higgs doublet $C_H = 3/4 gw^2 + 1/4 gy^2$ (@tab:invariants).

*Integrands.* (a) For the $Q$--$D$ loop the $H_j$ vertex is $overline(Q) D H_j$ and the $H^dagger_i$ vertex is $overline(D) Q H^dagger_i$ (@eq:V-yuk). With $k$ on the $D$ line (from the $H^dagger$ vertex to the $H$ vertex) and $k + p$ on the $Q$ line,
$ G_"(a)" = - 3 delta_(i j) intk tr[ (-i y_d PR) (i sl(k))/(k^2) (-i y_d PL) (i (sl(k) + sl(p)))/((k+p)^2) ]
  = - 3 y_d^2 thin delta_(i j) integral_k (tr[slb(k) (slb(k) + slb(p)) PR])/(k^2 (k+p)^2) , $ <eq:g0_SS-int-a>
where @eq:sandwich gave $PR sl(k) PL (sl(k) + sl(p)) PR = slb(k) (slb(k) + slb(p)) PR$: the numerator is four-dimensional, the denominators are #box[$d$-dimensional] (@eq:phys-prop). The $L$--$E$ loop is the same with $3 y_d^2 -> y_e^2$. For the $Q$--$U$ loop the projectors are interchanged ($PL$ at the $H_j$ vertex $overline(U) Q H_j$) and the colour--isospin factor is $3 epsilon_(l i) epsilon_(l j) = 3 delta_(i j)$; the trace is again $tr[slb(k)(slb(k) + slb(p)) PL]$. Since $tr[slb(a) slb(b) g5] = 0$, all three give
$ G_"(a)" = - Yt thin delta_(i j) integral_k (2 thin bar(k) dot (bar(k) + bar(p)))/(k^2 (k+p)^2) . $ <eq:g0_SS-int-a2>
(b) With $k$ the gauge-boson momentum flowing into the $H_j$ vertex and $k + p$ on the internal Higgs line, both vertices @eq:V-ssv carry $i (2 p + k)^mu$:
$ G_"(b)" & = integral_k i (2p + k)^nu (Tg^a_H)_(i l) thin (i)/((k+p)^2 - m^2) thin i (2 p + k)^mu (Tg^a_H)_(l j) thin (-i g_(mu nu))/(k^2) \
  & = - C_H thin delta_(i j) integral_k ((2 p + k)^2)/(k^2 ((k + p)^2 - m^2)) . $ <eq:g0_SS-int-b>
(c) The seagull @eq:V-ssv with the two gauge legs contracted:
$ G_"(c)" = 1/2 integral_k i g^(mu nu) {Tg^a_H, Tg^a_H}_(i j) thin (-i g_(mu nu))/(k^2) = d thin C_H thin delta_(i j) integral_k 1/(k^2) . $ <eq:g0_SS-int-c>
(d) The quartic vertex @eq:V-ssss, with the legs $H^dagger_i H_j H^dagger_n H_l$ and $H_l$, $H^dagger_n$ contracted:
$ G_"(d)" = integral_k (-2 i lambda) (delta_(i j) delta_(n l) + delta_(i l) delta_(n j)) (i delta_(l n))/(k^2 - m^2) = 6 lambda thin delta_(i j) integral_k 1/(k^2 - m^2) . $ <eq:g0_SS-int-d>

*$Rop$-operation.* At one loop $Rbar G = G$ and $Z(G) = -Kop Top(2) G$ (@eq:m-Z1). With @eq:uv-expansion, $1/k^2 = 1/(k^2 - MUV^2) - MUV^2/DM(2) + dots$ and, to degree two,
$ Top(2) 1/(k^2 ((k+p)^2 - m^2)) = 1/DM(2) - (2 k dot p + p^2 - m^2 + 2 MUV^2)/DM(3) + (4 (k dot p)^2)/DM(4) , $ <eq:g0_SS-T2>
with $k dot p = bar(k) dot bar(p) + hat(k) dot hat(p)$ and $m = 0$ for (a). The vacuum poles needed are (@eq:logpole and the cases $b = a + 1$)
$ Kop intk 1/(k^2 - MUV^2) = ilf MUV^2/ep, quad Kop intk (k^2)/DM(2) = ilf (2 MUV^2)/ep, \
  Kop intk 1/DM(2) = Kop intk (k^2)/DM(3) = Kop intk ((k^2)^2)/DM(4) = ilf 1/ep . $
The tensor reductions are @eq:tensor-red-split and
$ bar(k)^2 (k dot p)^2 -> ((k^2)^2)/(d(d+2)) (4 p^2 + 2 bar(p)^2), quad
  (bar(k) dot bar(p))(k dot p) -> (k^2)/d thin bar(p)^2, quad (k dot p)^2 -> (k^2)/d thin p^2 , $ <eq:g0_SS-red>
with $p^2 = bar(p)^2 + hat(p)^2$; the #box[$d$-dimensional] metric of @eq:tensor-red produces the full $p^2$ from the products of $bar(k) dot bar(p)$ and $hat(k) dot hat(p)$ of the denominator.

(a) The numerator $2 thin bar(k)^2 + 2 thin bar(k) dot bar(p)$ times @eq:g0_SS-T2 gives, term by term after reduction and $Kop$ (in units of $ilf \/ ep$; factors $4\/d$ are taken at $ep = 0$ because the poles are simple):
$ underbrace(2 bar(k)^2 / DM(2), 4 MUV^2) quad
  underbrace(- 2 bar(k)^2 (p^2 + 2 MUV^2)/DM(3), -2 p^2 - 4 MUV^2) quad
  underbrace(8 bar(k)^2 (k dot p)^2 / DM(4), (4 p^2 + 2 bar(p)^2)\/3) quad
  underbrace(-4 (bar(k) dot bar(p))(k dot p) / DM(3), - bar(p)^2) $
and the terms odd in $k$ vanish. The sum is $-2/3 p^2 - 1/3 bar(p)^2 = -bar(p)^2 - 2/3 hat(p)^2$. With the prefactor $-Yt$ of @eq:g0_SS-int-a2,
$ Top(2) G_"(a)" = ilf thin 1/ep N_0^"(a)" + dots, quad quad N_0^"(a)" = Yt thin (bar(p)^2 + 2/3 hat(p)^2) thin delta_(i j) . $ <eq:g0_SS-N-a>
No hatted index is contracted with itself (the only hatted loop components are $hat(k) dot hat(p)$), so no $eta$ is produced.

(b) The numerator $k^2 + 4 k dot p + 4 p^2$ times @eq:g0_SS-T2 gives the non-vanishing terms
$ underbrace(k^2/DM(2), 2 MUV^2) quad
  underbrace(- k^2 (p^2 - m^2 + 2 MUV^2)/DM(3), -p^2 + m^2 - 2 MUV^2) quad
  underbrace(4 k^2 (k dot p)^2/DM(4), p^2) quad
  underbrace(- 8 (k dot p)^2 / DM(3), -2 p^2) quad
  underbrace(4 p^2 / DM(2), 4 p^2) $
in units of $ilf\/ep$, with sum $2 p^2 + m^2$:
$ Top(2) G_"(b)" = ilf thin 1/ep N_0^"(b)" + dots, quad quad N_0^"(b)" = - C_H (2 p^2 + m^2) thin delta_(i j) . $ <eq:g0_SS-N-b>
Here all vectors are #box[$d$-dimensional] and $p^2$ is the full square.

(c) $Top(2) 1\/k^2 = 1\/(k^2 - MUV^2) - MUV^2 \/ DM(2)$, whose pole is $ilf (MUV^2 - MUV^2)\/ep = 0$. By @eq:tadpole the combination is finite ($Gamma(-1 + ep) + Gamma(ep) = Gamma(ep) thin ep\/(ep - 1)$ up to the common factor), so the factor $d$ does not create a pole:
$ N_0^"(c)" = 0 . $ <eq:g0_SS-N-c>
(d) $Top(2) 1\/(k^2 - m^2) = 1\/(k^2 - MUV^2) + (m^2 - MUV^2)\/DM(2)$, with pole $ilf thin m^2 \/ ep$:
$ N_0^"(d)" = 6 lambda m^2 thin delta_(i j) . $ <eq:g0_SS-N-d>
In every topology the $MUV^2$ terms cancel, and $N_1 = 0$.

*Result.* With the translation of @sec:m-objects (after @eq:m-dS-per-graph: $p_j -> i partial$ on leg $j$), the vertex functions $p^2 delta_(i j)$, $bar(p)^2 delta_(i j)$, $hat(p)^2 delta_(i j)$ and $delta_(i j)$ belong to $integral dif^d x thin (partial_mu H)^dagger partial^mu H$, $integral dif^d x thin (bar(partial)_mu H)^dagger bar(partial)^mu H$, $E_H$ (@eq:evanescent) and $Q_m = integral dif^d x thin H^dagger H$; the two-point part of $tilde(E)_H$ (@eq:evanescent) has the vertex function $-hat(p)^2 delta_(i j)$. $Q_(S S)$ and $Q_m$ denote the #box[$d$-dimensional] continuations of @tab:Qbasis, as in $S_0$.
#result[
  $ Gamma^((1))_(H^dagger_i H_j) |_"pole" = lf 1/ep N_0, quad
    N_0 = [ Yt (bar(p)^2 + 2/3 hat(p)^2) - (2 p^2 + m^2) sum_a Tg^a_H Tg^a_H + 6 lambda m^2 ]_(i j) $ <eq:g0_SS-pole>
  $ delta S^((1)) supset - lf 1/ep integral dif^d x thin [ & Yt thin (bar(partial)_mu H)^dagger bar(partial)^mu H + 2/3 Yt thin (hat(partial)_mu H)^dagger hat(partial)^mu H \
    & - 2 (partial_mu H)^dagger (sum_a Tg^a_H Tg^a_H) partial^mu H + (6 lambda - C_H) m^2 H^dagger H ] $ <eq:g0_SS-ct>
  In the bases of @sec:bases, with $integral dif^d x thin (bar(D)_mu H)^dagger bar(D)^mu H = integral dif^d x thin (D_mu H)^dagger D^mu H + tilde(E)_H$ (@eq:evanescent) for the fermion part (@sec:g0_SSV_SSVV: the fermion loops couple only to $bar(A)$, and the $H^dagger H A$ fermion pole has no hatted structure):
  $ delta S^((1)) supset - lf 1/ep [ (Yt - 2 C_H) Q_(S S) + Yt thin tilde(E)_H + 2/3 Yt thin E_H + (6 lambda - C_H) m^2 Q_m ] , $ <eq:g0_SS-ct-basis>
  restricted here to the terms without gauge fields. The bosonic loops multiply the #box[$d$-dimensional] $Q_(S S)$, the fermion loops the barred kinetic term $Q_(S S) + tilde(E)_H$; @eq:g0_SS-ct-basis is the Higgs two-point part of the assembled counterterm of @sec:st-solution. With the constants $delta Z_H$ and $delta m^2$ of @sec:st-solution ($H_0 = Z_H^(1\/2) H$, $m_0^2 = m^2 + delta m^2$, so that $-m_0^2 H_0^dagger H_0$ gives $-(delta m^2 + m^2 delta Z_H) Q_m$), equivalently $delta S^((1)) supset delta Z_H Q_(S S) - (delta m^2 + m^2 delta Z_H) Q_m - lf Yt\/ep thin (tilde(E)_H + 2/3 E_H)$ with
  $ delta Z_H = lf 1/ep (2 C_H - Yt), quad quad delta m^2 + m^2 delta Z_H = lf 1/ep (6 lambda - C_H) m^2 . $ <eq:g0_SS-dZ>
]

#figure(
  table(
    columns: 4,
    align: (left, center, center, center),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([contribution to $N_0 \/ delta_(i j)$], [$p^2$ (d-dim.)], [$bar(p)^2$ and $hat(p)^2$], [$m^2$]),
    [(a) $Q$--$U$, $Q$--$D$, $L$--$E$], [$0$], [$(3 y_u^2 + 3 y_d^2 + y_e^2)(bar(p)^2 + 2/3 hat(p)^2)$], [$0$],
    [(b) $W$], [$-3/2 gw^2$], [$0$], [$-3/4 gw^2$],
    [(b) $B$], [$-1/2 gy^2$], [$0$], [$-1/4 gy^2$],
    [(c) $W$, $B$], [$0$], [$0$], [$0$],
    [(d)], [$0$], [$0$], [$6 lambda$],
  ),
  caption: [Standard-Model values of the contributions to $N_0$. The gluon contributes nothing.],
) <tab:g0_SS-sm>

With these values,
$ delta Z_H = - lf 1/ep (3 y_u^2 + 3 y_d^2 + y_e^2 - 3/2 gw^2 - 1/2 gy^2), \
  delta m^2 = lf m^2/ep (6 lambda - 9/4 gw^2 - 3/4 gy^2 + 3 y_u^2 + 3 y_d^2 + y_e^2), $ <eq:g0_SS-sm>
and the evanescent counterterms are $delta S^((1)) supset - lf 1/ep [(3 y_u^2 + 3 y_d^2 + y_e^2) tilde(E)_H + (2 y_u^2 + 2 y_d^2 + 2/3 y_e^2) E_H]$.

#check[*Independence of $MUV$.* The $MUV^2$ poles cancel separately in each topology: $4 - 4$ in (a), $2 - 2$ in (b), $1 - 1$ in (c) and (d).]

#check[*Known limits.* $Pphys N_0 = (Yt - 2 C_H) bar(p)^2 + (6 lambda - C_H) m^2$ is the result with an anticommuting $g5$: the single trace in @eq:g0_SS-int-a has no $g5$ left after @eq:sandwich. The renormalisation of the mass parameter, $m_0^2 = m^2 + delta m^2 = m^2 [1 + lf (6 lambda + Yt - 3 C_H)\/ep]$ (@eq:g0_SS-dZ), gives $16 pi^2 thin dif m^2 \/ dif ln mu = m^2 (12 lambda + 2 Yt - 9/2 gw^2 - 3/2 gy^2)$, the standard one-loop result; it is independent of the gauge parameter, while $delta Z_H$ contains the Feynman-gauge value $(3 - xi) C_H |_(xi = 1) = 2 C_H$ of the gauge contribution to the Higgs field anomalous dimension.]

#check[*Symmetry.* $N_0$ is proportional to $delta_(i j)$ (rigid $"SU"(2)_L times "U"(1)_Y$ invariance) and real (hermiticity of $delta S^((1))$). The gluon does not couple, and the $"SU"(2)_L$ and $"U"(1)_Y$ parts enter only through $C_H$.]

#check[*Relations to other classes.* The bosonic part $-2 C_H Q_(S S)$ and the fermionic part $Yt (bar(partial) H)^dagger bar(partial) H$ are the two-point parts of covariant operators whose $H^dagger H A$ and $H^dagger H A A$ components must appear in @sec:g0_SSV_SSVV: #box[$d$-dimensional] $(D H)^dagger D H$ for the gauge loops, $(bar(D) H)^dagger bar(D) H$ for the fermion loops, which couple only to $bar(A)$. The two-point function alone fixes only the coefficient of $hat(p)^2$, i.e. the difference of the coefficients of $E_H$ and $tilde(E)_H$ (@eq:evanescent), whose two-point vertex functions are $hat(p)^2$ and $-hat(p)^2$; the split of @eq:g0_SS-ct-basis is the one required by the $H^dagger H A$ function. Since $tilde(E)_H$ and the covariant operators are $bop_d$-closed and $bop_d E_H != 0$ (@eq:b-evanescent), the evanescent pole fixes the pole of $(Dbrk_1 dot Gamma)^((1))_(gh H^dagger H)$ in @sec:g1_cSS_cSSV: $[(Dbrk_1 dot Gamma)^((1))]_"div" = - bop_d delta S^((1)) supset 2/3 lf Yt/ep thin bop_d E_H$, i.e. $2/3 lf Yt/ep thin i integral dif^d x [(hat(partial)_mu H)^dagger hat(partial)^mu (wm H) - hat(partial)_mu (H^dagger wm) hat(partial)^mu H]$. The invariant combination of $delta Z_H$ and $delta m^2$ is the field direction $cal(O)_H$ and the parameter direction $Q_m$ of @eq:closed; their coefficients are not fixed by the Slavnov--Taylor identity. The assembled counterterm, with the gauge-field parts of these operators from @sec:g0_SSV_SSVV, is given in @sec:st-solution.]

#general[At one loop the forest is empty, and $m^2$ enters only linearly through $Top(2)$. Beyond one loop the evanescent counterterms $-lf Yt\/ep thin (tilde(E)_H + 2/3 E_H)$ are kept in $delta S^((1))$ with its hatted structure: inserted into a Higgs line of a one-loop graph, $hat(partial)^2 -> -hat(k)^2$ gives $eta\/ep^2$, i.e. a two-loop pole and an evanescent-induced finite part (@sec:m-beyond). The bosonic subgraphs of @fig:g0_SS-diagrams (b)--(d) keep #box[$d$-dimensional] numerators at every order and produce no $eta$.]
