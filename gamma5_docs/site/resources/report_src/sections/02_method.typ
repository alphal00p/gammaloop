#import "../macros.typ": *
#import "../diagrams.typ": *

// local shorthands (this file only)
#let Psib = $overline(Psi)$
#let Rs = $overline(srcF)$
#let STq = $cal(S)_"q"$
#let Dl(l) = $Dbrk^([#l])$
#let Ke = $Kop_eta$
#let Kf = $Kop^((eta))$
#let sgn = $op("sgn")$
#let lit(..keys) = [_(literature: #keys.pos().map(k => cite(k)).join([, ]))_]
#let ours = [_(this report)_]

// tables and long figures of this section may break across pages
#show figure: set block(breakable: true)
#show table: set par(justify: false)

= Systematics of the method <sec:method>

This section states what is computed and translates every operator identity of the algebraic method into an operation on Feynman diagrams and Feynman rules. @tab:m-symbols lists the symbols of @sec:setup that are used here; symbols introduced in this section are defined where they first appear. Statements taken from the literature carry a citation. Statements marked #ours are rewritings for a diagrammatic pipeline made in this report; they follow from the cited ones and are checked in the explicit calculation.

#figure(
  text(size: 8.8pt, table(
    columns: (auto, 1fr, auto),
    align: (left, left, left),
    inset: 3.6pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { (bottom: 0.2pt + luma(220)) },
    table.header([symbol], [meaning], [defined in]),
    [$d = 4 - 2 ep$; $bar(x)$, $hat(x)$], [dimension; four-dimensional and $(d-4)$-dimensional parts of a vector, index, metric ($gmb$, $gmh$) or Dirac matrix ($gbar$, $ghat$); slashed: $slb(p)$, $slh(p)$], [@eq:split, @eq:metric],
    [$Pphys$, $Evan = 1 - Pphys$], [physical projection (hatted external quantities set to zero) and evanescent part], [@eq:Pproj],
    [$Psi_f = psi_f + chi_f$; $P_c$, $Pcb$], [Dirac field of multiplet $f$: physical component $psi_f = P_(c_f) Psi_f$, sterile partner $chi_f$; projector on the physical chirality $c_f$ and the opposite one], [@eq:dirac-pack],
    [$A^a_mu$, $gh^a$, $agh^a$, $b^a$, $H$], [gauge field, ghost, antighost, auxiliary field, Higgs doublet], [@sec:fields, @eq:Lgf],
    [$Tg^a_f$, $Tg^a_H$, $fg^(a b c)$], [generators and structure constants with the couplings included], [@eq:Tg, @eq:fg],
    [$srcA^mu_a$, $srcC_a$, $srcF_f$, $srcH$], [BRST sources coupled to $sop A^a_mu$, $sop gh^a$, $sop Psib_f$, $sop H^dagger$ ($Rs_f$ and $srcH^dagger$ couple to $sop Psi_f$, $sop H$)], [@eq:Lext],
    [$Phi^*$], [source paired with the field $Phi$: $A^* = srcA$, $gh^* = srcC$, $Psi^* = Rs$, $Psib^* = srcF$, $H^* = srcH^dagger$, $(H^dagger)^* = srcH$], [@eq:Phistar],
    [$S_0$], [regularised tree-level action in $d$ dimensions, sources included], [@eq:L0],
    [$sop$; $ST(cal(F))$], [BRST operator; Slavnov--Taylor operator, quadratic in $cal(F)$], [@eq:brst, @eq:ST],
    [$bop_d$, $bop$], [linearised Slavnov--Taylor operator of $S_0$ in $d$ dimensions, and its four-dimensional limit], [@eq:b-sources, @eq:b-action],
    [$Dbrk_1 = ST(S_0)$], [tree-level breaking: ghost number $1$, evanescent; Feynman rule @eq:V-delta], [@eq:Delta1],
    [$V_Y$], [vertex function of a local functional $Y$; its Feynman rule is $i V_Y$], [@sec:dictionary],
    [$MUV$, $Top(n)$, $dg(G)$], [auxiliary mass; Taylor operator of degree $n$; superficial degree of the graph $G$], [@sec:integrals],
    [$Kop$; $lf$], [pole part in $ep$ ($overline("MS")$); one-loop factor], [@sec:integrals],
    [$eta$; $N_0$, $N_1$], [tag of a closed evanescent trace, $gmh^mu_(space mu) = -2 eta$; coefficients of $1\/ep$ and $eta\/ep$ in a one-loop graph], [@eq:eta, @eq:N0N1],
    [$S^((r))$, $C^((r))_"ad"$; $D^(a b c)$], [$S^((r)) = S^((r))_"lh" + S^((r))_"rh"$, sum over left- and right-handed multiplets; adjoint Casimir; anomaly trace], [@eq:SLR, @eq:invariants, @eq:Dabc],
    [$Q_alpha$, $N_beta$, $E^((r))_A$], [basis of the restoring counterterm; basis of ghost-number-1 functionals; evanescent counterterm operator], [@tab:Qbasis, @tab:Nbasis, @eq:evanescent],
  )),
  caption: [Symbols used in this section and where they are defined in @sec:setup.],
) <tab:m-symbols>

== The objects <sec:m-objects>

*Renormalised action* #lit(<Thomsen:starterpack>, <Belusca-Maito:2023wah>).
$ S_"ren" = S_0 + delta S + Delta S, quad
  delta S = sum_(l >= 1) delta S^((l)), quad Delta S = sum_(l >= 1) Delta S^((l)), \
  S_("ren", l) = S_0 + sum_(L = 1)^l C^((L)), quad C^((L)) = delta S^((L)) + Delta S^((L)) . $ <eq:m-Sren>
$delta S^((l))$ (singular counterterm) is an integrated local functional of ghost number $0$ and dimension $<= 4$ whose coefficients are pure poles in $ep$; it contains evanescent and non-invariant operators. $Delta S^((l))$ (finite counterterm) is local, finite, of ghost number $0$ and dimension $<= 4$. In this report $Delta S^((l))$ is purely Slavnov--Taylor restoring (no finite invariant part, i.e. $overline("MS")$ for everything else) and four-dimensional, $Evan Delta S^((l)) = 0$: it is written with $gbar$, $bar(partial)$, $bar(A)$ only. This choice of its continuation to $d$ dimensions does not affect order $l$ and does affect order $l+1$. A vertex of $C^((L))$ counts as order $L$; a graph of order $l$ has $l - sum L$ loops, the sum running over its counterterm vertices.

*Effective action.* With incoming momenta and the overall $(2 pi)^d delta^((d))(sum p)$ stripped,
$ Gamma_(Phi_1 dots Phi_n)(p_1, dots, p_n) = (delta^n Gamma)/(delta Phi_1 (p_1) dots delta Phi_n (p_n)) |_("fields" = "sources" = 0), quad
  i thin Gamma_(Phi_1 dots Phi_n) = sum_(G in cal(G)(Phi_1 dots Phi_n)) G . $ <eq:m-Gamma>
$cal(G)(Phi_1 dots Phi_n)$ is the set of amputated 1PI graphs with these external legs, built with the Feynman rules of $S_"ren"$, and $G$ is the value of a graph including symmetry factor and fermion-loop sign. Fields and sources are both legs; a source leg never propagates. The single-vertex graphs are included, so the tree part of a vertex function is its Feynman rule divided by $i$. For Grassmann-odd legs, $Gamma_(Phi_1 dots Phi_n)$ is the component $V$ of @sec:dictionary: the coefficient of $theta_1 dots theta_n$ in $Gamma[sum_k theta_k e_k e^(-i p_k dot x)]$ with the polarisations $e_k$ stripped, i.e. left derivatives applied in the order written, $delta \/ delta Phi_1$ first. Exchanging two adjacent odd legs changes the sign, e.g. $Gamma_(srcA gh) = -Gamma_(gh srcA)$. A component with one $Psib$ (or $H^dagger$, $Rs$, $srcH^dagger$) leg and one $Psi$ (or $H$, $srcF$, $srcH$) leg is a matrix in spinor and multiplet space, sandwiched between the two polarisations. $Gamma = sum_l Gamma^((l))$ is the expansion in the order defined above, and
$Gamma_l = Gamma[S_("ren", l)]$ is the effective action built with the counterterms up to order $l$. $b^a$ is kept as a field, so $Gamma^((0))_(A A)$ is the transverse $-(q^2 g^(mu nu) - q^mu q^nu) delta^(a b)$ and the gauge fixing sits in $Gamma_(b A)$, which is not corrected (@eq:linear-ids).

*What "$delta S^((l)) = -Kop Rbar Gamma^((l))$" means per graph* #lit(<Thomsen:starterpack>, <Chetyrkin:1984xa>, <Herzog:2017bjx>). For a 1PI graph $G$,
$ Rbar G = G + sum_(F != emptyset) (product_(gamma in F) Z(gamma)) * G\/F, quad quad
  Z(gamma) = - Kop thin Top(dg(gamma)) thin Rbar gamma , $ <eq:m-Rbar>
where $F$ runs over the non-empty sets of pairwise disjoint, proper, superficially divergent 1PI subgraphs of $G$ ($dg(gamma) >= 0$), $G\/F$ is $G$ with every $gamma in F$ shrunk to a vertex, and $*$ inserts the local polynomial $Z(gamma)$ at that vertex. The recursion through $Rbar gamma$ generates the nested forests. The counterterm is then defined vertex by vertex:
$ i thin delta S^((l))_(Phi_1 dots Phi_n) = sum_(G in cal(G)^((l))(Phi_1 dots Phi_n)) Z(G) quad "for every" quad sum_j dim Phi_j <= 4 , $ <eq:m-dS-per-graph>
i.e. $Z(G)$ is a contribution to the Feynman rule of the counterterm vertex, and the local functional $delta S^((l))$ is recovered from the polynomial Feynman rule by $p_j -> i partial$ on leg $j$. The set $cal(G)^((l))$ contains the graphs of order $l$ built from the vertices of $S_0 + sum_(L < l) Delta S^((L))$: finite counterterms are ordinary vertices, singular counterterm insertions are not separate graphs but the forest terms of @eq:m-Rbar. $Kop$ is the pole part at $eta = ep$ (@sec:m-K); all numerator algebra is $d$-dimensional and evanescent structures are kept in $Z(gamma)$.

*Insertions* #lit(<Belusca-Maito:2023wah>). For an integrated local functional $cal(O)$ of ghost number $1$ (Grassmann odd), $cal(O) dot Gamma$ is the left derivative of $Gamma[S_"ren" + kappa cal(O)]$ with respect to a Grassmann-odd constant $kappa$ at $kappa = 0$. Diagrammatically, $i (cal(O) dot Gamma)_(Phi_1 dots Phi_n)$ is the sum of the 1PI graphs with exactly one vertex from $i cal(O)$ and all other vertices from $S_"ren"$, the single-vertex graph included; @eq:V-delta is the case $cal(O) = Dbrk_1$.

== The dictionary <sec:m-dictionary>

@tab:m-dictionary lists each operator statement with the operation it stands for on diagrams; the subsections give the formulas.

#figure(
  table(
    columns: (auto, 0.8fr, 1.6fr),
    align: (left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { (bottom: 0.2pt + luma(200)) },
    table.header([], [operator statement], [operation on diagrams / Feynman rules]),
    [(i)], [$ST(S_0) = sop S_0 = Dbrk_1$], [one new vertex $gh Psib Psi$ of ghost number $1$ (@eq:V-delta); it is the defect of the tree-level Ward identity of the gauge vertex, @eq:m-tree-WI],
    [(ii)], [$ST(Gamma) = Dbrk dot Gamma$], [left: products of two 1PI functions glued at a pair (source $Phi^*$, field $Phi$) without propagator, @eq:m-glue; right: 1PI graphs with exactly one $Dbrk$ vertex. An identity between regularised graphs at every $d$, for the sum of the graphs of a class],
    [(iii)], [$bop_d Y$, $bop Y$], [tree-level vertices of ghost number $1$ from the vertices of $Y$, with the ghost $gh^a$ of polarisation $u^a$ and incoming momentum $q$: an added gauge leg of polarisation $-i q_mu u^a$ (ghost-momentum contraction); each gauge, ghost and matter leg rotated by the ghost, its momentum shifted by $q$; with several external ghosts, each of them in turn as the ghost, and a ghost leg of $Y$ replaced by a pair of external ghosts; a source leg $Phi^*$ replaced by the amputated leg $Phi$ of a tree vertex of $S_0$ (equation of motion), @eq:m-b-rules],
    [(iv)], [$Top(dg(G))$, $Kop$, $Ke$, $Kf$, $eta$], [Taylor expansion of the integrand around the auxiliary mass $MUV$ to degree $dg(G)$; pole part at $eta = ep$ ($Kop$), at fixed $eta$ with $eta$ kept ($Kf$), or at fixed $eta$ followed by $eta -> ep$ ($Ke$); $Ke - Kop$ is the evanescent-induced local finite part, @eq:m-K],
    [(v)], [$lim_(ep -> 0) Pphys$], [external hatted momenta and indices set to zero after the loop integration and all contractions, then $ep -> 0$],
    [(vi)], [$Dl(l) = ST(S_("ren",l-1) + delta S^((l)))$], [vertices of ghost number $1$ up to order $l$: $Dbrk_1$, $bop_d delta S^((L <= l))$, $bop_d Delta S^((L < l))$ and the trees $STq(C)$, in which a counterterm vertex with a source leg $Phi^*$ is glued, without propagator, to any counterterm vertex with a leg $Phi$ (@eq:m-S-exact, @eq:m-Delta-l); each inserted once, a vertex of order $L$ into the $(l - L)$-loop graphs],
  ),
  caption: [Operator statements of the algebraic method and the corresponding operations on diagrams.],
) <tab:m-dictionary>

=== (i) The breaking vertex <sec:m-i>

$ST(S_0) = sop S_0$ (@eq:ST-S0) is obtained by applying $sop$ to every leg of every vertex of $S_0$. Only the fermion sector survives (@eq:Delta1). In momentum space let $V_2 (p) = sl(p)$ be the vertex function of the kinetic term (Feynman rule $i sl(p)$, so that $(i V_2 (p)) (i sl(p)) \/ p^2 = -1$), $V^mu = gbar^mu Pc Tg^a$ that of the gauge vertex (Feynman rule @eq:V-ffv) and
$V_(Dbrk_1) = -i Tg^a (slh(p) Pcb - slh(p)' Pc)$ that of the breaking vertex (Feynman rule @eq:V-delta; multiplet label $f$ suppressed). With $q = p' - p$ the incoming ghost momentum, #ours
$ -i q_mu V^mu [A^a_mu (q) Psib(p') Psi(p)]
  = V_(Dbrk_1)[gh^a (q) Psib(p') Psi(p)] + i [ Tg^a Pcb V_2 (p) - V_2 (p') Pc Tg^a ] , $ <eq:m-tree-WI>
In words: contracting the gauge vertex with the ghost momentum gives two inverse propagators, which cancel the adjacent propagators as in a gauge theory with an anticommuting $g5$, plus the $Dbrk_1$ vertex. @eq:m-tree-WI is the $gh Psib Psi$ vertex function of $sop S_0 = Dbrk_1$. Its left side is the gauge vertex with the gauge leg varied by $sop A^a_mu supset partial_mu gh^a$, i.e. $partial_mu -> -i q_mu$. The bracket is minus the variation of the two fermion legs of the kinetic term by $sop Psi = i gh^a Tg^a Pc Psi$ and $sop Psib = i Psib Pcb Tg^a gh^a$ (@eq:brst-matter). The two variations add up to $V_(Dbrk_1)$. The remaining variations, the term $fg^(a b c) A^b_mu gh^c$ of $sop A^a_mu$ and the variation of the fermion legs of the gauge vertex, give vertices with one more gauge leg; they cancel because $overline(psi) gbar^mu D_mu psi$ is gauge invariant in $d$ dimensions.

Inside a loop the $Dbrk_1$ vertex always sits next to a propagator of the sterile component. With @eq:phys-prop and $P_(overline(c)) sl(q) P_(overline(c)) = slh(q) P_(overline(c))$, $Pc sl(q) Pc = slh(q) Pc$, the Feynman rule @eq:V-delta and the sterile propagator combine as #ours
$ Tg^a slh(q) Pcb dot (i sl(q))/(q^2) dot Pcb = i thin Tg^a (hat(q)^2)/(q^2) thin Pcb , quad quad
  Pc dot (i sl(q))/(q^2) dot Pc dot (- slh(q) Pc Tg^a) = -i (hat(q)^2)/(q^2) thin Pc Tg^a . $ <eq:m-sterile>
The first form applies when the sterile line enters the vertex with momentum $q$, the second when it leaves it. The vertex together with its sterile propagator is a scalar insertion $plus.minus i Tg^a hat(q)^2 \/ q^2$ on a physical fermion line, carrying a ghost leg. All other fermion numerators of the graph are four-dimensional (@eq:phys-prop). In a one-loop graph through $Dbrk_1$ every Lorentz index of a bosonic vertex or propagator ends, along the index line, on a $gbar$ of the fermion line or on an external leg; loop momenta in bosonic numerators therefore enter only as $bar(k)$, contracted with an external index or momentum, or in the unsplit $d$-dimensional scalars $k^2$, $k dot p$, which are never split and give no $eta$ (@sec:m-K).

=== (ii) The regularised action principle <sec:m-ii>

*Statement* #lit(<Breitenlohner:1977hr>, <Belusca-Maito:2023wah>). For the regularised effective action $Gamma = Gamma[S_"ren"]$ in $d$ dimensions,
$ ST(Gamma) = Dbrk dot Gamma, quad Dbrk = ST(S_"ren") , $ <eq:m-QAP>
as an identity between dimensionally regularised diagrams, for every $d$, order by order.

*Left side.* $ST$ is quadratic (@eq:ST). Evaluated on the plane waves of @sec:m-objects it gives, for a component with external legs $Phi_1 dots Phi_n$ (ghost number $1$, no external $b$) #ours,
$ ST(Gamma)_(Phi_1 dots Phi_n) = sum_((Phi^*, Phi)) s_Phi sum_(J union.sq I = {1, dots, n}) sgn(J, I) thin Gamma_(Phi^* thin J)(-k) dot Gamma_(Phi thin I)(k), quad k = sum_(j in J) p_j . $ <eq:m-glue>
The pairs $(Phi^*, Phi)$ are the six source--field pairs of @eq:ST, with $Phi^*$ the source of @eq:Phistar; $J$, $I$ run over the splits of the legs into two groups, each kept in the original order. $Gamma_(Phi^* J)$ has the source leg first, $Gamma_(Phi I)$ the field leg first; the two extra legs are identified (same point, momentum $k$ flowing from $J$ into $I$, no propagator) and their indices are summed. $sgn(J, I) = plus.minus 1$ is the sign of the permutation of the Grassmann-odd legs that takes $(J, I)$ to $(1, dots, n)$. $s_Phi$ is the sign that relates the signed derivative $delta' cal(F) \/ delta Phi^*$ of @eq:ST-compact to the left derivative: $s_Phi = -1$ for $Phi = Psib_f$ (minus sign in @eq:ST) and $Phi = H^dagger$ (right derivative $delta_R \/ delta srcH$), $s_Phi = +1$ for $A$, $gh$, $Psi_f$, $H$. Five components follow, with $gh^a (q)$, $A^b_nu (p_1)$, $A^c_beta (p_2)$, $Psi$ and $H$ incoming with $p$, $Psib$ and $H^dagger$ outgoing with $p' = p + q$, and $Gamma_(Psib Psi)(p)$, $Gamma_(H^dagger H)(p)$ the two-point functions at incoming $Psi$, $H$ momentum $p$. The adjoint index $e$ is summed. Matrices are multiplied in the order written; this absorbs $sgn(J, I)$ of the matter pairs:
$ ST(Gamma)_(gh^a A^b_nu) = Gamma_(srcA^mu_e gh^a)(-q, q) thin Gamma_(A^e_mu A^b_nu)(q) , $ <eq:m-glue-wA>
$ ST(Gamma)_(gh^a A^b_nu A^c_beta) = Gamma_(srcA^mu_e gh^a) thin Gamma_(A^e_mu A^b_nu A^c_beta)
  + Gamma_(srcA^mu_e gh^a A^b_nu) thin Gamma_(A^e_mu A^c_beta)(p_2)
  + Gamma_(srcA^mu_e gh^a A^c_beta) thin Gamma_(A^e_mu A^b_nu)(p_1) , $ <eq:m-glue-wAA>
$ ST(Gamma)_(gh^a Psib Psi) = Gamma_(srcA^mu_e gh^a) thin Gamma_(A^e_mu Psib Psi)
  + Gamma_(Psib Psi)(p') thin Gamma_(Rs gh^a Psi)
  - Gamma_(gh^a Psib srcF) thin Gamma_(Psib Psi)(p) , $ <eq:m-glue-wFF>
$ ST(Gamma)_(gh^a H^dagger H) = Gamma_(srcA^mu_e gh^a) thin Gamma_(A^e_mu H^dagger H)
  + Gamma_(H^dagger H)(p') thin Gamma_(srcH^dagger gh^a H)
  + Gamma_(gh^a H^dagger srcH) thin Gamma_(H^dagger H)(p) . $ <eq:m-glue-wHH>
$ST(Gamma)_(gh^a Psib Psi A^b_nu)$ consists of the three terms of @eq:m-glue-wFF with $A^b_nu$ added, in turn, to each of their two factors as the last leg (six terms, e.g. $Gamma_(Psib Psi A^b_nu) Gamma_(Rs gh^a Psi) + Gamma_(Psib Psi)(p') Gamma_(Rs gh^a Psi A^b_nu)$), plus $Gamma_(srcA^mu_e gh^a Psib Psi) Gamma_(A^e_mu A^b_nu)$. All other terms vanish for these legs (quantum numbers, vanishing one-point functions).

*Tree level.* The tree vertex functions (Feynman rules of @sec:feynman-rules divided by $i$) are
$ & Gamma^((0))_(srcA^mu_e gh^a)(-q, q) = -i q_mu delta^(a e) = -Gamma^((0))_(gh^a srcA^mu_e), quad
  Gamma^((0))_(srcA^mu_e gh^a A^b_nu) = fg^(e b a) g_mu^(space nu), quad
  Gamma^((0))_(A^e_mu Psib Psi) = gbar^mu Pc Tg^e, \
  & Gamma^((0))_(Rs gh^a Psi) = i Tg^a Pc, quad
  Gamma^((0))_(gh^a Psib srcF) = i Pcb Tg^a, quad
  Gamma^((0))_(srcH^dagger gh^a H) = i Tg^a_H, quad
  Gamma^((0))_(gh^a H^dagger srcH) = -i Tg^a_H, \
  & Gamma^((0))_(Psib Psi)(p) = sl(p), quad Gamma^((0))_(H^dagger H)(p) = p^2 - m^2 . $
With them @eq:m-glue-wA, @eq:m-glue-wAA, @eq:m-glue-wHH and the $gh Psib Psi A$ component vanish. For @eq:m-glue-wFF, with $Pcb gbar^mu = gbar^mu Pc$ and $Pcb ghat^mu = ghat^mu Pcb$,
$ ST(S_0)_(gh^a Psib Psi) = -i slb(q) Pc Tg^a + i Tg^a (sl(p)' Pc - Pcb sl(p))
  = -i Tg^a (slh(p) Pcb - slh(p)' Pc) = V_(Dbrk_1) : $
the tree level of (ii) is (i). All signs of @eq:m-glue to @eq:m-glue-wHH were checked numerically against @eq:ST by evaluating the functionals on plane waves (@sec:bases), for $cal(F) = S_0 + P$ with $P$ a random non-invariant source-dependent local functional, and the tree values also in six dimensions.

At one loop the left side is linear in the loop functions, $ST(Gamma)^((1)) = bop_d Gamma^((1))$: one loop function glued to one tree vertex of $S_0$. For example, since $Gamma^((1))_(srcA^mu gh)(-q, q) prop q_mu$ and $Gamma^((0))_(A A)$ is transverse,
$ ST(Gamma)^((1))_(gh^a A^b_nu) = -i q_mu thin Gamma^((1) mu nu)_(A^a A^b)(q) : $ <eq:m-wA-1loop>
the ghost-momentum contraction of the one-loop gauge self-energy.

*Right side.* $Dbrk dot Gamma$ is the sum of 1PI graphs with exactly one vertex of $Dbrk$ (@sec:m-vi), e.g. at one loop the graphs of @fig:m-delta-graphs plus the tree vertex $bop_d C^((1))$.

*Why it holds for the sum of graphs* #lit(<tHooft:1972tcz>, <Breitenlohner:1977hr>, <Belusca-Maito:2023wah>). Contract the ghost momentum into the gauge leg (or apply the BRST factor to a matter or source leg) of a left-side graph and use the tree identities @eq:m-tree-WI, and their bosonic analogues which hold without defect, at the vertex where the leg is attached. Each inverse propagator produced cancels a propagator; the resulting graphs cancel pairwise along the loop after shifts of the loop momentum, or are scaleless and vanish in dimensional regularisation. What remains is the graph with the gauge vertex replaced by $Dbrk_1$. The cancellation $(i V_2 (k)) (i sl(k)) \/ k^2 = -1$ requires the propagator to be the inverse of the complete $d$-dimensional kinetic operator; this is why the sterile partners are needed (@sec:action). @eq:m-QAP at one loop is therefore a check of the Feynman rules on the sum of the graphs of a class, valid after integration (shifts of the loop momentum and vanishing scaleless integrals), not graph by graph.

#figure(
  diagrow(
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("V", $A^b_nu$, 0),))),
      (("F", $chi$, 1), ("F", none, 1)), caption: [(a) $gh A$]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("V", $A^b$, 0),)), (legs: (("V", $A^c$, 0),))),
      (("F", $chi$, 1), ("F", none, 1), ("F", none, 1)), caption: [(b) $gh A A$]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("F", $Psib$, -1),)), (legs: (("F", $Psi$, 1),))),
      (("F", none, 1), ("V", none, 0), ("F", $chi$, 1)), caption: [(c) $gh Psib Psi$]),
  ),
  caption: [One-loop graphs with one $Dbrk_1$ vertex (crossed circle). The line marked $chi$ is the sterile component next to the vertex; together they form the insertion @eq:m-sterile.],
) <fig:m-delta-graphs>

#check[@eq:m-QAP at one loop on $gh^a (q) A^b_nu (-q)$, one gauge factor $r$, with the results of @sec:g0_VV, @sec:g0_ghost_sources and @sec:g1_cV; $S^((r)) = S^((r))_"lh" + S^((r))_"rh"$. $C^((1))$ enters both sides as $bop_d C^((1))$, so the identity is @eq:m-wA-1loop $=$ graph (a) of @fig:m-delta-graphs. Poles: the bosonic ($d$-dimensional) and the four-dimensional fermion poles are transverse; the evanescent pole gives $-i q_mu (-1/3 S^((r)) hat(q)^2 gmb^(mu nu)) lf 1/ep = i/3 S^((r)) hat(q)^2 bar(q)^nu lf 1/ep$, the pole of graph (a); $bop_d$ of the counterterm $-1/3 S^((r)) lf 1/ep E^((r))_A$ (@eq:evanescent) cancels both. Finite parts: the bosonic loops are exactly transverse; the physical finite part of the fermion loop contains $-1/3 S^((r)) bar(q)^2 gmb^(mu nu) lf$, so $Pphys [-i q_mu Gamma^((1) mu nu)_(A A)]_"fin" = i/3 S^((r)) bar(q)^2 bar(q)^nu lf = X_(1; gh A)$, i.e. $X_1 supset -1/3 S^((r)) lf N^((r))_3$ (@tab:Nbasis). The fermion-loop identity holds exactly at $O(ep^0)$, $log(-q^2)$ terms included. The restoring counterterm $-1/3 S^((r)) lf Q^((r))_L$ (@tab:Qbasis) adds $1/3 S^((r)) lf bar(q)^mu bar(q)^nu$ and makes the finite $Gamma_(A A)$ transverse.]

=== (iii) The linearised operator on a counterterm <sec:m-iii>

For a local functional $Y$ of ghost number $0$, @eq:b-action gives
$ bop_d Y = sop Y + integral dif^d x sum_Phi (delta'_Y Phi) (delta S_0)/(delta Phi), quad quad sop Y = integral dif^d x sum_Phi (sop Phi) (delta Y)/(delta Phi) , $ <eq:m-bY>
where $sop Y$ includes $sop agh^a = b^a$ and $delta'_Y Phi$ are the source derivatives of @eq:b-action ($delta'_Y Psib_f = -delta Y \/ delta srcF_f$, $delta'_Y H^dagger = delta_R Y \/ delta srcH$, the others $delta Y \/ delta Phi^*$); $bop$ is the same at $d = 4$. In vertex functions (@eq:dict; Feynman rule $i V_Y$; the ghost $gh^a (u, q)$ of $bop_d Y$ listed first) $sop Y$ acts leg by leg #ours:
$ "added leg" A^a_mu (q) &: quad e_mu = -i q_mu u^a quad "(ghost-momentum contraction)", \
  A^b, gh^b &: quad e^b -> fg^(b c a) e^c u^a quad "(gauge rotation of the leg)", \
  Psi, Psib &: quad V_Y -> i thin V_Y thin Pc Tg^a u^a quad "and" quad V_Y -> -i thin u^a Tg^a Pcb thin V_Y , \
  H, H^dagger &: quad V_Y -> i thin V_Y thin Tg^a_H u^a quad "and" quad V_Y -> -i thin u^a Tg^a_H thin V_Y , $ <eq:m-b-rules>
each rotated leg with momentum $p_i + q$, summed over the legs. These are @eq:dict-i and @eq:dict-ii; the matrix forms act on the polarisation of the leg, and the minus sign on $Psib$ is $s_i = -1$ of @eq:dict-ii. With several external ghosts, each of them in turn is the ghost of these rules (times the sign of moving it to the first position), and a ghost leg of $Y$ is replaced by any two external ghosts $gh(u_i, p_i) gh(u_j, p_j)$, $i < j$: one leg of polarisation $-fg^(b c e) u_i^c u_j^e$ and momentum $p_i + p_j$, listed first, times the sign of the permutation that takes $(i, j, "others")$ to the original order. With a single further ghost this is the rotation $gh^b$ of the second line. The first and third lines applied to the kinetic and gauge vertices of $S_0$ give back @eq:m-tree-WI.

*Source legs.* The second term of @eq:m-bY acts only on source-dependent $Y$ (@eq:dict, (iii)). If $Y$ has a source leg $Phi^*(k)$, multiply the vertex function of $Y$ by that of a tree vertex of $S_0$ with a leg $Phi(-k)$, drop the two legs $Phi^*$ and $Phi$ (no propagator between them) and keep the sign of @eq:b-sources. The remaining legs of the $S_0$ vertex become legs of $bop_d Y$. Example: the source term $i integral Rs wm psi$ glued to the $d$-dimensional kinetic vertex gives $i sl(p)' Tg^a Pc$, with $p'$ the outgoing momentum of the $Psib$ leg and $sl(p)'$ unsplit. Example 3 of @sec:dictionary uses the four-dimensional limit $bop$ and gives the barred $i slb(p)' Tg^a_f P_(c_f)$.

*Grassmann-odd functionals.* For a local functional $X$ of ghost number $1$ (Grassmann odd), $bop_d X$ and $bop X$ are given by the same rules, because $bop$ acts as a left derivation: the ghost of $bop$ stands to the left of all legs of $X$. $V_(bop X)(gh(u_0, q_0); Phi_1, dots, Phi_n)$ is the sum over the external ghosts, each in turn taken as the ghost of @eq:m-b-rules, times the sign of the permutation that moves it to the first position, of @eq:m-b-rules applied to the vertex function of $X$ on the remaining legs (in their order); there is no further sign. Example: for $X = integral gh^B (partial dot B)^2$, $V_X (gh(u, q), B(e_1, p_1), B(e_2, p_2)) = -2 u (p_1 dot e_1)(p_2 dot e_2)$ and $V_(bop X)(gh(u_0, q_0), gh(u_1, q_1), B(e, p)) = 2 i (p dot e) u_0 u_1 (q_0^2 - q_1^2)$: the term with $gh(u_0, q_0)$ as the ghost of $bop$ minus the term with $gh(u_1, q_1)$. This rule was checked numerically against the action of $bop$ on plane waves (@sec:bases), on this example and on $X = integral gh^B overline(psi)_U i sl(partial) psi_U$ with external gluon and hypercharge ghosts, where the fermion legs are rotated. $bop X = 0$ is then the statement that these ghost-number-2 vertex functions vanish in every class; it is the Wess--Zumino check of @sec:st-solution.

*Example: the fermion source counterterm.* For $Y_R = z^R_f integral [Rs_f sop Psi_f - (sop Psib_f) srcF_f]$, $sop Y_R = 0$ ($sop^2 = 0$) and the second term is $z^R_f$ times the variation of $S_0$ under the fermion BRST transformation alone. Its $gh Psib Psi$ vertex function is, from @eq:m-tree-WI,
$ V_(bop_d Y_R)[gh^a (q) Psib(p') Psi(p)] = z^R_f thin V_(Dbrk_1) + i z^R_f thin slb(q) Pc Tg^a , $ <eq:m-bYR>
with further $gh Psib Psi A$ and $gh Psib Psi H$ vertices. A divergent $Y_R$ therefore produces a pole times the evanescent $Dbrk_1$ structure; this term cancels the divergent part of the $gh Psib Psi$ graphs with a $Dbrk_1$ vertex #lit(<Belusca-Maito:2020ala>, <Lazopoulos:2026two>). The complete $bop$ variations of the operator bases used for $delta S^((1))$ and $Delta S^((1))$ are given in @sec:btable.

=== (iv) Taylor operator, pole parts and the tag $eta$ <sec:m-K>

$Top(dg(G))$ is the integrand-level expansion @eq:uv-expansion to degree $dg(G)$ (equivalently: rescale every momentum and mass external to $G$, and $MUV$ in the numerators, by $t$, expand in $t$ to order $t^(dg(G))$, set $t=1$). $Top(dg(G)) G$ is a sum of one-scale vacuum integrals with mass $MUV$ on every line and has the same pole part as $G$. Every one-loop graph then has the form @eq:N0N1. Three pole operators are used on a Laurent series in $ep$ whose coefficients are polynomials in $eta$:
- $Kf$ keeps every term $eta^r ep^(-k)$ with $k > 0$, for all $r >= 0$, and keeps $eta$ as a symbol #lit(<Lazopoulos:2026two>)\;
- $Ke = [Kf thin dot thin]_(eta = ep)$ is the same followed by $eta -> ep$;
- $Kop$ sets $eta = ep$ first and then takes the pole part ($overline("MS")$).

At one loop
$ Kop thin Top(dg(G)) G = ilf N_0 / ep, quad quad
  Kf thin Top(dg(G)) G = ilf (N_0 + eta N_1) / ep, \
  Ke thin Top(dg(G)) G - Kop thin Top(dg(G)) G = ilf thin N_1 . $ <eq:m-K>
$Ke$ and $Kop$ differ by the evanescent-induced local finite part $ilf N_1$. In general a term $eta^r ep^(-k)$ ($k > 0$) is kept by $Kf$ and $Ke$ for all $r$ and by $Kop$ only for $r < k$: after $eta = ep$ the difference consists of finite local terms ($r = k$) and terms of positive order in $ep$ ($r > k$), which matter only when inserted into higher loops. $Kf$ and $eta -> ep$ do not commute: $Kf (eta \/ ep)|_(eta = ep) = 1$, while $Kop (eta \/ ep|_(eta = ep)) = 0$.

The tag is unambiguous only with a fixed rule for where $d - 4$ is written as $-2 eta$: here, the traces $gmh^mu_(space mu)$ of @eq:eta, including those produced by the reduction of hatted loop momenta @eq:tensor-red-split, $hat(k)^2 -> -2 eta k^2 \/ d$. A hatted loop-momentum component is never rewritten as $hat(k)^2 = k^2 - bar(k)^2$, which would turn the tagged term $eta N_1 \/ ep$ into an untagged finite term, discarded together with the $MUV$-dependent finite part. With this rule the one-loop ghost-number-0 graphs have $N_1 = 0$: by @eq:phys-prop their fermion numerators are four-dimensional, bosonic numerators are $d$-dimensional and never split, and hatted loop momenta arise only from $hat(k) dot hat(p)$ in the denominators, giving $hat(p)$-dependent (evanescent) structures but no trace $gmh^mu_(space mu)$. Hence $Kop = Ke = Kf$ on these graphs and $delta S^((1))$ does not depend on the choice.

*The tag at any loop order.* The rule above is the one-loop case of the following. The tensor reduction of the complete integrand is done in $d$ dimensions with the full metric ($g^mu_(space mu) = d$, never tagged). Afterwards a metric contracted with a hatted object becomes $gmh$; each closed cycle of hatted metrics gives $gmh^mu_(space mu) = -2 eta$, each closed cycle of $gmb$ gives $4$. The $ep$ dependence of reduction coefficients and vacuum integrals is never tagged #lit(<Lazopoulos:2026two>). The decomposition of a tensor integral into products of metrics and momenta is unique for generic $d$. A subgraph that is integrated first and then inserted into its cograph therefore carries the same power of $eta$ as in the complete graph, which is what a subtraction at fixed $eta$ requires (@sec:m-beyond).

=== (v) Physical projection and the limit <sec:m-v>

$lim_(ep -> 0) Pphys X$ for a sum $X$ of graphs means: carry out the loop integration and all index contractions with $d$-dimensional external kinematics, set every remaining hatted external quantity ($hat(p)_i$, $gmh$ and $ghat$ with external indices) to zero (@eq:Pproj), then take $ep -> 0$. $Pphys$ does not commute with the loop integration: inside the loop $hat(k)^mu hat(k)_mu$ produces $eta$, which survives $Pphys$. $Pphys$ is applied to the complete expression after the gluing contractions of @eq:m-glue. In the counterterms nothing is projected: $delta S^((l))$ and $Delta S^((l))$ enter higher orders in their $d$-dimensional form.

=== (vi) The insertion operator and the restoring condition at order $l$ <sec:m-vi>

Split @eq:ST as $ST(cal(F)) = STq(cal(F)) + integral b^a delta cal(F) \/ delta agh^a$, with $STq$ the sum of the source--field products. Since $ST$ is quadratic, for any local $C$ of ghost number $0$ (exact)
$ ST(S_0 + C) = Dbrk_1 + bop_d C + STq(C) . $ <eq:m-S-exact>
$STq(C) = integral sum_Phi (delta C \/ delta Phi^*)(delta C \/ delta Phi)$ (signs as in @eq:ST) is a tree in which a vertex of $C$ with a source leg $Phi^*$ is glued, without propagator, to a vertex of $C$ with a leg $Phi$; it vanishes unless $C$ depends on sources, and only the first factor needs a source. The starter pack's operator $Dbrk_l = ST(S_("ren", l-1) + delta S^((l)))$ #lit(<Thomsen:starterpack>) is written here $Dl(l)$ to keep $Dbrk_1$ for the tree breaking @eq:Delta1. Its parts up to order $l$ are
$ Dl(l) = Dbrk_1 + sum_(L = 1)^l bop_d delta S^((L)) + sum_(L = 1)^(l - 1) bop_d Delta S^((L))
  + sum_(L = 2)^l [STq(C)]^((L)) + ("terms of order" > l) , $ <eq:m-Delta-l>
where $[STq(C)]^((L))$ collects the pairs of $C^((L_1))$, $C^((L_2))$ with $L_1 + L_2 = L$, $L_1, L_2 >= 1$ (so both factors have order $< l$). Each term is a set of tree-level vertices of ghost number $1$ and dimension $<= 5$, obtained with @eq:m-b-rules and @eq:m-glue. A vertex of order $L$ is inserted once into the $(l - L)$-loop 1PI graphs with ghost-number-1 external legs; $L = l$ is the vertex itself.

*Restoring condition* #lit(<Thomsen:starterpack>, <Belusca-Maito:2023wah>). $Delta S^((l))$ is fixed by the requirement that the renormalised effective action satisfy the Slavnov--Taylor identity in four dimensions at order $l$ (eqs. (2.10)--(2.12) of #cite(<Thomsen:starterpack>)):
$ lim_(ep -> 0) Pphys [ST(Gamma_l)]^((l)) = 0 , quad quad Gamma_l = Gamma[S_("ren", l)] . $ <eq:m-restore>
In each glued product of @eq:m-glue the two factors share only the indices of the identified leg. After the loop integrations each factor depends on hatted quantities only through external hatted momenta and through $gmh$ and $ghat$ with external indices, and the contraction of the indices of one leg cannot close them into a trace $gmh^mu_(space mu)$. Hence $Pphys$ of a glued product is the product of the projected factors contracted in four dimensions, and @eq:m-restore is $ST(Gamma_"phys") = 0$ at order $l$, with $Gamma_"phys" = lim_(ep -> 0) Pphys Gamma_l$ and the four-dimensional gluing. To evaluate it, apply @eq:m-QAP to $S_("ren", l)$. With @eq:m-S-exact, $ST(S_("ren", l)) = Dl(l) + bop_d Delta S^((l)) + ("terms of order" > l)$, and a vertex of order $l$ contributes at order $l$ only as the tree vertex itself, so
$ [ST(Gamma_l)]^((l)) = (Dl(l) dot Gamma_l)^((l)) + bop_d Delta S^((l)) . $
$Delta S^((l))$ is four-dimensional, so $Pphys bop_d Delta S^((l)) = bop Delta S^((l))$. @eq:m-restore then reads #lit(<Thomsen:starterpack>)
$ X_l = lim_(ep -> 0) Pphys thin (Dl(l) dot Gamma_l)^((l)), quad quad bop thin Delta S^((l)) = - X_l . $ <eq:m-Xl>
$X_l$ is local, finite, of ghost number $1$ and dimension $<= 5$; it satisfies $bop X_l = 0$, and $X_l$ is $bop$-exact when the gauge anomaly cancels, $D^(a b c) = 0$ (@eq:Dabc) #lit(<Belusca-Maito:2023wah>). @eq:m-Xl defines $X_l$ through the complete sum of graphs; it needs no tag $eta$.

== The algorithm at loop order $l$ <sec:m-algorithm>

The steps below are the order-$l$ step of the construction; @sec:m-oneloop gives their content at $l = 1$. "New" lists what a QCD calculation with an anticommuting $g5$ does not need. "gammaloop" names the component of the gammaloop tool chain (with its diagram generator HepKIT, the Dirac and colour algebra Idenso and the vacuum-integral library Vakint) that does the step or would do it, with the missing items M$n$ and requirements R$n$ of @sec:automation, where a table gives the status of each step. In this report all steps were carried out by hand and with stand-alone codes (@sec:supporting).

+ *Model.* In: $S_("ren", l-1)$ and the ghost-number-1 vertices of $Dl(l)$ known before step 5. Out: propagators @eq:propagators; vertices of $S_0$ (@sec:feynman-rules); vertices of $Delta S^((L))$, $L < l$, with a coupling order $L$; ghost-number-1 vertices with a separate order "$Dbrk$" required to be exactly $1$. New: BRST sources as legs that never propagate, with statistics opposite to their field; the Dirac field $Psi_f$ containing the sterile component; hatted Lorentz structures in vertices ($slh(p)$ in @eq:V-delta, $hat(p)^2$ in evanescent counterterms); a vertex of ghost number $1$. _gammaloop:_ the model read by HepKIT; statistics independent of spin, non-propagating sources, typed index spaces and product gauge groups are missing (M1--M4).
+ *Generate the ghost-number-0 graphs.* In: the ghost-number-0 classes of @sec:g0 (external fields and sources with $sum dim <= 4$; no external $b$, no external antighost, by @eq:linear-ids). Out: $cal(G)^((l))$ of @eq:m-dS-per-graph. New: source legs; two Higgs fields $H$, $H^dagger$ and sources $srcH$, $srcH^dagger$ of opposite statistics. _gammaloop:_ the HepKIT diagram generator, with coupling orders; Grassmann signs of ghost and source legs missing (M1, M2, M13).
+ *$Rop$-operation.* In: each graph. Out: $Rbar G$ with the forest of @eq:m-Rbar, then $Top(dg(G)) Rbar G$. New: for ghost-number-0 graphs nothing in the forest; child counterterms $Z(gamma)$ must keep their evanescent structures and hatted external momenta, which become loop momenta of the parent. _gammaloop:_ the forest and the Taylor operator with auxiliary mass of its UV module; blocked for chiral vertices by the $g5$ gate (M7); hatted child counterterms (R2).
+ *Numerator algebra.* In: integrands. Out: Laurent series in $ep$ with the tag $eta$. New: every index and momentum split into barred and hatted parts (@eq:split); $g5$ anticommutes only with $gbar$ (@eq:g5comm); traces with $g5$ and @eq:tr5; vacuum tensor reduction with split metrics (@eq:tensor-red-split); $eta$ from $gmh^mu_(space mu)$ only (@sec:m-K). _gammaloop:_ Idenso and the HepKIT tensor reducer; BMHV algebra and split-metric reduction missing (M5, M6, bug B1).
+ *Pole part and operator matching.* In: $Kop Top(dg(G)) Rbar G$, summed per class. Out: $delta S^((l))$ as coefficients of a basis of local operators of ghost number $0$ and dimension $<= 4$, by matching momentum polynomials (linear system). Check: $MUV$-independence. New: the basis contains non-invariant and evanescent operators (e.g. $bar(A) hat(partial)^2 bar(A)$); only physical chiralities are renormalised. _gammaloop:_ Vakint vacuum integrals and the integrated counterterms of the UV module (pole part); basis construction and matching missing (M8, M9).
+ *Build $Dl(l)$.* In: $delta S^((l))$, the stored $bop_d C^((L))$ and $C^((L))$ for $L < l$. Out: vertices of @eq:m-Delta-l via @eq:m-b-rules and @eq:m-S-exact. New: entirely. _gammaloop:_ no component (M10).
+ *Generate the insertion graphs.* In: the ghost-number-1 classes of @sec:g1 (ghost number $1$, $sum dim <= 5$). Out: graphs of order $l$ with exactly one vertex of $Dl(l)$. New: the coupling-order constraint "one $Dbrk$ vertex"; external ghosts and sources together with fermions. _gammaloop:_ the HepKIT generator with coupling order "$Dbrk$" $= 1$ (M1, M2).
+ *Evaluate and extract $X_l$.* In: insertion graphs. Out: their divergent part, which must cancel within each class against the tree vertices $bop_d delta S^((l)) + [STq(C)]^((l))$ (check), and $X_l$ of @eq:m-Xl. Steps 3--5 apply with the subtraction rule of @tab:m-bookkeeping (@sec:m-beyond). In bookkeeping B, in this order: integrate at fixed $eta$, take $Kf Top(dg(G)) Rbar_eta G$ per graph, sum over the graphs of the class, apply $Pphys$, set $eta = ep$, keep the $ep^0$ term, add the finite part of the trees. Checks: the poles cancel, and the readout is independent of $MUV$; at $l = 2$ so is each coefficient $eta^r ep^(-k)$. New: the finite part is kept; at one loop it is the coefficient $N_1$ of $eta\/ep$, projected with $Pphys$. _gammaloop:_ returns poles only; $Kf$, $Ke$, $Pphys$ and a finite-part output are missing (M8), and beyond one loop a subtraction rule keyed by ghost number (R1, R3).
+ *Solve $bop Delta S^((l)) = -X_l$.* In: $X_l$ per class; an ansatz $Delta S^((l)) = sum_alpha r_alpha Q_alpha$ over four-dimensional local operators of ghost number $0$ and dimension $<= 4$. Out: $r_alpha$ from matching the vertex functions of $bop Q_alpha$ (@eq:m-b-rules at $d = 4$) to $-X_l$ in every class: an overdetermined linear system, solvable if $D^(a b c) = 0$ and $X_l$ satisfies the existence conditions of the chosen normalisation (@sec:normalisation; (N4$'$) if the $gh^B$ condition fails). The solution is unique up to $bop$-closed (invariant) combinations, which are fixed by convention. New: entirely. _gammaloop:_ no component (M9, M11).
+ *Update the model.* Add $delta S^((l))$ and $Delta S^((l))$ as vertices in their $d$-dimensional form, store $bop_d C^((l))$ and $C^((l))$ for order $l+1$, and for bookkeeping B also $hat(Dbrk)^((l))$ and $E^((l))$ of @eq:m-LPl. _gammaloop:_ counterterm vertices with their own coupling order (HepKIT, four-dimensional vertices); evanescent vertices and a counterterm store missing (M3, M12, R4).

== The one-loop limit <sec:m-oneloop>

At $l = 1$ only the empty forest contributes, $Rbar G = G$, and the graphs are built from $S_0$ alone.

*Singular counterterms.* For every one-loop graph $G$ of a ghost-number-0 class (@sec:g0),
$ Z(G) = -Kop thin Top(dg(G)) thin G = -ilf thin (N_0 (G))/ep , quad quad
  i thin delta S^((1))_(Phi_1 dots Phi_n) = sum_(G) Z(G) , $ <eq:m-Z1>
with $Kop = Ke$ on these graphs (@sec:m-K).

*Breaking.* $Dl(1) = Dbrk_1 + bop_d delta S^((1))$ up to terms of order $> 1$, so $(Dl(1) dot Gamma_1)^((1)) = (Dbrk_1 dot Gamma)^((1))_"loop" + bop_d delta S^((1))$: one-loop graphs with one $Dbrk_1$ vertex built from $S_0$ (@fig:m-delta-graphs), plus the tree vertex $bop_d delta S^((1))$. The poles cancel (check below), and
#eqbox[$ i thin X_(1; Phi_1 dots Phi_n) = sum_(G in cal(G)_(Dbrk_1)(Phi_1 dots Phi_n)) Pphys [ ilf thin N_1 (G) ]_(ep = 0), quad quad
  V_(bop Delta S^((1)))[Phi_1 dots Phi_n] = - X_(1; Phi_1 dots Phi_n) , $ <eq:m-X1>]
for every ghost-number-1 class (@sec:g1), where $cal(G)_(Dbrk_1)$ are the one-loop graphs with one $Dbrk_1$ vertex and $N_1 (G)$ is read off @eq:N0N1. $X_(1; Phi_1 dots Phi_n)$ is a component of $X_1$ in the sense of @eq:m-Gamma, i.e. the vertex function of the local functional $X_1$; the right equation is the restoring condition @eq:m-Xl in vertex functions.

*Why $X_1$ is local, finite and independent of $MUV$* #ours (the statement is #cite(<Belusca-Maito:2023wah>), from the Bonneau identities #cite(<Bonneau:1980zp>)). $Dbrk_1$ is evanescent: each term of @eq:V-delta has one hatted momentum, on the sterile leg. If the sterile leg is external, the graph is proportional to $hat(p)$ of that leg and $Pphys$ removes it. If it is internal, @eq:m-sterile turns the vertex and the sterile propagator into $hat(q)^2\/q^2$ with $q = k + p_0$ ($k$ the loop momentum, $p_0$ a sum of external momenta). Of $hat(q)^2 = hat(k)^2 + 2 hat(k) dot hat(p)_0 + hat(p)_0^2$ only $hat(k)^2$ survives $Pphys$: the other terms carry an external hatted momentum, and by the last paragraph of @sec:m-i every other loop momentum enters as $bar(k)$, contracted with external quantities, or in unsplit scalars. After integration $hat(k)^2$ gives one trace $gmh^mu_(space mu) = -2 eta$, so
$ Pphys thin G = eta thin Pphys thin B_G + "(evanescent)", quad quad lim_(ep -> 0) Pphys thin G |_(eta = ep) = Pphys thin "Res"_(ep = 0) B_G = Pphys [ilf thin N_1 (G)] , $ <eq:m-PG>
where $B_G$ is a one-loop integral with the external legs and the superficial degree $dg(G)$ of $G$, regular at $eta = 0$. The finite part of $B_G$ multiplies $eta = ep$ and drops out. The residue of a one-loop integral is a polynomial in the external momenta and masses of degree $dg(G)$ (local), it is the coefficient of $eta\/ep$ in $Kf Top(dg(G)) G$, and that coefficient is independent of $MUV$ identically in $eta$. Consequences: only graphs with $dg(G) >= 0$ and an internal sterile line next to $Dbrk_1$ contribute to $X_1$; the $MUV$-dependent finite terms of @eq:N0N1 are never needed; no reading point enters.

#check[$ [ (Dbrk_1 dot Gamma)^((1)) ]_"div" + bop_d delta S^((1)) = 0 , quad "i.e." quad
  sum_(G in cal(G)_(Dbrk_1)(Phi_1 dots Phi_n)) Kop thin Top(dg(G)) G + i thin V_(bop_d delta S^((1)))[Phi_1 dots Phi_n] = 0 $
for every ghost-number-1 class, with the complete $d$-dimensional structure (evanescent poles included) #lit(<Belusca-Maito:2023wah>). The second term is @eq:m-b-rules applied to the counterterm vertices @eq:m-Z1: the poles of the ghost-number-1 graphs are fixed by those of the ghost-number-0 graphs. A second check is @eq:m-QAP at one loop before expansion, e.g. @eq:m-wA-1loop against @fig:m-delta-graphs (a).]

== Beyond one loop <sec:m-beyond>

*General-order condition* #lit(<Belusca-Maito:2023wah>, <Thomsen:starterpack>). With @eq:m-Delta-l,
$ (Dl(l) dot Gamma_l)^((l)) = (Dbrk_1 dot Gamma_(l-1))^((l)) + sum_(L = 1)^(l-1) ([Dl(l)]^((L)) dot Gamma_(l-1))^((l - L)) + [Dl(l)]^((l)) , $ <eq:m-general>
with $[Dl(l)]^((L)) = bop_d delta S^((L)) + bop_d Delta S^((L)) + [STq(C)]^((L))$ for $L < l$ and $[Dl(l)]^((l)) = bop_d delta S^((l)) + [STq(C)]^((l))$. This is the condition "LIM" of #cite(<Belusca-Maito:2023wah>) with $bop Delta S^((l))$ separated. Its divergent part vanishes: $Gamma_l$ is finite at order $l$ by the definition of $delta S^((l))$, so the left side $[ST(Gamma_l)]^((l))$ of @eq:m-QAP is finite, and so is the right side. Equivalently, $bop_d delta S^((l))$ is fixed by the ghost-number-0 poles. A remaining pole signals an error in the counterterms or in the insertion graphs. The finite part gives $X_l$ by @eq:m-Xl.

*Classes at two loops* #lit(<Belusca-Maito:2021lnk>, <Lazopoulos:2026two>). At $l = 2$, @eq:m-general consists of
(a) two-loop graphs with one $Dbrk_1$ vertex, with the ghost-number-0 subdivergences subtracted by $Rbar$ (this includes the insertions of $delta S^((1))$, evanescent parts included);
(b) one-loop graphs with one $Dbrk_1$ vertex and one vertex of $Delta S^((1))$;
(c) one-loop graphs with one vertex of $bop_d delta S^((1)) + bop_d Delta S^((1))$;
(d) the tree vertices $bop_d delta S^((2)) + [STq(C)]^((2))$, $C = delta S^((1)) + Delta S^((1))$.
By ghost-number conservation a subgraph of a graph of (a) has ghost number $1$ exactly when it contains the $Dbrk_1$ vertex. These one-loop subgraphs are not subtracted by the forest of @eq:m-Rbar; their divergences cancel against (c), because $i V_(bop_d delta S^((1))) = -sum_G Kop Top(dg(G)) G$ (check of @sec:m-oneloop) #ours.

*Two bookkeepings.* The sum (a)--(d) can be organised in two ways, summarised in @tab:m-bookkeeping together with the one-loop rules. They cannot be combined term by term: each comes with its own subtraction of the ghost-number-1 subgraphs, and individual contributions are comparable only together with that subtraction #lit(<Lazopoulos:2026two>).

#figure(
  text(size: 8.6pt, {
    set table.cell(breakable: false)
    table(
    columns: (auto, 1fr, 1fr, 1.1fr),
    align: (left, left, left, left),
    inset: 3.6pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { (bottom: 0.2pt + luma(200)) },
    table.header([case], [subtracted, with which operator], [inserted], [read off]),
    [one loop, ghost number 0], [overall: $Z(G) = -Kop Top(dg(G)) G$; $Kop = Ke = Kf$ since $N_1 = 0$], [vertices of $S_0$], [$delta S^((1))$ from the poles, @eq:m-Z1 (this report)],
    [one loop, ghost number 1], [nothing; the poles $Kop Top(dg(G)) G$ cancel $i V_(bop_d delta S^((1)))$], [$Dbrk_1$; tree $bop_d delta S^((1))$], [$i X_1 = Pphys [Ke Top(dg(G)) G - Kop Top(dg(G)) G] = Pphys [ilf N_1]$ summed over graphs, @eq:m-X1 (proven, @eq:m-PG)],
    [order $l$, bookkeeping A], [ghost-number-0 subgraphs $gamma$ with $Kop Top(dg(gamma)) Rbar gamma$ (i.e. $delta S^((L))$); ghost-number-1 subgraphs not subtracted], [$Dbrk_1$; vertices of $Delta S^((L))$, $bop_d C^((L))$ ($L < l$); trees $bop_d delta S^((l))$, $[STq(C)]^((L <= l))$], [@eq:m-Xl: the complete finite part of the sum at $eta = ep$. At fixed $eta$ the tagged pole coefficients are non-local, so UV-expanded integrals do not suffice],
    [order $l$, bookkeeping B], [every proper subgraph $gamma$ with $dg(gamma) >= 0$ at fixed $eta$: ghost number $0$ with $delta S^((L))_eta$, ghost number $1$ with $Kf Top(dg(gamma)) Rbar_eta gamma$; $eta$ kept until the last integration], [$Dbrk_1$; $hat(Dbrk)^((L))$, vertices of $Delta S^((L))$ and $E^((L))$ ($L < l$, @eq:m-LPl; at $l = 2$: $hat(Dbrk)^((1))$ of @eq:m-LP, $E^((1)) = 0$); trees $bop_d delta S^((l))$, $[STq(C)]^((l))$], [@eq:m-X2B: the terms $eta^r ep^(-r)$ of $Kf Top(dg(G)) Rbar_eta G$, after $Pphys$ and $eta = ep$, plus the finite part of the trees (this report; prescription of #cite(<Lazopoulos:2026two>) at $l = 2$)],
  )}),
  caption: [Subtraction, insertions and extraction of the finite breaking at one loop and in the two bookkeepings beyond one loop. $Kf$: pole part at fixed $eta$, $eta$ kept; $Ke$: the same followed by $eta -> ep$; $Kop$: pole part at $eta = ep$ (@sec:m-K).],
) <tab:m-bookkeeping>

*Bookkeeping A* (classes (a)--(d) as listed, the organisation of #cite(<Belusca-Maito:2021lnk>)). $X_2 = lim Pphys$ of the complete sum is exact by @eq:m-Xl and needs no tag. The ghost-number-1 subgraphs $gamma$ of (a) are compensated by $bop_d delta S^((1))$ in (c), i.e. subtracted with $Kop$ at $eta = ep$. At fixed $eta$ the pole $ilf eta N_1 (gamma) \/ ep$ of $gamma$ is left over and multiplies the complete cograph, which is non-local. The $eta$-tagged pole coefficients of the sum are therefore non-local functions of the external momenta. Their counterpart is the untagged finite part of the one-loop graphs of (c) with the four-dimensional vertex $-X_1$. Read off tagged terms of UV-expanded integrals, bookkeeping A gives a wrong, $MUV$-dependent $X_2$ (toy check below). It requires the complete finite parts of all classes, with physical masses and momenta.

*Bookkeeping B* #lit(<Lazopoulos:2026two>). Let $Lambda_1$ be the local functional whose Feynman rules are the local parts (poles and $N_1$) of the one-loop insertion graphs, $i V_(Lambda_1) = sum_(G in cal(G)_(Dbrk_1)) Ke Top(dg(G)) G$. The one-loop results (@eq:m-X1 and the check of @sec:m-oneloop) give, as an identity between local functionals #ours,
$ bop_d delta S^((1)) + bop_d Delta S^((1)) = - Lambda_1 + hat(Dbrk)^((1)), quad quad
  V_(hat(Dbrk)^((1))) = Evan [lf sum_(G in cal(G)_(Dbrk_1)) N_1 (G)] + Evan thin V_(bop_d Delta S^((1))) . $ <eq:m-LP>
Inserting $-Lambda_1$ into one-loop graphs is the same as subtracting the ghost-number-1 subgraphs $gamma$ of (a) with $Ke Top(dg(gamma)) gamma = [Kf Top(dg(gamma)) gamma]_(eta = ep)$. The prescription of #cite(<Lazopoulos:2026two>) subtracts them instead with $Kf Top(dg(gamma)) gamma$, $eta$ kept until after the last integration. The two subtractions differ by $(eta\/ep - 1) ilf N_1 (gamma)$ times the cograph. That difference vanishes at $eta = ep$ in the complete value but moves terms between the $eta$-tagged and the untagged parts; only the $eta$-kept form removes the subdivergence at fixed $eta$, which reading $X_2$ off $eta$-tagged terms requires (step (3) below). Class (c) is thereby replaced by this subtraction plus the insertion of the purely evanescent vertex $hat(Dbrk)^((1))$; the physical finite breaking of the subgraph and $bop Delta S^((1))$ cancel and are not inserted. Write $Rbar_eta G$ for $Rbar G$ of @eq:m-Rbar with the ghost-number-0 subgraphs subtracted by $Z(gamma)$ and the ghost-number-1 subgraphs by $-Kf Top(dg(gamma)) gamma$. Then, with (c$'$) the one-loop graphs with one vertex of $hat(Dbrk)^((1))$,
$ i thin X_2 = lim_(ep -> 0) Pphys { sum_(G in "(a), (b), (c')") [ Kf thin Top(dg(G)) thin Rbar_eta G ]_(eta = ep)
  + i thin V_(bop_d delta S^((2)) + [STq(C)]^((2))) } , $ <eq:m-X2B>
where the poles cancel (check) and the finite terms are the $eta^r ep^(-r)$ terms of the loop graphs and the products of finite coefficients in $[STq(C)]^((2))$. For the one-loop restoring counterterms of #cite(<Lazopoulos:2026two>), $Evan bop_d Delta S^((1))$ comes from the source term $Rs sop Psi$: by @eq:m-bYR it is $z^R_f$ times $Dbrk_1$. @eq:m-LP fixes $hat(Dbrk)^((1))$ uniquely; the "local evanescent part" of the one-loop breaking in #cite(<Lazopoulos:2026two>) (sec. 6) is this functional.

*Why @eq:m-X2B is exact* #ours. Let $S_B (eta, ep)$ be the sum of the loop graphs (a), (b), (c$'$), with $eta$ kept; the trees of (d) are added separately.
+ _Values._ $S_B (ep, ep)$ plus the trees equals the sum (a)--(d) of bookkeeping A. At fixed $eta$ the two differ by the cograph with the vertex $ilf [(eta\/ep) N_1 - N_1]$, in which the second $N_1$ enters through its four-dimensional part ($-X_1$, class (c)) and its evanescent part ($hat(Dbrk)^((1))$) separately; at $eta = ep$ this vertex vanishes.
+ _No untagged physical term._ Every term of a loop graph of $S_B$, subtraction terms included, contains one ghost-number-1 vertex: $Dbrk_1$, $hat(Dbrk)^((1))$ or the subtraction vertex $ilf (N_0 + eta N_1)\/ep$, whose $N_0$ is evanescent. Each carries a hatted object or an explicit $eta$. After the integrations a hatted index line ends on an external hatted quantity, which $Pphys$ removes, or on a barred object, which gives zero, or closes into a cycle, which gives $-2 eta$ (@sec:m-K). Hence $Pphys S_B = sum_(r >= 1) eta^r s_r (ep)$.
+ _Locality at fixed $eta$._ In $Rbar_eta G$ every subdivergence is subtracted at fixed $eta$ (for the ghost-number-0 subgraphs $Kop = Kf$ at this order, since $N_1 = 0$), and all other vertices have coefficients independent of $ep$. This is the ordinary $overline("MS")$ subtraction of the $d$-dimensional tensor integral with the hatted metrics and momenta as fixed external tensors. Its overall poles are local, and by the uniqueness of the tensor decomposition (@sec:m-K) so is the pole coefficient of each power of $eta$. $Kf Top(dg(G)) Rbar_eta G$ gives them with UV-expanded integrals, independent of $MUV$.
By (1) and (2), $X_2 = lim_(ep -> 0) Pphys S_B (ep, ep) + "trees"$, and in $eta^r s_r$ at $eta = ep$ only the coefficient of $ep^(-r)$ in $s_r$ survives, since $r >= 1$. By (3) these coefficients are those of @eq:m-X2B. The argument is the Bonneau statement #lit(<Bonneau:1980zp>, <Belusca-Maito:2023wah>) in this setting: the finite breaking is a sum of pole coefficients, hence local. At $l = 1$ it reduces to @eq:m-PG.

#check[Toy graph with the structure of @fig:m-delta-graphs (a) inside a second loop (Euclidean, measure $dif^d k \/ pi^(d\/2)$). Subgraph $gamma(k) = integral_l (hat(l) + hat(k))^2 \/ (l^2 (l + k)^2)$, with $Kf Top(2) gamma = (hat(k)^2\/3 + eta k^2\/6)\/ep$, i.e. $N_0 prop hat(k)^2$, $N_1 prop k^2$. Cograph $integral_k (k dot p) gamma(k) \/ ((k^2)^2 (k + p)^2)$; (c) and (c$'$) as in @tab:m-bookkeeping with $Evan bop_d Delta S = 0$. Exact results in $ep$ and $eta$: the sums A and B agree at $eta = ep$; $Pphys S_B$ has no $eta^0$ term; its tagged coefficients, $1\/72$ at $eta\/ep^2$ and $-5\/432$ at $eta\/ep$ (times $bar(p)^2$), contain no $log p^2$, and $-5\/432$ is the exact finite part of the sum at $eta = ep$. In A at fixed $eta$ the coefficient of $eta\/ep$ is $(log p^2 + gamma_E)\/12 - 89\/432$: non-local, and its readout misses the exact value by $(log p^2 + gamma_E)\/12 - 7\/36$. The difference A $-$ B at fixed $eta$ is the cograph with $(eta\/ep - 1) N_1$, the barred and hatted parts of $N_1$ inserted separately. A test on published two-loop results, outside the model of this report: the abelian $C(R) C_2(R)$ part of the two-loop $gh A$ breaking of #cite(<Lazopoulos:2026two>) (bookkeeping B) and the corresponding result of #cite(<Belusca-Maito:2021lnk>) for an abelian chiral theory (bookkeeping A) are both $(1\/(6 ep^2) - 17\/(72 ep)) hat(p)^2 bar(p)^mu - 11\/24 thin bar(p)^2 bar(p)^mu$, in units of $g^4\/(16 pi^2)^2$ times the group factor and up to an overall phase convention.]

*Order $l$* #ours. The construction extends without new ingredients. At order $L$, with $eta$ kept, let $cal(G)^((L))_0$ be the ghost-number-0 graphs of order $L$ built from $S_0$ and the vertices of $Delta S^((L'))$ and $E^((L'))$, and $cal(G)^((L))_1$ the ghost-number-1 graphs of order $L$ with one vertex $Dbrk_1$ or $hat(Dbrk)^((L'))$ and otherwise the same vertices ($L' < L$). In both sets $Rbar_eta$ subtracts every proper subgraph $gamma$ by $-Kf Top(dg(gamma)) Rbar_eta gamma$; for the ghost-number-0 subgraphs this is the vertex $delta S^((L'))_eta$, which replaces $Z(gamma)$ of @eq:m-Rbar (at $L = 2$ the two agree). Then
$ & i thin V_(Lambda_L) = sum_(G in cal(G)^((L))_1) Kf Top(dg(G)) Rbar_eta G, quad
  i thin V_(delta S^((L))_eta) = - sum_(G in cal(G)^((L))_0) Kf Top(dg(G)) Rbar_eta G, quad
  E^((L)) = delta S^((L)) - delta S^((L))_eta, \
  & hat(Dbrk)^((L)) = bop_d delta S^((L)) + bop_d Delta S^((L)) + [STq(C)]^((L)) + Lambda_L |_(eta -> ep "in the terms" eta^r ep^(-k), thin r <= k) . $ <eq:m-LPl>
Bookkeeping B at order $l$ subtracts every proper ghost-number-0 subgraph with the vertex $delta S^((L))_eta$ and inserts $E^((L))$ as an ordinary vertex; it subtracts every proper ghost-number-1 subgraph with $Kf Top(dg(gamma)) Rbar_eta gamma$ and inserts $hat(Dbrk)^((L))$; $X_l$ is read off as in @eq:m-X2B. The three steps hold, with a weaker form of (3). (1): $delta S_eta + E = delta S$ holds identically in $eta$, and $hat(Dbrk)^((L)) - Lambda_L = bop_d C^((L)) + [STq(C)]^((L))$ at $eta = ep$. (2): the pole check and the restoring condition at order $L$ make $Pphys$ of the $r <= k$ part of $hat(Dbrk)^((L))$ vanish identically, so $Pphys hat(Dbrk)^((L))$ consists of the $r > k$ terms of $Lambda_L$, each with a factor $eta$; $E^((L))$ and $Delta S^((L))$ are ghost-number-0 vertices of graphs that contain one evanescent vertex. (3): every forest subtraction is the pole part at fixed $eta$, with one tag rule for graphs and subgraphs, and every inserted vertex is free of poles at $eta = ep$: $hat(Dbrk)^((L))$ by the pole check, $E^((L))$ because at $eta = ep$ it is minus the terms $r >= k$ of $delta S^((L))_eta$. Such a vertex changes the $ep^0$ term at $eta = ep$ only through pole coefficients, at fixed $eta$, of the graphs it enters, which are local. Since the subtraction by $delta S^((L))_eta$ together with the insertion of $E^((L))$ equals the subtraction by $delta S^((L))$, individual coefficients $eta^r ep^(-k)$ can be non-local from $l = 3$ on; only their combination in @eq:m-X2B is local and independent of $MUV$. At $l = 2$: $E^((1)) = 0$, $Lambda_1$ has no terms with $r > k$, and @eq:m-LPl is @eq:m-LP. $delta S^((L))_eta$ and $E^((L))$ are auxiliary; the action contains $delta S^((L))$, the poles at $eta = ep$, in both bookkeepings.

*Status.* Proven in this report: the one-loop rows of @tab:m-bookkeeping (@eq:m-PG). Exact by @eq:m-Xl and @eq:m-general, without further assumption: bookkeeping A at every order. Shown here by steps (1)--(3), checked on the toy graph and consistent with the published abelian two-loop $gh A$ breaking (check above): bookkeeping B at $l = 2$, @eq:m-X2B. It is the prescription of #cite(<Lazopoulos:2026two>) (sec. 4 for $eta$ and $Kf$, sec. 6 for $hat(Dbrk)^((1))$ and the two-loop condition), where it is stated without this argument. Argued here but not tested on any graph: bookkeeping B at $l > 2$ in the form @eq:m-LPl. Not done here: a two-loop calculation in the model.

*Example: the quadratic term $[STq(C)]^((2))$.* Let $C^((1))$ contain $z^R_f integral [Rs_f sop Psi_f - (sop Psib_f) srcF_f]$ (the source term $Y_R$ of @eq:m-bYR) and $k_f integral overline(psi)_f i sl(partial) psi_f$. In this report $delta S^((1))$ contains both terms. $Delta S^((1))_"ST"$ contains the source term (monomials $Q^(f r)_R$, coefficients $r[Q^(f r)_R] = -1/4 C^((r))_"ad"$ for $r = G, W$, @eq:sts-r-R) but no fermion kinetic term: the normalisation (N1) sets the coefficient $r[Q^f_K]$ of the kinetic monomial $Q^f_K$ (@tab:Qbasis) to zero (@sec:normalisation, @sec:Sren1). Hence $z^R_f$ has a pole and a finite part, and $k_f$ has only a pole. The representatives of #cite(<Belusca-Maito:2020ala>) and #cite(<Lazopoulos:2026two>) have a finite $k_f$ as well (@eq:sts-lit). In $[STq(C)]^((2))$ on $gh^a (q) Psib(p') Psi(p)$ the source vertex is glued to the kinetic vertex through the two matter terms of @eq:m-glue-wFF, with $Gamma_(Psib Psi)(p) -> k_f slb(p) Pc$, $Gamma_(Rs gh^a Psi) -> i z^R_f Tg^a Pc$ and $Gamma_(gh^a Psib srcF) -> i z^R_f Pcb Tg^a$:
$ V_([STq(C)]^((2)))[gh^a (q) Psib(p') Psi(p)] supset k_f slb(p)' Pc dot i z^R_f Tg^a Pc - i z^R_f Pcb Tg^a dot k_f slb(p) Pc = i z^R_f k_f thin slb(q) Tg^a Pc . $
This is $z^R_f$ times the variation of $k_f integral overline(psi)_f i sl(partial) psi_f$ under the fermion BRST transformation alone, the analogue of the second term of $bop_d Y_R$. The product of the poles enters the pole check. In a representative with a finite $k_f$ the product of the finite parts enters $X_2$. In the normalisation of this report $z^R_f k_f$ has no finite part; the finite-times-finite product arises instead on the legs $gh Psib Psi A$, where the finite source vertex is glued to the vertex of the finite fermion--gauge monomial $Q^(f r)_V$ of $Delta S^((1))_"ST"$ (@tab:Qbasis, coefficient @eq:sts-r-V).

*What must be kept from lower orders.*
+ The evanescent singular counterterms of $delta S^((L))$ (e.g. $bar(A) hat(partial)^2 bar(A)\/ep$) with their hatted structure: inserted into a loop, $hat(partial)^2 -> -hat(k)^2$ gives $eta\/ep^2$, i.e. a pole and a finite part.
+ The finite restoring counterterms $Delta S^((L))$ in the $d$-dimensional form chosen at order $L$; a different evanescent continuation changes $X_(L+1)$ #lit(<Belusca-Maito:2023wah>, <Lazopoulos:2026two>).
+ For bookkeeping B, the remainders $hat(Dbrk)^((L))$ and $E^((L))$ of @eq:m-LPl ($hat(Dbrk)^((1))$ of @eq:m-LP: invisible at order $1$, needed at order $2$), with $eta$ kept in the terms with $r > k$.
+ From two loops on, singular non-invariant counterterms can be non-evanescent #lit(<Belusca-Maito:2023wah>), and $STq(C)$ is non-zero when $C$ depends on BRST sources (non-abelian case) #lit(<Lazopoulos:2026two>).

#general[The one-loop statements of @sec:m-oneloop are the case $l = 1$ of @eq:m-Xl and @eq:m-general: the forest is empty, $Dl(1) = Dbrk_1 + bop_d delta S^((1))$, and the finite part is local because the only evanescent object is the single vertex $Dbrk_1$ (@eq:m-PG). From two loops on, the finite part of the sum of classes is local only after the sum. It is obtained either from complete finite parts (bookkeeping A) or, with UV-expanded integrals only, through the $eta$-tagged terms of bookkeeping B: all subtractions at fixed $eta$, the remainders of @eq:m-LPl inserted, $eta = ep$ after $Kf$ and $Pphys$. Steps (1)--(3) show that this is exact at $l = 2$, where it is the prescription of #cite(<Lazopoulos:2026two>). For $l > 2$ they are argued, not tested.]
