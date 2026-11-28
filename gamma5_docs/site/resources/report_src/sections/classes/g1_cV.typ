#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Breaking of $omega A$ <sec:g1_cV>

External legs $gh^a (q)$ and $A^b_nu (-q)$, both momenta incoming, with $a, b$ any two of the twelve gauge indices. The sum of the diagrams below is $i thin (Dbrk_1 dot Gamma)^((1))_(gh^a A^b_nu)(q)$. The external momentum and the index $nu$ are $d$-dimensional; $nu$ is barred because fermions couple only to $bar(A)$. Degree $dg(G) = 3$. All gauge factors and all pairs (diagonal and mixed) are covered, for the five multiplets $Q, U, D, L, E$.

*Diagrams.* The only one-loop graph with one $Dbrk_1$ vertex and these legs is the fermion bubble of @fig:g1_cV-diagrams: the ghost leg sits on the vertex @eq:V-delta, the gauge leg on @eq:V-ffv. There is one graph per multiplet. The two terms of @eq:V-delta make the sterile component propagate on the line that enters the $Dbrk_1$ vertex (term $slh(p) Pcb$) or on the line that leaves it (term $-slh(p)' Pc$); they are the two rows (a1), (a2) of @tab:g1_cV-diagrams, parts of one graph.

#figure(
  oneloop(
    ((legs: (("G", [$gh^a (q)$#h(5mm)], 1),), mark: "delta"), (legs: (("V", [#h(6mm)$A^b_nu (-q)$], 0),))),
    (("F", $Psi_f (k + q)$, 1), ("F", $Psi_f (k)$, 1)), r: 0.75, leg: 0.9, caption: [(a)]),
  caption: [The topology of $(Dbrk_1 dot Gamma)^((1))_(gh A)$. The crossed circle is the vertex @eq:V-delta. The fermion line enters it with momentum $k$ and leaves it with $k + q$.],
) <fig:g1_cV-diagrams>

#figure(
  {
  set par(justify: false)
  table(
    columns: 5,
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [symmetry factor], [sign], [group factor]),
    [(a1)], [$Psi_f$; sterile part on the line $k$ (entering $Dbrk_1$)], [$1$], [$-1$ (fermion loop)], [$Tr_f (Tg^a_f Tg^b_f) = S_f^((a)) delta^(a b)$],
    [(a2)], [$Psi_f$; sterile part on the line $k + q$ (leaving $Dbrk_1$)], [$1$], [$-1$ (fermion loop)], [$Tr_f (Tg^a_f Tg^b_f) = S_f^((a)) delta^(a b)$],
  )},
  caption: [Field assignments of @fig:g1_cV-diagrams, $f in {Q, U, D, L, E}$. (a1) and (a2) are the two terms of one graph. The group factor vanishes for $a$ and $b$ in different gauge factors, because the generators of $"SU"(3)_c$ and $"SU"(2)_L$ are traceless on every multiplet.],
) <tab:g1_cV-diagrams>

*Integrands.* The loop momentum $k$ runs along the fermion arrow on the line from the gauge vertex to the $Dbrk_1$ vertex; $k + q$ runs from the $Dbrk_1$ vertex to the gauge vertex. In @eq:V-delta, $p = k$ and $p' = k + q$. The generators commute with the Dirac matrices, so the group trace $Tr_f$ over the colour and $"SU"(2)_L$ indices factorises from the Dirac trace $tr$. With the $Dbrk_1$ vertex leftmost in the trace and $Pc = P_(c_f)$,
$ G_"(a)" = - Tr_f (Tg^a_f Tg^b_f) intk tr[ (slh(k) Pcb - (slh(k) + slh(q)) Pc) (i sl(k))/(k^2) (i gbar^nu Pc) (i (sl(k) + sl(q)))/((k + q)^2) ] . $ <eq:g1_cV-int>
With $Tr_f (Tg^a_f Tg^b_f) = S_f^((a)) delta^(a b)$ (@eq:invariants), $-i^3 = i$ and the cyclicity of the trace:
$ G_"(a)" = i thin S_f^((a)) delta^(a b) intk (tr[ gbar^nu Pc (sl(k) + sl(q)) (slh(k) Pcb - (slh(k) + slh(q)) Pc) sl(k) ])/(k^2 (k + q)^2) . $ <eq:g1_cV-int2>

*$Rop$-operation.* At one loop $Rbar G = G$. _Dirac algebra._ With @eq:sandwich, $Pc (sl(k) + sl(q)) = (slb(k) + slb(q)) Pcb + (slh(k) + slh(q)) Pc$, and $Pcb slh(k) = slh(k) Pcb$, $Pcb Pc = 0$:
$ Pc (sl(k) + sl(q)) (slh(k) Pcb - (slh(k) + slh(q)) Pc) = (slb(k) + slb(q)) slh(k) Pcb - (hat(k) + hat(q))^2 Pc . $
The first term is (a1), the second (a2); they are the two forms of @eq:m-sterile. Multiplying by $sl(k)$ on the right, $Pcb sl(k) = slb(k) Pc + slh(k) Pcb$. The trace with a single hatted matrix vanishes (it factorises into a four-dimensional trace times $tr(ghat^mu) = 0$, @sec:bmhv), $slh(k) slh(k) = hat(k)^2$, and $tr(gbar^nu slb(v) P_(upright(L),upright(R))) = 2 bar(v)^nu$. The trace of @eq:g1_cV-int2 is
$ tr[dots] = 2 [ hat(k)^2 (bar(k) + bar(q))^nu - (hat(k) + hat(q))^2 thin bar(k)^nu ]
  = 2 [ hat(k)^2 bar(q)^nu - 2 (hat(k) dot hat(q)) thin bar(k)^nu - hat(q)^2 thin bar(k)^nu ] . $ <eq:g1_cV-num>
$g5$ has dropped out: the four-dimensional traces contain at most two barred matrices, so @eq:tr5 does not contribute and left- and right-handed multiplets give the same numerator. The terms cubic in $k$ cancel between (a1) and (a2), so the integrand is at most quadratically divergent although $dg(G) = 3$:
$ G_"(a)" = 2 i thin S_f^((a)) delta^(a b) intk (hat(k)^2 bar(q)^nu - 2 (hat(k) dot hat(q)) thin bar(k)^nu - hat(q)^2 thin bar(k)^nu)/(k^2 (k + q)^2) . $ <eq:g1_cV-G>

_Expansion._ With $D_MUV = 1\/(k^2 - MUV^2)$ and @eq:uv-expansion,
$ 1/k^2 & = D_MUV (1 - MUV^2 D_MUV + dots), \
  1/((k + q)^2) & = D_MUV (1 + (-2 k dot q - q^2 - MUV^2) D_MUV + 4 (k dot q)^2 D_MUV^2 + dots) . $ <eq:g1_cV-T>
$Top(3)$ keeps the terms of total degree $<= 3$ in $q$ and $MUV$. Each term of @eq:g1_cV-G carries one explicit factor $bar(q)^nu$ or $hat(q)$, so the loop integrals needed are the tensor bubbles to degree $2$ and $1$:
$ Top(2) intk (k^alpha k^beta)/(k^2 (k + q)^2) & -> intk k^alpha k^beta [D_MUV^2 - 2 MUV^2 D_MUV^3 - q^2 D_MUV^3 + 4 (k dot q)^2 D_MUV^4], \
  Top(1) intk (k^alpha)/(k^2 (k + q)^2) & -> intk k^alpha (-2 k dot q) D_MUV^3 . $ <eq:g1_cV-T2>
By @eq:tensor-red and the vacuum poles below @eq:logpole ($intk k^2 D_MUV^2 -> ilf 2 MUV^2\/ep$, $intk (k^2)^a D_MUV^(a+2) -> ilf 1\/ep$), with $d -> 4$ in the reduction coefficients,
$ Kop Top(2) intk (k^alpha k^beta)/(k^2 (k + q)^2) & = ilf 1/ep [ (1/2 MUV^2 - 1/2 MUV^2 - 1/4 q^2 + 1/6 q^2) g^(alpha beta) + 1/3 q^alpha q^beta ] \
  & = ilf 1/ep [ 1/3 q^alpha q^beta - 1/12 q^2 g^(alpha beta) ] , $ <eq:g1_cV-bub2>
$ Kop Top(1) intk (k^alpha)/(k^2 (k + q)^2) = - ilf 1/(2 ep) q^alpha . $ <eq:g1_cV-bub1>
The term $4 (k dot q)^2$ gives $4 q_gamma q_delta (g^(alpha beta) g^(gamma delta) + g^(alpha gamma) g^(beta delta) + g^(alpha delta) g^(beta gamma)) \/ 24 = (q^2 g^(alpha beta) + 2 q^alpha q^beta)\/6$. The two $MUV^2$ terms multiply the same tensor $g^(alpha beta)$ and cancel identically, also after the contraction with $gmh_(alpha beta)$ below.

_Tensor reduction with the tag._ The three terms of @eq:g1_cV-G contract @eq:g1_cV-bub2 and @eq:g1_cV-bub1 with split metrics. The tag arises only in the first term, from $gmh_(alpha beta) g^(alpha beta) = gmh^alpha_(space alpha) = -2 eta$ (@eq:eta); in the second, $gmh^beta_(space alpha) gmb^(nu alpha) = 0$:
$ hat(k)^2 bar(q)^nu &-> ilf 1/ep [ 1/3 hat(q)^2 + 1/6 eta thin q^2 ] bar(q)^nu, \
  -2 (hat(k) dot hat(q)) thin bar(k)^nu &-> ilf 1/ep [ -2/3 hat(q)^2 ] bar(q)^nu, \
  -hat(q)^2 thin bar(k)^nu &-> ilf 1/ep [ 1/2 hat(q)^2 ] bar(q)^nu . $ <eq:g1_cV-red>
In the first line, $gmh_(alpha beta) (1/3 q^alpha q^beta - 1/12 q^2 g^(alpha beta)) = 1/3 hat(q)^2 + 1/6 eta thin q^2$; $q^2 = bar(q)^2 + hat(q)^2$ is $d$-dimensional. The sum is $(1/6 hat(q)^2 + 1/6 eta thin q^2) bar(q)^nu$. With the prefactor $2 i S_f^((a)) delta^(a b)$ of @eq:g1_cV-G,
$ Top(3) G_"(a)" = ilf ( 1/ep N_0^"(a)" + eta/ep N_1^"(a)" ) + dots, quad
  N_0^"(a)" = i/3 S_f^((a)) delta^(a b) thin hat(q)^2 bar(q)^nu, quad
  N_1^"(a)" = i/3 S_f^((a)) delta^(a b) thin q^2 bar(q)^nu . $ <eq:g1_cV-N>
(a1) and (a2) are not separated further: their cubic terms cancel only in the sum. $N_1^"(a)"$ comes entirely from the term $hat(k)^2 bar(q)^nu$ of @eq:g1_cV-G. $Pphys N_1^"(a)" = i/3 S_f^((a)) delta^(a b) bar(q)^2 bar(q)^nu$; the evanescent part $Evan N_1^"(a)" = i/3 S_f^((a)) delta^(a b) hat(q)^2 bar(q)^nu$ is not used at one loop.

*Result.* Write $sum_f$ for the sum over $Q, U, D, L, E$; $sum_f Tr_f (Tg^a_f Tg^b_f) = S^((r)) delta^(a b)$ for $a, b$ in the factor $r$, with $S^((r))$ of @sec:g1 (zero for $a$, $b$ in different factors).
#result[
  $ (Dbrk_1 dot Gamma)^((1))_(gh^a A^b_nu)(q) |_"pole" = lf 1/ep N_0, quad quad
    N_0 = i/3 sum_f Tr_f (Tg^a_f Tg^b_f) thin hat(q)^2 bar(q)^nu quad "(evanescent)", $ <eq:g1_cV-pole>
  the vertex function of the local functional
  $ (Dbrk_1 dot Gamma)^((1))|_("pole", gh A) & = - lf 1/(3 ep) sum_f integral dif^d x thin Tr_f (wm thin hat(partial)^2 bar(partial)^mu bar(Am)_mu) \
    & = lf 1/(3 ep) sum_(r = G, W, B) S^((r)) thin hat(N)^((r))_A = lf 1/(3 ep) sum_(r = G, W, B) S^((r)) thin bop_d E^((r))_A |_(gh A) , $ <eq:g1_cV-pole-op>
  with $hat(N)^((r))_A$ of @eq:Nhat.
  $ X_(1; gh^a A^b_nu)(q) = lf thin Pphys N_1, quad quad
    Pphys N_1 = i/3 sum_f Tr_f (Tg^a_f Tg^b_f) thin bar(q)^2 bar(q)^nu , $ <eq:g1_cV-X>
  $ X_1 supset - lf 1/3 sum_f integral dif^4 x thin Tr_f (wm thin partial^2 partial^mu Am_mu)
    = - lf 1/3 sum_(r = G, W, B) S^((r)) thin N^((r))_3 . $ <eq:g1_cV-X-op>
]
Translation to vertex functions (@sec:dictionary): with $gh^a (q)$ first and $A^b_nu (-q)$, $partial -> -i(-q) = i q$ on the gauge leg, so $N^((r))_3 = integral gh^a partial^2 partial dot A^a$ has the vertex function $(i q)^2 (i q^nu) delta^(a b) = -i q^2 q^nu delta^(a b)$, and $integral gh^a hat(partial)^2 bar(partial) dot bar(A)^a$ has $-i hat(q)^2 bar(q)^nu delta^(a b)$. Hence $-1/3 S^((r)) N^((r))_3$ has the vertex function $i/3 S^((r)) q^2 q^nu delta^(a b)$, whose four-dimensional part is @eq:g1_cV-X. By @eq:b-evanescent and @eq:Nhat, after integration by parts on $bar(partial)^mu$, $hat(N)^((r))_A = bop_d E^((r))_A |_(gh A) = - integral gh^a hat(partial)^2 bar(partial) dot bar(A)^a$, with vertex function $i hat(q)^2 bar(q)^nu delta^(a b)$. No term $N^((r))_1$ ($integral gh^a partial dot A^a$) arises: the fermions are massless.

#figure(
  table(
    columns: 7,
    align: (left, center, center, center, center, center, center),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$r$], [$Q$], [$U$], [$D$], [$L$], [$E$], [sum: coefficient of $N^((r))_3$ \ in $16 pi^2 X_1$]),
    [$G$], [$-gs^2\/3$], [$-gs^2\/6$], [$-gs^2\/6$], [$0$], [$0$], [$-2/3 gs^2$],
    [$W$], [$-gw^2\/2$], [$0$], [$0$], [$-gw^2\/6$], [$0$], [$-2/3 gw^2$],
    [$B$], [$-gy^2\/18$], [$-4 gy^2\/9$], [$-gy^2\/9$], [$-gy^2\/6$], [$-gy^2\/3$], [$-10/9 gy^2$],
    [mixed], [$0$], [$0$], [$0$], [$0$], [$0$], [$0$],
  ),
  caption: [Standard-Model values of @eq:g1_cV-X-op per multiplet, $-S_f^((r))\/3$ with $S_f^((r))$ of @tab:invariants. The pole @eq:g1_cV-pole has per multiplet the coefficient $i S_f^((r))\/3$ of $hat(q)^2 bar(q)^nu delta^(a b)$, in total $2 i gs^2\/3$, $2 i gw^2\/3$, $10 i gy^2\/9$ for $r = G, W, B$; equivalently the coefficient of $hat(N)^((r))_A = bop_d E^((r))_A |_(gh A)$ in @eq:g1_cV-pole-op is $+S^((r))\/3 = 2/3 gs^2, 2/3 gw^2, 10/9 gy^2$ in units of $1\/(16 pi^2 ep)$.],
) <tab:g1_cV-sm>

#check[*Cancellation of the pole.* @sec:g0_VV gives $delta S^((1)) supset -1/3 S^((r)) E^((r))_A \/ (16 pi^2 ep)$, so $bop_d delta S^((1)) |_(gh A) = -1/3 S^((r)) hat(N)^((r))_A \/ (16 pi^2 ep)$, which cancels @eq:g1_cV-pole-op exactly, as required by the check of @sec:m-oneloop. No other ghost-number-0 counterterm has a $gh A$ vertex under $bop_d$ with a pole: the four-dimensional and $bop_d$-closed parts of @sec:g0_VV are transverse.]

#check[*Regularised Slavnov--Taylor identity before expansion.* By @eq:m-QAP and @eq:m-wA-1loop, $G_"(a)" = -i q_mu thin G^(mu nu)_"(a)"$ with $G_"(a)"^(mu nu)$ the fermion bubble of @sec:g0_VV (its graph (a), with $p = q$). Contracting its numerator with $bar(q)_mu$ gives $2 (bar(k) dot bar(q)) bar(k)^nu + bar(q)^2 bar(k)^nu - bar(q)^nu bar(k)^2$. At $hat(q) = 0$, $2 bar(k) dot bar(q) = (k + q)^2 - k^2 - q^2$, and the terms that cancel a propagator are scaleless; with $bar(k)^2 = k^2 - hat(k)^2$ what remains is $bar(q)^nu hat(k)^2$, so $-i q_mu G^(mu nu)_"(a)" = 2 i S^((r))_f delta^(a b) bar(q)^nu intk hat(k)^2 \/ (k^2 (k + q)^2)$, which is @eq:g1_cV-G at $hat(q) = 0$. The pole of the $d$-dimensional contraction, summed over $f$, is the non-transverse evanescent pole $i/3 S^((r)) hat(q)^2 bar(q)^nu delta^(a b)$ found in @sec:g0_VV; it equals @eq:g1_cV-pole.]

#check[*Auxiliary mass, chirality, group structure.* $MUV$ cancels in @eq:g1_cV-bub2 before and after the contraction with $gmh_(alpha beta)$. $g5$ drops out of @eq:g1_cV-num, so $X_1$ is the same for left- and right-handed multiplets and involves only $Tr_f (Tg^a Tg^b)$; no Levi-Civita tensor and no anomaly trace $D^(a b c)$ can appear in a two-point function. The mixed functions $G W$, $G B$, $W B$ vanish diagram by diagram. The Dirac trace @eq:g1_cV-num was also checked numerically with an explicit eight-dimensional representation of the six-dimensional Clifford algebra (four barred, two hatted directions) for both chiralities. In the solution of @sec:st-solution, $N^((r))_3$ is produced only by $bop Q^((r))_L = -N^((r))_3 - N^((r))_(f 2)$ (@eq:btable-gauge); its partner $N^((r))_(f 2)$ lives in the class $gh A A$ (@sec:g1_cVV). With @eq:g1_cV-X-op, $bop Delta S^((1)) = -X_1$ requires $r[Q^((r))_L] = -1/3 S^((r))$ for the coefficient of $Q^((r))_L$ in $16 pi^2 Delta S^((1))_"ST"$.]

#check[*Literature.* For a single gauge group in the same regularisation @Belusca-Maito:2020ala, translated to our conventions,
  $ X_1 & supset lf 1/3 sum_f integral Tr_f (bar(partial)^mu wm thin bar(partial)^2 bar(Am)_mu), \
    (Dbrk_1 dot Gamma)^((1))|_("pole", gh A) & = lf 1/(3 ep) sum_f integral Tr_f (bar(partial)^mu wm thin hat(partial)^2 bar(Am)_mu) . $ <eq:g1_cV-lit>
  After integration by parts these are @eq:g1_cV-X-op and @eq:g1_cV-pole-op. In the product group only the diagonal traces $S_f^((r))$ enter.]

#general[At one loop $X_1$ is the residue of $hat(k)^2$ times a bubble: the sterile line next to $Dbrk_1$ turns the vertex into the scalar insertion $hat(k)^2 \/ k^2$ (@eq:m-sterile), and $hat(k)^2 -> -2 eta k^2 \/ d$ yields the local term $eta \/ ep$. At two loops the same graph appears as a ghost-number-1 subgraph of class (a) of @sec:m-beyond, and its evanescent finite part $Evan [ilf N_1] = ilf thin i/3 sum_f Tr_f (Tg^a Tg^b) hat(q)^2 bar(q)^nu$, invisible here, enters through @eq:m-LP.]
