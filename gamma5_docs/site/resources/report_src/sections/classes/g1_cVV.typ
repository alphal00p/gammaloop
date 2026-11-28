#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Breaking of $omega A A$ and the anomaly <sec:g1_cVV>

External legs $gh^a (q)$, $A^b_mu (p_1)$, $A^c_nu (p_2)$, all momenta incoming, $q = -p_1 - p_2$; $a, b, c$ run over all twelve gauge indices. The sum of the diagrams below is $i thin (Dbrk_1 dot Gamma)^((1))_(gh^a A^b_mu A^c_nu)$. The free indices $mu, nu$ sit on $gbar$ matrices and are four-dimensional; the external momenta are $d$-dimensional. Degree $dg(G) = 4 + 4 - 2 dot 3 = 2$ (numerator degree $4$, three propagators). All combinations of gauge factors are covered; only $G G G$, $W W W$ (parity even) and the anomaly traces $D^(a b c)$ of @eq:Dabc for $B B B$, $B W W$, $B G G$, $G G G$ (parity odd) can be non-zero.

*Diagrams.* $gh^a$ couples to fermions only through @eq:V-delta, and the gauge fields couple to fermions only through @eq:V-ffv. The one-loop graphs are the fermion triangles of @fig:g1_cVV-diagrams with the two orientations of the loop, for each multiplet. Graphs with an internal gauge, ghost or Higgs line contain no $Dbrk_1$ vertex.

#figure(
  diagrow(
    oneloop(
      ((legs: (("G", [$gh^a (q)$#h(5mm)], 1),), mark: "delta"), (legs: (("V", $A^b_mu (p_1)$, 0),)), (legs: (("V", $A^c_nu (p_2)$, 0),))),
      (("F", $k$, -1), ("F", [#h(7mm)$k - p_1$], -1), ("F", [$k + q$#h(3mm)], -1)), caption: [(a)]),
    oneloop(
      ((legs: (("G", [$gh^a (q)$#h(5mm)], 1),), mark: "delta"), (legs: (("V", $A^b_mu (p_1)$, 0),)), (legs: (("V", $A^c_nu (p_2)$, 0),))),
      (("F", [$k + q$#h(3mm)], 1), ("F", [#h(7mm)$k - p_2$], 1), ("F", $k$, 1)), caption: [(b)]),
  ),
  caption: [Topologies of $(Dbrk_1 dot Gamma)^((1))_(gh A A)$: fermion triangles $Psi_f$ with the breaking vertex (crossed circle). The line labels are the momenta along the arrows. Each line next to the crossed circle is in turn the sterile component (@eq:m-sterile).],
) <fig:g1_cVV-diagrams>

#figure(
  table(
    columns: 5,
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [symmetry factor], [sign], [group factor]),
    [(a)], [$Psi_f$, $f in {Q, U, D, L, E}$], [$1$], [$-1$ (fermion loop)], [$Tr_f (Tg^a_f Tg^b_f Tg^c_f)$],
    [(b)], [$Psi_f$, $f in {Q, U, D, L, E}$], [$1$], [$-1$ (fermion loop)], [$Tr_f (Tg^a_f Tg^c_f Tg^b_f)$],
  ),
  caption: [Field assignments of @fig:g1_cVV-diagrams. (b) is (a) with $(b, mu, p_1) <-> (c, nu, p_2)$. A multiplet contributes only if it is charged under the three gauge factors of $a, b, c$.],
) <tab:g1_cVV-diagrams>

*Integrands.* In (a) the loop momentum $k$ flows from the $A^b$ vertex into the $Dbrk_1$ vertex, $k - p_1$ from the $A^c$ vertex to the $A^b$ vertex, and $k + q$ from the $Dbrk_1$ vertex to the $A^c$ vertex. In @eq:V-delta the incoming fermion momentum is $k$ and the outgoing one $k + q$. Writing the $Dbrk_1$ vertex leftmost and following the loop against the arrow,
$ G_"(a)" = - intk Tr_f tr [ & Tg^a (slh(k) Pcb - (slh(k) + slh(q)) Pc) (i sl(k))/(k^2) (i gbar^mu Pc Tg^b) \
  & times (i (sl(k) - sl(p)_1))/((k - p_1)^2) (i gbar^nu Pc Tg^c) (i (sl(k) + sl(q)))/((k + q)^2) ] , $ <eq:g1_cVV-int-a>
with $Tg = Tg_f$, $Pc = P_(c_f)$, and $G_"(b)" = G_"(a)" |_((b, mu, p_1) <-> (c, nu, p_2))$. Write $Pc = (1 + s_f g5)\/2$, with $s_f = -1$ for $f = Q, L$ and $s_f = +1$ for $f = U, D, E$ (a sign local to this subsection), and
$ Theta^(mu nu)(u, v) = tr[gbar^mu slb(u) gbar^nu slb(v) Pcb] = 2 (bar(u)^mu bar(v)^nu + bar(u)^nu bar(v)^mu - gmb^(mu nu) bar(u) dot bar(v)) + 2 i s_f thin epsilon^(mu alpha nu beta) bar(u)_alpha bar(v)_beta $ <eq:g1_cVV-Theta>
by @eq:tr5; $epsilon^(mu alpha nu beta)$ is the four-dimensional Levi-Civita tensor. The projector algebra is that of @eq:m-sterile. In the term $slh(k) Pcb$ the line $k$ is sterile: $Pcb sl(k) Pcb = slh(k) Pcb$ gives $hat(k)^2$. In the term $-(slh(k) + slh(q)) Pc$ the line $k + q$ is sterile: $Pc (sl(k) + sl(q)) Pc = (slh(k) + slh(q)) Pc$ gives $-(hat(k) + hat(q))^2$. With @eq:sandwich for the remaining projectors and $i^5 dot (-1) = -i$,
$ G_"(a)" & = -i Tr_f (Tg^a Tg^b Tg^c) intk (hat(k)^2 thin Theta^(mu nu)(k - p_1, k + q) - (hat(k) + hat(q))^2 thin Theta^(mu nu)(k - p_1, k))/(k^2 (k - p_1)^2 (k + q)^2) \
  & = -i Tr_f (Tg^a Tg^b Tg^c) intk (hat(k)^2 thin Theta^(mu nu)(k - p_1, q) - (2 hat(k) dot hat(q) + hat(q)^2) thin Theta^(mu nu)(k - p_1, k))/(k^2 (k - p_1)^2 (k + q)^2) . $ <eq:g1_cVV-num-a>
The second line uses that $Theta$ is bilinear. The two terms are the two sterile lines. The numerator is four-dimensional apart from the scalars $hat(k)^2$, $hat(k) dot hat(q)$ and $hat(q)^2$.

*$Rop$-operation.* At one loop $Rbar G = G$. Every term of the numerator in @eq:g1_cVV-num-a carries at least one external momentum, so $Top(2)$ needs the denominators only to first order. With $D_MUV = 1\/(k^2 - MUV^2)$ and @eq:uv-expansion,
$ 1/(k^2 (k - p_1)^2 (k + q)^2) = D_MUV^3 + 2 thin k dot w thin D_MUV^4 + dots, quad w = p_1 - q = 2 p_1 + p_2 , $ <eq:g1_cVV-T>
with $k dot w = bar(k) dot bar(w) + hat(k) dot hat(w)$. No $MUV^2$ enters the terms kept. Four terms contribute. Their poles follow from @eq:logpole, with @eq:tensor-red-split and its four-momentum analogue, at $d = 4$ in the coefficients:
$ hat(k)^2 Theta^(mu nu)(-p_1, q) thin D_MUV^3 & -> -eta/2 thin Theta^(mu nu)(-p_1, q), \
  2 hat(k)^2 (bar(k) dot bar(w)) Theta^(mu nu)(k, q) thin D_MUV^4 & -> -eta/6 thin Theta^(mu nu)(w, q), \
  -hat(q)^2 Theta^(mu nu)(k, k) thin D_MUV^3 & -> hat(q)^2 gmb^(mu nu), \
  -4 (hat(k) dot hat(q))(hat(k) dot hat(w)) Theta^(mu nu)(k, k) thin D_MUV^4 & -> 2/3 thin hat(q) dot hat(w) thin gmb^(mu nu) , $ <eq:g1_cVV-terms>
where $->$ means: the pole of $intk$ of the left-hand side, in units of $ilf 1\/ep$. In the first two entries $hat(k)^2 -> -2 eta k^2 \/ d$ and $hat(k)^2 bar(k)^alpha bar(k)^beta -> -2 eta (k^2)^2 gmb^(alpha beta) \/ (d (d + 2))$, with $intk k^2 D_MUV^3 -> ilf \/ ep$ and $intk (k^2)^2 D_MUV^4 -> ilf \/ ep$. In the last two, $Theta(k, k) = 2 (2 bar(k)^mu bar(k)^nu - gmb^(mu nu) bar(k)^2)$ (its Levi-Civita part vanishes) gives $-k^2 gmb^(mu nu)$ under $k^2 D_MUV^3$, and $(hat(k) dot hat(q))(hat(k) dot hat(w)) Theta(k, k) -> -1/6 (k^2)^2 thin hat(q) dot hat(w) thin gmb^(mu nu)$. The factors $eta$ come only from $hat(k)^2$ (@eq:eta): in the first two entries the hatted trace closes inside the loop; in the last two the hatted loop components are paired with $hat(q)$ and $hat(w)$ and give no $eta$. All other terms are odd in $bar(k)$ or in $hat(k)$. Summing, with $q = -p_1 - p_2$,
$ -eta/2 thin Theta(-p_1, q) - eta/6 thin Theta(2 p_1 + p_2, q) & = -eta/6 thin Theta^(mu nu)(p_1 - p_2, p_1 + p_2), \
  hat(q)^2 + 2/3 thin hat(q) dot (2 hat(p)_1 + hat(p)_2) & = -1/3 (hat(p)_1^2 - hat(p)_2^2) . $ <eq:g1_cVV-sum>
Hence $Top(2) G_"(a)" = ilf (N_0^"(a)" \/ ep + eta N_1^"(a)" \/ ep) + dots$ with
$ N_0^"(a)" = i/3 Tr_f (Tg^a Tg^b Tg^c) thin gmb^(mu nu) (hat(p)_1^2 - hat(p)_2^2), quad
  N_1^"(a)" = i/6 Tr_f (Tg^a Tg^b Tg^c) thin Theta^(mu nu)(p_1 - p_2, p_1 + p_2) , $ <eq:g1_cVV-N-a>
$ N_1^"(a)" = Tr_f (Tg^a Tg^b Tg^c) [ (2 i)/3 (bar(p)_1^mu bar(p)_1^nu - bar(p)_2^mu bar(p)_2^nu) - i/3 gmb^(mu nu) (bar(p)_1^2 - bar(p)_2^2) + 2/3 s_f thin epsilon^(mu nu alpha beta) bar(p)_(1 alpha) bar(p)_(2 beta) ] . $ <eq:g1_cVV-N1-a>
$N_1^"(a)"$ has no hatted external momentum, so $Pphys N_1^"(a)" = N_1^"(a)"$. (b) follows by $(b, mu, p_1) <-> (c, nu, p_2)$. Under this exchange the parity-even kinematic factors of @eq:g1_cVV-N-a and @eq:g1_cVV-N1-a change sign and $epsilon^(mu nu alpha beta) bar(p)_(1 alpha) bar(p)_(2 beta)$ is invariant. Therefore the even part of (a) + (b) is proportional to $Tr_f (Tg^a [Tg^b, Tg^c]) = i fg^(a b c) S_f^((r))$, and the Levi-Civita part to $Tr_f (Tg^a {Tg^b, Tg^c})$.

*Result.* Let $sum_f$ run over $Q, U, D, L, E$, and let $S^((r)) = sum_f S_f^((r))$ (@sec:g1) for $a, b, c$ in the factor $r$.
#result[
  $ (Dbrk_1 dot Gamma)^((1))_(gh^a A^b_mu A^c_nu) |_"pole" = lf 1/ep N_0, quad
    N_0 & = i/3 sum_f Tr_f (Tg^a [Tg^b, Tg^c]) thin gmb^(mu nu) (hat(p)_1^2 - hat(p)_2^2) \
    & = -1/3 S^((r)) fg^(a b c) gmb^(mu nu) (hat(p)_1^2 - hat(p)_2^2) , $ <eq:g1_cVV-pole>
  which is evanescent. It is the vertex function of the local functional
  $ (Dbrk_1 dot Gamma)^((1))|_("pole", gh A A) & = -1/3 lf 1/ep sum_(r = G, W) S^((r)) integral dif^d x thin fg^(a b c) gh^a bar(A)^b_mu hat(partial)^2 bar(A)^(c mu) \
    & = lf 1/(3 ep) sum_(r = G, W) S^((r)) thin hat(N)^((r))_(A A) , $ <eq:g1_cVV-pole-op>
  with $hat(N)^((r))_(A A)$ of @eq:Nhat.
  The finite breaking is $X_(1; gh^a A^b_mu A^c_nu) = lf thin Pphys N_1$ with
  $ Pphys N_1 & = sum_f Tr_f (Tg^a [Tg^b, Tg^c]) [ (2 i)/3 (bar(p)_1^mu bar(p)_1^nu - bar(p)_2^mu bar(p)_2^nu) - i/3 gmb^(mu nu) (bar(p)_1^2 - bar(p)_2^2) ] - 2/3 D^(a b c) epsilon^(mu nu alpha beta) bar(p)_(1 alpha) bar(p)_(2 beta) \
    & = S^((r)) fg^(a b c) [ 1/3 gmb^(mu nu) (bar(p)_1^2 - bar(p)_2^2) - 2/3 (bar(p)_1^mu bar(p)_1^nu - bar(p)_2^mu bar(p)_2^nu) ] - 2/3 D^(a b c) epsilon^(mu nu alpha beta) bar(p)_(1 alpha) bar(p)_(2 beta) , $ <eq:g1_cVV-X>
  with $D^(a b c)$ of @eq:Dabc. The sign $s_f$ of @eq:g1_cVV-N1-a turns $sum_f s_f Tr_f (Tg^a {Tg^b, Tg^c})$ into $-D^(a b c)$. As a local functional, with the matrix fields of @eq:matrix-fields,
  $ X_1 supset lf sum_f integral dif^4 x [ & -i/3 Tr_f (wm [Am_mu, partial^2 Am^mu - 2 partial^mu partial^nu Am_nu]) \
    & + 2/3 s_f thin epsilon^(mu nu rho sigma) Tr_f (wm thin partial_mu Am_nu thin partial_rho Am_sigma) ] , $ <eq:g1_cVV-X-matrix>
  i.e. in the basis of @tab:Nbasis
  $ X_1 supset lf [ & sum_(r = G, W) S^((r)) (1/3 N^((r))_(f 1) - 2/3 N^((r))_(f 2)) \
    & - 1/3 ( D_(G G G) N_(epsilon d) + D_(B B B) N_(epsilon B) + sum_(r = G, W) D_(B r r) (N^((r))_(epsilon B r r) + 2 N^((r))_(epsilon r B)) ) ] . $ <eq:g1_cVV-X-basis>
]
Here $D_(G G G)$, $D_(B B B)$, $D_(B r r)$ are the pure numbers of @sec:st-solution, defined by $D^(A C E) = D_(G G G) gs^3 d^(A C E)$, $D^(B B B) = D_(B B B) gy^3$ and $D^(B a b) = D_(B r r) c_r^2 gy delta^(a b)$ for $a, b$ in the factor $r = G, W$ with coupling $c_r$; $D^(a b c)$ is the trace of @eq:Dabc. The other index orderings follow because $D^(a b c)$ is totally symmetric; the orderings $(a, B, b)$ and $(a, b, B)$ give the two equal terms $2 N^((r))_(epsilon r B)$. The parity-even part has no symmetric ($d$-type, $Tr_f (Tg^a {Tg^b, Tg^c})$) component, so $S_i [t]$ and $N^((r))_(r B i)$ of @tab:Nbasis receive nothing from this class. The vertex functions used are $V(N^((r))_(f 1)) = fg^(a b c) g^(mu nu) (p_1^2 - p_2^2)$, $V(N^((r))_(f 2)) = fg^(a b c) (p_1^mu p_1^nu - p_2^mu p_2^nu)$, and $2 t^(a b c) epsilon^(mu nu alpha beta) p_(1 alpha) p_(2 beta)$ for $integral t^(a b c) epsilon^(mu nu rho sigma) gh^a partial_mu A^b_nu partial_rho A^c_sigma$ with $t$ totally symmetric.

#figure(
  table(
    columns: 4,
    align: (left, left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([factors], [pole $N_0$ (units \ $gmb^(mu nu) (hat(p)_1^2 - hat(p)_2^2)$)], [$16 pi^2 X_1$, parity even], [anomaly number, $sum_"lh" - sum_"rh"$]),
    [$G G G$], [$-2/3 gs^2 fg^(A B C)$], [$2/3 gs^2 N^((G))_(f 1) - 4/3 gs^2 N^((G))_(f 2)$], [$D_(G G G) = 1/2 [2 - (1 + 1)] = 0$],
    [$W W W$], [$-2/3 gw^2 fg^(I J K)$], [$2/3 gw^2 N^((W))_(f 1) - 4/3 gw^2 N^((W))_(f 2)$], [$Tr (tau^I {tau^J, tau^K}) = 0$: \ no Levi-Civita term],
    [$B B B$], [$0$], [$0$], [$D_(B B B) = 2 [6 (1/6)^3 + 2 (-1/2)^3$ \ $quad - 3 (2/3)^3 - 3 (-1/3)^3 - (-1)^3]$ \ $quad = 2 (-2/9 + 2/9) = 0$],
    [$B W W$], [$0$], [$0$], [$D_(B W W) = 3 dot 1/6 + 1 dot (-1/2) = 0$],
    [$B G G$], [$0$], [$0$], [$D_(B G G) = 2 dot 1/6 - [2/3 + (-1/3)] = 0$],
    [other], [$0$], [$0$], [$0$ (traceless generators)],
  ),
  caption: [Standard-Model values. Pole and even part: $S^((G)) = 2 gs^2$ (per multiplet $Q$: $gs^2$, $U$, $D$: $gs^2\/2$), $S^((W)) = 2 gw^2$ ($Q$: $3 gw^2\/2$, $L$: $gw^2\/2$). Anomaly numbers, with $N^c_f$ and $N^w_f$ the colour and $"SU"(2)_L$ multiplicities of $f$ ($N^c_f = 3$ for $Q, U, D$ and $1$ for $L, E$; $N^w_f = 2$ for $Q, L$ and $1$ for $U, D, E$) and $Y_f$ the hypercharge: $Tr_f (Tg^A {Tg^C, Tg^E}) = 1/2 gs^3 N^w_f d^(A C E)$ gives $D_(G G G) = 1/2 (sum_"lh" - sum_"rh") N^w_f$ over the triplets; $Tr_f (Tg^B {Tg^B, Tg^B}) = 2 gy^3 N^c_f N^w_f Y_f^3$ gives $D_(B B B) = 2 (sum_"lh" - sum_"rh") N^c_f N^w_f Y_f^3$; $Tr_f (Tg^B {Tg^I, Tg^J}) = gy gw^2 N^c_f Y_f delta^(I J)$ gives $D_(B W W) = sum_"lh" N^c_f Y_f$ over the doublets; $Tr_f (Tg^B {Tg^A, Tg^C}) = gy gs^2 N^w_f Y_f delta^(A C)$ gives $D_(B G G) = (sum_"lh" - sum_"rh") N^w_f Y_f$ over the triplets. Every Levi-Civita coefficient of @eq:g1_cVV-X-basis vanishes; "other" is $G G W$, $G W W$, $G B B$, $W B B$, $G W B$.],
) <tab:g1_cVV-sm>

#check[*Pole cancellation* (@sec:m-oneloop). The singular counterterm of @sec:g0_VV contains $-1/3 S^((r)) E^((r))_A \/ (16 pi^2 ep)$. By @eq:b-evanescent and @eq:Nhat, $bop_d E^((r))_A |_(gh A A) = hat(N)^((r))_(A A) = integral fg^(a b c) bar(A)^(b mu) gh^c hat(partial)^2 bar(A)^a_mu$, whose vertex function on $gh^a (q) A^b_mu (p_1) A^c_nu (p_2)$ is $-fg^(a b c) gmb^(mu nu) (hat(p)_1^2 - hat(p)_2^2)$. Its contribution $+1/3 S^((r)) fg^(a b c) gmb^(mu nu) (hat(p)_1^2 - hat(p)_2^2) \/ (16 pi^2 ep)$ cancels @eq:g1_cVV-pole exactly. This is the three-point counterpart of the two-point check in @sec:g0_VV.]

#check[*Auxiliary mass, Bose symmetry, Furry.* $MUV$ does not enter (@eq:g1_cVV-T); an independent evaluation with the full $Top(2)$ expansion to second order in all denominators confirms that the $MUV^2$ terms cancel. The result is symmetric under $(b, mu, p_1) <-> (c, nu, p_2)$: $fg^(a b c)$ times an antisymmetric kinematic factor, plus $D^(a b c)$ times a symmetric one. The parity-even part of every multiplet is proportional to $Tr_f (Tg^a [Tg^b, Tg^c])$, so the symmetric even traces cancel between (a) and (b), as for a vector-like fermion loop with three vector couplings. The chiral reduction leading to @eq:g1_cVV-num-a was checked numerically with explicit $8 times 8$ Dirac matrices in $d = 6$ (four barred and two hatted directions), for both chiralities. The poles satisfy the regularised action principle @eq:m-glue-wAA orientation by orientation and for both chiralities, with the fermion loops of $Gamma^((1))_(A A A)$ and $Gamma^((1))_(A A)$ computed separately, evanescent parts included. $Pphys N_1$ was recomputed with $4 times 4$ Dirac matrices and Feynman parameters, and the group factors with the explicit generators of @eq:Tg-explicit over all $12^3$ index values.]

#check[*Anomaly and relations to other classes.* By @eq:g1_cVV-X-matrix, each Weyl multiplet contributes $plus.minus 2/3 lf = plus.minus 1\/(24 pi^2)$ times $epsilon^(mu nu rho sigma) Tr_f (wm thin partial_mu Am_nu thin partial_rho Am_sigma)$. This is the magnitude of the consistent (Bardeen) anomaly of one Weyl fermion. The coefficient of the totally symmetric trace $D^(a b c)$ cannot be removed by any $bop Q_alpha$: the only parity-odd monomial with a $gh A A$ vertex, $Q^((r))_(B epsilon)$, generates $N^((r))_(epsilon B r r) - N^((r))_(epsilon r B)$ in this class (@eq:btable-gauge), not the combination $N^((r))_(epsilon B r r) + 2 N^((r))_(epsilon r B)$ of @eq:g1_cVV-X-basis. It vanishes for the content of @tab:fields (@tab:g1_cVV-sm). The parity-even part lies in $bop("span"{Q_alpha})$. Under the normalisation conditions of @sec:normalisation, only $Q^((r))_(3 f)$ and $Q^((r))_L$ produce $N^((r))_(f 1)$ and $N^((r))_(f 2)$ (@eq:btable-gauge). In this class $bop Delta S^((1)) = -X_1$ then gives $Delta S^((1))_"ST" supset -1/3 lf S^((r)) (Q^((r))_(3 f) + Q^((r))_L)$. Its $N^((r))_3$ component is the consistency condition $X_1 supset -1/3 lf S^((r)) N^((r))_3$ of @sec:g1_cV, which holds. Its $N^((r))_(f f)$ component, from $bop Q^((r))_(3 f) supset N^((r))_(f f)$, is not fixed by this class: $Q_(G 4 c)$ and $Q_(W 4 b)$ also produce $N^((r))_(f f)$ (@eq:btable-gauge), and the balance with $X_1 supset lf (-4/9 gs^2 N^((G))_(f f) - 1/3 gw^2 N^((W))_(f f))$ of @sec:g1_cVVV_cVVVV fixes their coefficients in @sec:st-solution.]

#general[At one loop the finite breaking is the coefficient $N_1$ of $eta\/ep$ from the single closed hatted trace $hat(k)^2$ on the sterile line, @eq:m-PG. Beyond one loop, the triangle of @fig:g1_cVV-diagrams is a ghost-number-1 subgraph of the two-loop graphs of class (a) in @sec:m-beyond. Its local part, the evanescent pole @eq:g1_cVV-pole together with $ilf N_1$ at fixed $eta$, is subtracted with the pole operator at fixed $eta$ (@eq:m-LP). The anomaly condition $D^(a b c) = 0$ enters each order only after the sum over multiplets; for a single multiplet the Levi-Civita part of $X_l$ is not $bop$-exact.]
