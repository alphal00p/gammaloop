#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Gauge two-point functions $A A$ <sec:g0_VV>

External legs $A^a_mu (p)$ and $A^b_nu (-p)$, both momenta incoming, with $a, b$ any two of the twelve gauge indices. The sum of the diagrams below is $i thin Gamma^((1) mu nu)_(A^a A^b) (p)$; the tree-level value is $Gamma^((0) mu nu)_(A^a A^b) = -(p^2 g^(mu nu) - p^mu p^nu) delta^(a b)$ (@sec:m-objects). External momentum and indices are $d$-dimensional. The superficial degree of divergence is $dg(G) = 2$ for every topology. All pairs $G G$, $W W$, $B B$ and the mixed pairs $W B$, $G B$, $G W$ are covered.

*Diagrams.* The one-loop 1PI graphs are the four bubbles and the two tadpoles of @fig:g0_VV-diagrams. There is no $A A H$ vertex in the unbroken phase, so no mixed boson loop occurs, and fermions and ghosts have no four-point vertex with two gauge fields.

#figure(
  stack(dir: ttb, spacing: 5mm,
    diagrow(
      oneloop(
        ((legs: (("V", $A^a_mu$, 0),)), (legs: (("V", $A^b_nu$, 0),))),
        (("F", $Psi_f$, 1), ("F", $Psi_f$, 1)), caption: [(a) fermion loop]),
      oneloop(
        ((legs: (("V", $A^a_mu$, 0),)), (legs: (("V", $A^b_nu$, 0),))),
        (("S", $H$, 1), ("S", $H$, 1)), caption: [(b1) Higgs bubble]),
      oneloop(
        ((legs: (("V", $A^a_mu$, 0), ("V", $A^b_nu$, 0))),),
        (("S", $H$, 1),), caption: [(b2) Higgs seagull]),
    ),
    diagrow(
      oneloop(
        ((legs: (("V", $A^a_mu$, 0),)), (legs: (("V", $A^b_nu$, 0),))),
        (("V", $A^d$, 0), ("V", $A^c$, 0)), caption: [(c1) gauge bubble]),
      oneloop(
        ((legs: (("V", $A^a_mu$, 0), ("V", $A^b_nu$, 0))),),
        (("V", $A^c$, 0),), caption: [(c2) gauge tadpole]),
      oneloop(
        ((legs: (("V", $A^a_mu$, 0),)), (legs: (("V", $A^b_nu$, 0),))),
        (("G", $gh$, 1), ("G", $gh$, 1)), caption: [(d) ghost loop]),
    ),
  ),
  caption: [Topologies of $Gamma^((1))_(A A)$. $A^a_mu (p)$ is on the left vertex, $A^b_nu (-p)$ on the right one.],
) <fig:g0_VV-diagrams>

#figure(
  table(
    columns: 5,
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [symmetry factor], [sign], [group factor]),
    [(a)], [$Psi_f$, $f in {Q, U, D, L, E}$], [$1$], [$-1$ (fermion loop)], [$Tr_f (Tg^a_f Tg^b_f) = S_f^((a)) delta^(a b)$],
    [(b1)], [$H$, $H$], [$1$], [$+1$], [$Tr (Tg^a_H Tg^b_H) = S_H^((a)) delta^(a b)$],
    [(b2)], [$H$], [$1$], [$+1$], [$Tr {Tg^a_H, Tg^b_H} = 2 S_H^((a)) delta^(a b)$],
    [(c1)], [$A^c$, $A^d$], [$1\/2$], [$+1$], [$fg^(a c d) fg^(b c d) = C_"ad"^((a)) delta^(a b)$],
    [(c2)], [$A^c$], [$1\/2$], [$+1$], [$fg^(a c e) fg^(b c e) = C_"ad"^((a)) delta^(a b)$, twice],
    [(d)], [$gh^c$, $agh^e$], [$1$], [$-1$ (ghost loop)], [$fg^(e a c) fg^(c b e) = -C_"ad"^((a)) delta^(a b)$],
  ),
  caption: [Field assignments of @fig:g0_VV-diagrams. The Higgs line is a complex doublet, so (b1) and (b2) have no factor $1\/2$. Each multiplet $f$ in (a) is one Dirac line $Psi_f$; only its physical chirality couples.],
) <tab:g0_VV-diagrams>

Every group factor is proportional to $delta^(a b)$ within one gauge factor, and all of them vanish for $a$ and $b$ in different factors: the generators of $"SU"(3)_c$ and $"SU"(2)_L$ are traceless, so $Tr_f (Tg^A_f Tg^I_f) = Tr_f (Tg^A_f Tg^B_f) = Tr_f (Tg^I_f Tg^B_f) = 0$ for every $f$ (and the same for $H$), and $fg^(a c d) = 0$ unless all three indices lie in one non-abelian factor. Hence $Gamma^((1))_(W B) = Gamma^((1))_(G B) = Gamma^((1))_(G W) = 0$, diagram by diagram, before any integration.

*Integrands.* The loop momentum $k$ flows into the left vertex on the lower line and $k + p$ leaves it on the upper line; in (a), (b1), (d) $k$ is along the arrow. Writing the left vertex leftmost in the fermion trace and taking the group matrices out of the Dirac trace $tr$,
$ G_"(a)" = - Tr_f (Tg^a_f Tg^b_f) intk tr[ (i gbar^mu P_(c_f)) (i sl(k))/(k^2) (i gbar^nu P_(c_f)) (i (sl(k) + sl(p)))/((k + p)^2) ] . $ <eq:g0_VV-int-a> With @eq:V-ssv ($H$ enters the left vertex with $k$ and leaves with $k + p$),
$ G_"(b1)" & = intk Tr[ i (2 k + p)^mu Tg^a_H dot i (2 k + p)^nu Tg^b_H ] (i)/(k^2 - m^2) (i)/((k + p)^2 - m^2), \
  G_"(b2)" & = intk Tr[ i g^(mu nu) {Tg^a_H, Tg^b_H} ] (i)/(k^2 - m^2) . $ <eq:g0_VV-int-b>
In (c1) the left vertex @eq:V-vvv has the legs $A^a_mu (p)$, $A^c_rho (k)$, $A^d_sigma (-k-p)$ and the right one $A^b_nu (-p)$, $A^c_(rho') (-k)$, $A^d_(sigma') (k + p)$. Write the bracket of @eq:V-vvv as $cal(V)^(lambda_1 lambda_2 lambda_3)(q_1, q_2, q_3) = g^(lambda_1 lambda_2) (q_1 - q_2)^(lambda_3) + g^(lambda_2 lambda_3) (q_2 - q_3)^(lambda_1) + g^(lambda_3 lambda_1) (q_3 - q_1)^(lambda_2)$, odd in the momenta, and abbreviate $cal(V)^(mu rho sigma) = cal(V)^(mu rho sigma)(p, k, -k-p) = g^(mu rho) (p - k)^sigma + g^(rho sigma) (2 k + p)^mu - g^(sigma mu) (k + 2 p)^rho$. The left vertex is $fg^(a c d) cal(V)^(mu rho sigma)$ and the right one is $fg^(b c d) cal(V)^(nu rho' sigma')(-p, -k, k + p) = -fg^(b c d) cal(V)^(nu rho' sigma')$, so
$ G_"(c1)" & = 1/2 intk fg^(a c d) cal(V)^(mu rho sigma) (-fg^(b c d) cal(V)^(nu rho' sigma')) (-i g_(rho rho'))/(k^2) (-i g_(sigma sigma'))/((k + p)^2) \
  & = C_"ad"^((a)) delta^(a b) / 2 intk (cal(V)^(mu rho sigma) cal(V)^nu_(space rho sigma))/(k^2 (k + p)^2) . $ <eq:g0_VV-int-c1>
In (c2) the Feynman rule @eq:V-vvvv is taken with the legs $A^a_mu$, $A^b_nu$, $A^c_rho$, $A^c_sigma$, summed over $c$ and contracted with the propagator $-i g_(rho sigma) \/ k^2$. Its first line vanishes ($fg^(c c e) = 0$); each of the other two lines gives $fg^(a c e) fg^(b c e) (d - 1) g^(mu nu)$ after the contraction with $g_(rho sigma)$, so that the rule gives $-2 i (d - 1) C_"ad"^((a)) delta^(a b) g^(mu nu)$ and
$ G_"(c2)" = 1/2 intk [-2 i (d - 1) C_"ad"^((a)) delta^(a b) g^(mu nu)] (-i)/(k^2) = -(d - 1) C_"ad"^((a)) delta^(a b) g^(mu nu) intk 1/(k^2) . $ <eq:g0_VV-int-c2> In (d) the ghost enters the left vertex with $k$ and the antighost leaves it with $k + p$ (@eq:V-ghost: $-fg^(e a c) (k + p)^mu$), and at the right vertex the antighost leaves with $k$ ($-fg^(c b e) k^nu$):
$ G_"(d)" = - intk (-fg^(e a c) (k + p)^mu) (-fg^(c b e) k^nu) (i)/(k^2) (i)/((k + p)^2) = - C_"ad"^((a)) delta^(a b) intk ((k + p)^mu k^nu)/(k^2 (k + p)^2) . $ <eq:g0_VV-int-d>

*$Rop$-operation.* At one loop $Rbar G = G$. _Dirac algebra of (a)._ By @eq:phys-prop and @eq:sandwich the Dirac trace of @eq:g0_VV-int-a is $tr[gbar^mu slb(k) gbar^nu (slb(k) + slb(p)) P_(c_f)]$, a purely four-dimensional trace. The $g5$ part is $plus.minus 2 i epsilon^(mu alpha nu beta) bar(k)_alpha bar(p)_beta$ (@eq:tr5); after integration $bar(k)_alpha$ becomes a multiple of $bar(p)_alpha$ and it vanishes. With $i^4 = 1$,
$ G_"(a)" = - 2 S_f^((a)) delta^(a b) intk (2 bar(k)^mu bar(k)^nu + bar(k)^mu bar(p)^nu + bar(p)^mu bar(k)^nu - gmb^(mu nu) (bar(k)^2 + bar(k) dot bar(p)))/(k^2 (k + p)^2) . $ <eq:g0_VV-num-a>
The numerator contains no hatted object; the denominators are $d$-dimensional, $k dot p = bar(k) dot bar(p) + hat(k) dot hat(p)$. _Numerators of the bosonic loops._ With $i^4 = 1$ in (b1) and $i^2 = -1$ in (b2),
$ G_"(b1)" = S_H^((a)) delta^(a b) intk ((2 k + p)^mu (2 k + p)^nu)/((k^2 - m^2)((k + p)^2 - m^2)), quad
  G_"(b2)" = - 2 S_H^((a)) delta^(a b) g^(mu nu) intk 1/(k^2 - m^2) . $ <eq:g0_VV-num-b>
In (c1), with $u = p - k$, $v = 2 k + p$, $w = -(k + 2 p)$ the three momenta of $cal(V)^(mu rho sigma) = g^(mu rho) u^sigma + g^(rho sigma) v^mu + g^(sigma mu) w^rho$ and $g^rho_(space rho) = d$,
$ cal(V)^(mu rho sigma) cal(V)^nu_(space rho sigma) = g^(mu nu) (u^2 + w^2) + d thin v^mu v^nu + u^mu v^nu + v^mu u^nu + u^mu w^nu + w^mu u^nu + v^mu w^nu + w^mu v^nu . $ <eq:g0_VV-num-c1>
All bosonic numerators are $d$-dimensional and never split.

_Expansion._ With @eq:uv-expansion and $D = 1\/(k^2 - MUV^2)$,
$ 1/k^2 = D (1 - MUV^2 D + dots), \
  1/((k + p)^2 - m^2) = D (1 + (-2 k dot p - p^2 + m^2 - MUV^2) D + 4 (k dot p)^2 D^2 + dots) , $ <eq:g0_VV-T>
and $Top(2)$ keeps the terms of total degree $<= 2$ in $p$, $m$ and $MUV$ (numerator included). The vacuum poles needed are those of @eq:logpole and below it: $intk (k^2)^a D^(a+2) -> ilf 1\/ep$ and $intk (k^2)^a D^(a+1) -> ilf (a + 1) MUV^2 \/ ep$. Tensor reduction @eq:tensor-red, with $d -> 4$ in its coefficients for $N_0$, gives the bubble integrals
$ Kop Top(2) intk {1, thin k^alpha, thin k^alpha k^beta}/((k^2 - m^2)((k + p)^2 - m^2))
  = ilf 1/ep {1, thin -1/2 p^alpha, thin 1/3 p^alpha p^beta - 1/12 p^2 g^(alpha beta) + 1/2 m^2 g^(alpha beta)} , $ <eq:g0_VV-masters>
$ Kop Top(2) intk 1/(k^2 - m^2) = ilf m^2/ep, quad quad Kop Top(2) intk 1/(k^2) = 0 . $ <eq:g0_VV-tad>
In the tensor integral the degree-$0$ term $k^alpha k^beta D^2$ gives $+1/2 MUV^2 g^(alpha beta)$ and the term $-MUV^2 D$ of @eq:g0_VV-T (from either line) gives $-1/2 MUV^2 g^(alpha beta)$; in the tadpoles $D$ gives $MUV^2$ and $(m^2 - MUV^2) D^2$ gives $m^2 - MUV^2$. $MUV$ cancels in each topology separately. For the split numerator of (a), the projections @eq:tensor-red-split of @eq:g0_VV-masters at $m = 0$ are
$ bar(k)^mu bar(k)^nu -> 1/3 bar(p)^mu bar(p)^nu - 1/12 p^2 gmb^(mu nu), quad
  bar(k)^mu -> -1/2 bar(p)^mu, quad
  bar(k)^2 = gmb_(alpha beta) k^alpha k^beta -> 1/3 bar(p)^2 - 1/3 p^2 = -1/3 hat(p)^2 , $ <eq:g0_VV-split>
with $gmb_(alpha beta) g^(alpha beta) = 4$ and $p^2 = bar(p)^2 + hat(p)^2$ $d$-dimensional: the hatted part of $p^2$ comes from $hat(k) dot hat(p)$ in the denominator, paired with another $hat(k) dot hat(p)$. No trace $gmh^mu_(space mu)$ occurs in any topology, so $N_1 = 0$ throughout (@sec:m-K).

_Pole residues._ Inserting @eq:g0_VV-masters and @eq:g0_VV-split in @eq:g0_VV-num-a to @eq:g0_VV-num-c1, every topology takes the form $Top(2) G = ilf N_0 \/ ep + dots$ with
$ N_0^"(a)" = - S_f^((a)) delta^(a b) [ 2/3 (bar(p)^2 gmb^(mu nu) - bar(p)^mu bar(p)^nu) + 1/3 hat(p)^2 gmb^(mu nu) ], $ <eq:g0_VV-N-a>
$ N_0^"(b1)" = S_H^((a)) delta^(a b) [ -1/3 (p^2 g^(mu nu) - p^mu p^nu) + 2 m^2 g^(mu nu) ], quad N_0^"(b2)" = - 2 m^2 S_H^((a)) delta^(a b) g^(mu nu), $ <eq:g0_VV-N-b>
$ N_0^"(c1)" = C_"ad"^((a)) delta^(a b) (19/12 p^2 g^(mu nu) - 11/6 p^mu p^nu), quad N_0^"(c2)" = 0, quad
  N_0^"(d)" = C_"ad"^((a)) delta^(a b) (1/12 p^2 g^(mu nu) + 1/6 p^mu p^nu) . $ <eq:g0_VV-N-cd>
For (a): $2 dot 2 (1/3 bar(p) bar(p) - 1/12 p^2 gmb) - 2 dot 2 dot 1/2 bar(p) bar(p) - 2 gmb (-1/3 hat(p)^2 - 1/2 bar(p)^2)$ is the bracket of @eq:g0_VV-N-a after $p^2 = bar(p)^2 + hat(p)^2$. For (b1): $4 (1/3 p p - 1/12 p^2 g + 1/2 m^2 g) - 2 p p + p p$. For (d): $-(1/3 p p - 1/12 p^2 g - 1/2 p p)$. For (c1), @eq:g0_VV-masters at $m = 0$ (with $d -> 4$ in the factor $d$ of $d thin v v$) gives term by term
$ g^(mu nu) (u^2 + w^2) -> 4 p^2 g^(mu nu), quad d thin v^mu v^nu -> 4/3 p^mu p^nu - 4/3 p^2 g^(mu nu), quad u^mu v^nu + v^mu u^nu -> -1/3 p^mu p^nu + 1/3 p^2 g^(mu nu), \
  u^mu w^nu + w^mu u^nu -> -13/3 p^mu p^nu - 1/6 p^2 g^(mu nu), quad v^mu w^nu + w^mu v^nu -> -1/3 p^mu p^nu + 1/3 p^2 g^(mu nu) , $ <eq:g0_VV-c1-terms>
where $u^2 + w^2 = 2 k^2 + 2 k dot p + 5 p^2$ and $intk k^2 \/ (k^2 (k + p)^2)$ has no pole (@eq:g0_VV-tad). The sum is $19/6 p^2 g^(mu nu) - 11/3 p^mu p^nu$, times $1\/2$.

*Result.* Let $a, b$ lie in one gauge factor $r$ (the mixed functions vanish), and write $S^((r)) = S_"lh"^((r)) + S_"rh"^((r)) = sum_f S_f^((r))$ (@eq:SLR) for the sum over the five multiplets, as in @sec:st-solution.
#result[
  $ Gamma^((1) mu nu)_(A^a A^b) (p) |_"pole" = lf 1/ep N_0, quad
    N_0 = & [5/3 fg^(a c d) fg^(b c d) - 1/3 Tr(Tg^a_H Tg^b_H)] (p^2 g^(mu nu) - p^mu p^nu) \
    & - sum_f Tr_f (Tg^a_f Tg^b_f) [ 2/3 (bar(p)^2 gmb^(mu nu) - bar(p)^mu bar(p)^nu) + 1/3 hat(p)^2 gmb^(mu nu) ] $ <eq:g0_VV-pole>
  with $fg^(a c d) fg^(b c d) = C_"ad"^((a)) delta^(a b)$, $Tr(Tg^a_H Tg^b_H) = S_H^((a)) delta^(a b)$, $sum_f Tr_f (Tg^a_f Tg^b_f) = S^((a)) delta^(a b)$. The counterterm is minus the local functional whose vertex function is $lf N_0 \/ ep$ (quadratic parts; $F$ $d$-dimensional, $bar(F)$ with $bar(A)$ and barred indices):
  $ delta S^((1)) supset - lf 1/ep integral dif^d x thin [ & -(5/3 C_"ad"^((r)) - 1/3 S_H^((r))) (-1/4 F^a_(mu nu) F^(a mu nu)) \
    & + 2/3 S^((r)) (-1/4 bar(F)^a_(mu nu) bar(F)^(a mu nu)) + 1/3 S^((r)) thin 1/2 bar(A)^a_mu hat(partial)^2 bar(A)^(a mu) ] , $ <eq:g0_VV-ct>
  summed over $a$ in the factor $r$. The bosonic loops multiply the $d$-dimensional operator with unsplit indices, the fermion loop the barred operator. With $Q^((r))_T$ continued with unsplit $d$-dimensional indices (@sec:evanescent), $-1/4 F F$ is $Q^((r))_T$ at quadratic order, and $-1/4 bar(F) bar(F) = -1/4 F F + tilde(E)^((r))_F$ by @eq:evanescent. Hence
  $ delta S^((1)) supset delta Z^((r))_A thin Q^((r))_T - lf 1/ep [ 2/3 S^((r)) thin tilde(E)^((r))_F |_(A A) + 1/3 S^((r)) thin E^((r))_A ], $ <eq:g0_VV-basis>
  $ 16 pi^2 ep thin delta Z^((r))_A = 5/3 C_"ad"^((r)) - 1/3 S_H^((r)) - 2/3 S^((r)) , $ <eq:g0_VV-zA>
  where $tilde(E)^((r))_F |_(A A)$ is the part of $tilde(E)^((r))_F$ quadratic in the gauge field. $delta Z^((r))_A$ is the gauge-field renormalisation constant of the table of constants in @sec:st-solution ($A_0 = (Z^((r))_A)^(1\/2) A$, $Z^((r))_A = 1 + delta Z^((r))_A$). The coefficient $-2/3 S^((r))$ of $tilde(E)^((r))_F$ (in units of $1\/(16 pi^2 ep)$) is the same in @sec:g0_VVV and @sec:g0_VVVV, and @eq:g0_VV-basis is the two-point part of the assembled counterterm. $N_1 = 0$. In terms of the tree vertex function $Gamma^((0))_(A A)$ and its barred form $overline(Gamma)^((0))_(A A) = -(bar(p)^2 gmb^(mu nu) - bar(p)^mu bar(p)^nu) delta^(a b)$, the physical part of the pole is $N_0 = z^((r))_"bos" Gamma^((0))_(A A) + z^((r))_"ferm" overline(Gamma)^((0))_(A A)$ plus the evanescent term $-1/3 S^((r)) hat(p)^2 gmb^(mu nu) delta^(a b)$, with $z^((r))_"bos" = -5/3 C_"ad"^((r)) + 1/3 S_H^((r))$ (bosonic loops) and $z^((r))_"ferm" = 2/3 S^((r))$ (fermion loops). The same two names are used in @sec:g0_VVV and @sec:g0_VVVV, with $z^((r))_"ferm" = 2/3 S^((r))$ in all three classes.
]
The vertex function of $-1/4 F F$ is $-(p^2 g - p p) delta^(a b)$, that of $-1/4 bar(F) bar(F)$ is $-(bar(p)^2 gmb - bar(p) bar(p)) delta^(a b)$, that of $tilde(E)^((r))_F$ is the difference $[(p^2 g - p p) - (bar(p)^2 gmb - bar(p) bar(p))] delta^(a b) = (p^2 gmh^(mu nu) + hat(p)^2 gmb^(mu nu) - bar(p)^mu hat(p)^nu - hat(p)^mu bar(p)^nu - hat(p)^mu hat(p)^nu) delta^(a b)$, and that of $E^((r))_A$ is $-hat(p)^2 gmb^(mu nu) delta^(a b)$ ($hat(partial)^2 -> -hat(p)^2$, two assignments times $1\/2$).

#figure(
  table(
    columns: 6,
    align: (left, center, center, center, center, center),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$r$], [gauge + ghost \ $5/3 C_"ad"^((r))$], [Higgs \ $-1/3 S_H^((r))$], [fermions, four-dim. \ $-2/3 sum_f S_f^((r))$], [evanescent $hat(p)^2 gmb$ \ $-1/3 sum_f S_f^((r))$], [$16 pi^2 ep thin delta Z^((r))_A$]),
    [$G$], [$5 gs^2$], [$0$], [$-4/3 gs^2$], [$-2/3 gs^2$], [$11/3 gs^2$],
    [$W$], [$10/3 gw^2$], [$-1/6 gw^2$], [$-4/3 gw^2$], [$-2/3 gw^2$], [$11/6 gw^2$],
    [$B$], [$0$], [$-1/6 gy^2$], [$-20/9 gy^2$], [$-10/9 gy^2$], [$-43/18 gy^2$],
    [mixed], [$0$], [$0$], [$0$], [$0$], [$0$],
  ),
  caption: [Standard-Model values of the coefficients of @eq:g0_VV-pole, in units of $delta^(a b) \/ (16 pi^2 ep)$. The gauge-plus-ghost and Higgs columns multiply the $d$-dimensional $(p^2 g^(mu nu) - p^mu p^nu)$, the four-dimensional fermion column multiplies $(bar(p)^2 gmb^(mu nu) - bar(p)^mu bar(p)^nu)$, the evanescent column multiplies $hat(p)^2 gmb^(mu nu)$; the last column is the sum of the gauge-plus-ghost, Higgs and four-dimensional fermion columns. Per multiplet the fermion entries are $-2/3 S_f^((r))$ and $-1/3 S_f^((r))$ with $S_f^((r))$ of @tab:invariants: for $G$, $Q$: $gs^2$, $U$, $D$: $gs^2\/2$ each; for $W$, $Q$: $3 gw^2\/2$, $L$: $gw^2\/2$; for $B$, $Q$: $gy^2\/6$, $U$: $4 gy^2\/3$, $D$: $gy^2\/3$, $L$: $gy^2\/2$, $E$: $gy^2$. In @eq:g0_VV-basis the coefficient of $tilde(E)^((r))_F |_(A A)$ is the four-dimensional fermion column and that of $E^((r))_A$ the evanescent column.],
) <tab:g0_VV-sm>

#check[*Auxiliary mass and $m^2$.* $MUV$ drops out of $Kop Top(2)$ of each topology separately (@eq:g0_VV-masters, @eq:g0_VV-tad). The Higgs mass enters only through $+2 m^2 S_H g^(mu nu)$ in (b1) and $-2 m^2 S_H g^(mu nu)$ in (b2); their sum vanishes, so no $Q^((r))_M$ (gauge-boson mass) counterterm is generated. The massless tadpole (c2) has no pole.]

#check[*Transversality.* The bosonic part of @eq:g0_VV-pole is transverse in $d$ dimensions and the four-dimensional fermion part is transverse, $p_mu (bar(p)^2 gmb^(mu nu) - bar(p)^mu bar(p)^nu) = 0$; so no $Q^((r))_L$ counterterm arises, as required by @eq:linear-ids. The only non-transverse pole is the evanescent one: $-i p_mu Gamma^((1) mu nu)_(A^a A^b)|_"pole" = i/3 S^((a)) hat(p)^2 bar(p)^nu delta^(a b) \/ (16 pi^2 ep)$. By @eq:m-wA-1loop and @eq:m-QAP this is the pole of $(Dbrk_1 dot Gamma)^((1))_(gh^a (p) A^b_nu (-p))$, to be found in @sec:g1_cV. It is cancelled by $bop_d$ of the $E^((r))_A$ counterterm: with @eq:b-evanescent the vertex function of $bop_d E^((r))_A$ on $gh^a (p) A^b_nu (-p)$ is $i hat(p)^2 bar(p)^nu delta^(a b)$, and its coefficient in @eq:g0_VV-basis is $-1/3 S^((r)) \/ (16 pi^2 ep)$. $tilde(E)^((r))_F$ is $bop_d$-closed (@sec:evanescent) and does not enter.]

#check[*Known limits.* Per Weyl fermion the four-dimensional coefficient $-2/3 S_f$ is half the Dirac value $-4/3 T_F$ ($T_F$ the index of the representation), the scalar value $-1/3 S_H$ is that of one complex scalar, and gauge plus ghost give the Feynman-gauge value $5/3 C_"ad"$ (with $19/12, -11/6$ and $1/12, 1/6$ for (c1) and (d) separately). For the abelian factor $Z_(c_B) (Z^((B))_A)^(1\/2) = 1$, with $Z_(c_B)$ the renormalisation constant of $c_B = gy$ (@sec:st-solution), so $-16 pi^2 ep thin delta Z^((B))_A = 43/18 gy^2$ gives the one-loop running $16 pi^2 thin dif gy \/ dif ln mu = 43/18 thin gy^3$ for one family: $20/9$ per family from fermions and $1/6$ from $H$, i.e. the known coefficient $41/6$ for three families. The evanescent pole $bar(A) hat(partial)^2 bar(A)$, proportional to the fermion index and without bosonic contribution, has the form found for chiral Yang--Mills theory in @Belusca-Maito:2020ala. The assembled counterterm, with $delta Z^((r))_A$ and the constants of @sec:g0_VVV, @sec:g0_VVVV and @sec:g0_ghost_sources, is given in @sec:st-solution.]

#general[At one loop the forest is empty and the evanescent pole is a pure $hat(p)^2$ term of the external momentum. From two loops on, the counterterm $-1/3 S^((r)) E^((r))_A \/ (16 pi^2 ep)$ is inserted with its hatted structure into the one-loop subgraphs; there $hat(partial)^2 -> -hat(k)^2$ gives $eta \/ ep^2$ terms, i.e. both a pole and an evanescent-induced finite part of the two-loop gauge two-point function (@sec:m-beyond).]
