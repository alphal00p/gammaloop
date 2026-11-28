#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Breaking of $omega A A A$ and $omega A A A A$ <sec:g1_cVVV_cVVVV>

External legs $gh^a (q)$, $A^(b_j)_(mu_j) (p_j)$ with $j = 1, 2, 3$ (class $gh A A A$) or $j = 1, dots, 4$ (class $gh A A A A$), all momenta incoming, $q + sum_j p_j = 0$. The sum of the diagrams below is $i thin (Dbrk_1 dot Gamma)^((1))$ for these legs. Fermions couple only to $bar(A)$, so all $mu_j$ are barred; the momenta are $d$-dimensional. Degrees $dg(G) = 1$ ($gh A A A$) and $dg(G) = 0$ ($gh A A A A$). All gauge factors and all combinations of them are covered, for the five multiplets $Q, U, D, L, E$. Local symbol: $s_f = -1$ for $f = Q, L$ and $s_f = +1$ for $f = U, D, E$, so that $P_(c_f) = (1 + s_f g5)\/2$.

*Diagrams.* The $Dbrk_1$ vertex @eq:V-delta has two fermion legs. In a one-loop 1PI graph without external fermions both are internal, so the loop is a closed $Psi_f$ line through the $Dbrk_1$ vertex and through one gauge vertex @eq:V-ffv per external gauge leg. Graphs with a bosonic self-coupling attached to the loop are not 1PI. There are two topologies, the box (a) and the pentagon (b) of @fig:g1_cVVV_cVVVV-diagrams. A graph is fixed by the order in which the fermion line meets the gauge legs after leaving $Dbrk_1$; the two orientations of the loop are different orderings. Each term of @eq:V-delta puts the sterile component on one of the two lines adjacent to the vertex.

#figure(
  diagrow(
    oneloop(
      ((legs: (("G", [$gh^a (q)$#h(5mm)], 1),), mark: "delta"), (legs: (("V", $A^(b_z)$, 0),)), (legs: (("V", [#h(3mm)$A^(b_y)$], 0),)), (legs: (("V", $A^(b_x)$, 0),))),
      (("F", $k$, 1), ("F", none, 1), ("F", none, 1), ("F", [$k - q$#h(4mm)], 1)), r: 0.75, leg: 0.8, caption: [(a) box, $gh A A A$]),
    oneloop(
      ((legs: (("G", [$gh^a (q)$#h(5mm)], 1),), mark: "delta"), (legs: (("V", $A^(b_w)$, 0),)), (legs: (("V", $A^(b_z)$, 0),)), (legs: (("V", $A^(b_y)$, 0),)), (legs: (("V", $A^(b_x)$, 0),))),
      (("F", $k$, 1), ("F", none, 1), ("F", none, 1), ("F", none, 1), ("F", [$k - q$#h(4mm)], 1)), r: 0.8, leg: 0.8, caption: [(b) pentagon, $gh A A A A$]),
  ),
  caption: [Topologies of the class. The crossed circle is the vertex @eq:V-delta. The fermion line leaves it with momentum $k$ and returns with $k - q$; the line next to the vertex carries the sterile component in one of the two terms of @eq:V-delta.],
) <fig:g1_cVVV_cVVVV-diagrams>

#figure(
  {
  set par(justify: false)
  table(
    columns: (auto, 1fr, auto, auto, auto),
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [sym.], [sign], [group factor]),
    [(a)], [$Psi_f$, $f in {Q, U, D, L, E}$; 6 graphs, one per permutation $(x, y, z)$ of $(1, 2, 3)$], [$1$], [$-1$ (fermion loop)], [$Tr_f (Tg^a Tg^(b_x) Tg^(b_y) Tg^(b_z))$],
    [(b)], [$Psi_f$, $f in {Q, U, D, L, E}$; 24 graphs, one per permutation $(w, x, y, z)$ of $(1, 2, 3, 4)$], [$1$], [$-1$ (fermion loop)], [$Tr_f (Tg^a Tg^(b_w) Tg^(b_x) Tg^(b_y) Tg^(b_z))$],
  )},
  caption: [Field assignments; "sym." is the symmetry factor. The permutation lists the gauge legs in the order of the Dirac string, i.e. against the fermion flow starting at $Dbrk_1$: in (a) the line leaving $Dbrk_1$ meets $z$ first, then $y$, then $x$. All $Tg$ carry the label $f$.],
) <tab:g1_cVVV_cVVVV-diagrams>

*Integrands.* The loop momentum $k$ runs along the fermion arrow on the line from $Dbrk_1$ to the first gauge vertex; then $k + p_z$, $k + p_z + p_y$, and $k - q$ on the line that returns to $Dbrk_1$. In @eq:V-delta, $p = k - q$ and $p' = k$. Write $F_j = i gbar^(mu_j) Pc$ for the Dirac part of the Feynman rule @eq:V-ffv of the gauge leg $j$ (its generator $Tg^(b_j)$ goes into the group trace) and $S(l) = i sl(l) \/ l^2$ for the propagator. The generators commute with the Dirac matrices, so the group trace $Tr_f$ factorises from the Dirac trace $tr$. Writing the Dirac string against the fermion flow, with $Pc = P_(c_f)$,
$ G_"(a)"^(x y z) = & - Tr_f (Tg^a Tg^(b_x) Tg^(b_y) Tg^(b_z)) \
  & times intk tr[ ((slh(k) - slh(q)) Pcb - slh(k) Pc) S(k - q) thin F_x thin S(k + p_y + p_z) thin F_y thin S(k + p_z) thin F_z thin S(k) ] . $ <eq:g1_cVVV_cVVVV-int-a>
$G_"(b)"^(w x y z)$ is the same with one more factor $F S$ and the momenta $k + p_z + p_y + p_x$, $k + p_z + p_y$, $k + p_z$.

_The sterile chain._ Write $i V_(Dbrk_1)(p', p) = Tg^a (slh(p) Pcb - slh(p)' Pc)$ for the Feynman rule @eq:V-delta. The vertex and its two neighbouring propagators always occur in the combination $Pc S(p') thin i V_(Dbrk_1) S(p) Pcb$: the gauge vertex to the left supplies $gbar Pc$, the one to the right $Pcb gbar$. By the two forms of @eq:m-sterile, the first term of @eq:V-delta turns $S(p)$ into $i hat(p)^2 \/ p^2$ and the second turns $S(p')$ into $-i hat(p)'^2 \/ p'^2$; the remaining propagator is physical by @eq:phys-prop. Hence
$ Pc thin S(p') thin i V_(Dbrk_1)(p', p) thin S(p) thin Pcb = Tg^a (hat(p)'^2 thin slb(p) - hat(p)^2 thin slb(p)')/(p^2 thin p'^2) thin Pcb . $ <eq:g1_cVVV_cVVVV-chain>
All other numerator factors are four-dimensional. Collecting the factors $i$ ($i^3$ from the gauge vertices, $i^2$ from the two physical propagators) and moving all projectors to the right,
$ G_"(a)"^(x y z) = - i thin Tr_f (Tg^a Tg^(b_x) Tg^(b_y) Tg^(b_z)) intk (tr[ N_"ch" thin gbar^(mu_x) (slb(k) + slb(p)_y + slb(p)_z) gbar^(mu_y) (slb(k) + slb(p)_z) gbar^(mu_z) Pc ])/(k^2 (k - q)^2 (k + p_y + p_z)^2 (k + p_z)^2) , $ <eq:g1_cVVV_cVVVV-Ga>
$ N_"ch" = hat(k)^2 thin (slb(k) - slb(q)) - (hat(k) - hat(q))^2 thin slb(k) = - hat(k)^2 thin slb(q) + (2 hat(k) dot hat(q) - hat(q)^2) thin slb(k) . $ <eq:g1_cVVV_cVVVV-Nch>
$G_"(b)"^(w x y z)$ has the same form with the string $gbar^(mu_w) (slb(k) + dots) gbar^(mu_x) dots gbar^(mu_z) Pc$, the group factor $Tr_f (Tg^a Tg^(b_w) Tg^(b_x) Tg^(b_y) Tg^(b_z))$ and five denominators.

_Zero ghost momentum._ $N_"ch"$ vanishes at $q = 0$, so $G_"(a)"$ and $G_"(b)"$ vanish identically at $q = 0$, for all $p_j$ and before any expansion. At the level of the vertex this is the Dirac identity
$ Pc thin sl(k) thin i V_(Dbrk_1)(k, k) thin sl(k) thin Pcb = 0, quad quad i V_(Dbrk_1)(k, k) = Tg^a slh(k) (Pcb - Pc) . $ <eq:g1_cVVV_cVVVV-q0>
Proof: with @eq:sandwich, $Pcb sl(k) Pcb = slh(k) Pcb$ and $Pc sl(k) Pc = slh(k) Pc$, and $slh(k)$ commutes with both projectors. The first term gives $Pc sl(k) slh(k) slh(k) Pcb = hat(k)^2 thin Pc sl(k) Pcb$, the second $Pc sl(k) slh(k) Pc sl(k) Pcb = slh(k) Pc slh(k) sl(k) Pcb = hat(k)^2 thin Pc sl(k) Pcb$; they cancel. The identity holds for any denominators, because both neighbouring lines carry the same momentum.

*$Rop$-operation.* At one loop $Rbar G = G$.

_Box, $dg(G) = 1$._ Numerator degree $5$ ($1$ from $Dbrk_1$, $4$ from the propagators), four propagators: $dg(G) = 4 + 5 - 8 = 1$. Take $q, p_y, p_z$ as independent momenta. Since $G_"(a)"$ vanishes at $q = 0$ for all $p_y, p_z$, every term of $Top(1) G_"(a)"$ is linear in $q$ and taken at $p_y = p_z = 0$. The factor $q$ sits in $N_"ch"$; the term $hat(q)^2 slb(k)$ has degree $2$ and is dropped, and the expansion @eq:uv-expansion of the denominators contributes only at higher degree, because it multiplies $N_"ch" |_(q = 0) = 0$. With $D_MUV = 1 \/ (k^2 - MUV^2)$,
$ Top(1) G_"(a)"^(x y z) = - i thin Tr_f (dots) intk D_MUV^4 thin tr[ (- hat(k)^2 slb(q) + 2 (hat(k) dot hat(q)) slb(k)) gbar^(mu_x) slb(k) gbar^(mu_y) slb(k) gbar^(mu_z) Pc ] . $ <eq:g1_cVVV_cVVVV-T1>
_Tensor reduction with the tag._ The second term contains $hat(k)^alpha bar(k)^beta bar(k)^gamma bar(k)^delta$; every pairing of @eq:tensor-red joins the hatted component with a barred one and vanishes (@eq:tensor-red-split). In the first term $hat(k)^2 bar(k)^beta bar(k)^gamma = gmh_(alpha alpha') k^alpha k^(alpha') bar(k)^beta bar(k)^gamma$, and of the three pairings only $g^(alpha alpha') g^(beta gamma)$ survives:
$ hat(k)^2 thin bar(k)^beta bar(k)^gamma -> ((k^2)^2)/(d (d + 2)) thin gmh^alpha_(space alpha) thin gmb^(beta gamma) = - (2 eta (k^2)^2)/(d (d + 2)) thin gmb^(beta gamma) . $ <eq:g1_cVVV_cVVVV-red>
Then $gbar_beta gbar^(mu_y) gbar^beta = -2 gbar^(mu_y)$ (@eq:contr2), and $intk (k^2)^2 D_MUV^4 -> ilf \/ ep$ (@eq:logpole); with $d (d + 2) -> 24$,
$ Top(1) G_"(a)"^(x y z) = - i (- 1)(- 2 eta)(- 2)/24 thin ilf 1/ep thin tr[slb(q) gbar^(mu_x) gbar^(mu_y) gbar^(mu_z) Pc] thin Tr_f (dots) + dots
  = ilf thin eta/ep thin N_1^("(a)" x y z) + dots . $ <eq:g1_cVVV_cVVVV-Ta>
The only integral is logarithmic and $MUV$-independent. Hence
$ N_0^("(a)" x y z) = 0, quad quad N_1^("(a)" x y z) = i/6 thin Tr_f (Tg^a Tg^(b_x) Tg^(b_y) Tg^(b_z)) thin tr[slb(q) gbar^(mu_x) gbar^(mu_y) gbar^(mu_z) Pc] . $ <eq:g1_cVVV_cVVVV-Na>
$N_0 = 0$ holds with $d$-dimensional external momenta: the only hatted external structure, $hat(q)^2$, has degree $2 > dg(G)$, and $hat(k) dot hat(q)$ is removed by the reduction. $N_1$ is four-dimensional, $Pphys N_1 = N_1$. The Dirac trace follows from @eq:tr5 and $Pc = (1 + s_f g5)\/2$:
$ tr[slb(q) gbar^(mu_x) gbar^(mu_y) gbar^(mu_z) Pc] = 2 (bar(q)^(mu_x) gmb^(mu_y mu_z) - bar(q)^(mu_y) gmb^(mu_x mu_z) + bar(q)^(mu_z) gmb^(mu_x mu_y)) - 2 i s_f thin epsilon^(alpha mu_x mu_y mu_z) bar(q)_alpha . $ <eq:g1_cVVV_cVVVV-tr>

_Pentagon, $dg(G) = 0$._ Numerator degree $6$, five propagators: $dg(G) = 4 + 6 - 10 = 0$. $Top(0)$ sets all external momenta to zero (no $MUV$ term has degree $0$), so $q = 0$ and the integrand is proportional to @eq:g1_cVVV_cVVVV-q0:
$ Top(0) G_"(b)"^(w x y z) = 0 quad "identically", quad quad N_0^"(b)" = N_1^"(b)" = 0 . $ <eq:g1_cVVV_cVVVV-Nb>
The non-local part of $G_"(b)"$ at $q != 0$ is finite and not used.

*Result.* $sum_f$ runs over $Q, U, D, L, E$, $sum_((x,y,z))$ over the six permutations of $(1, 2, 3)$, and $D^(a b c)$ is @eq:Dabc.
#result[
  $ (Dbrk_1 dot Gamma)^((1))_(gh^a A^(b_1)_(mu_1) A^(b_2)_(mu_2) A^(b_3)_(mu_3)) |_"pole" = lf 1/ep N_0, quad N_0 = 0; quad quad
    (Dbrk_1 dot Gamma)^((1))_(gh^a A^(b_1)_(mu_1) dots A^(b_4)_(mu_4)) |_"pole" = 0 . $ <eq:g1_cVVV_cVVVV-pole>
  $ & X_(1; gh^a A^(b_1)_(mu_1) A^(b_2)_(mu_2) A^(b_3)_(mu_3)) = lf thin Pphys N_1 , \
    & Pphys N_1 = sum_f sum_((x,y,z)) Tr_f (Tg^a Tg^(b_x) Tg^(b_y) Tg^(b_z)) { i/3 (bar(q)^(mu_x) gmb^(mu_y mu_z) - bar(q)^(mu_y) gmb^(mu_x mu_z) + bar(q)^(mu_z) gmb^(mu_x mu_y)) \
    & wide + s_f/3 thin epsilon^(alpha mu_x mu_y mu_z) bar(q)_alpha } , $ <eq:g1_cVVV_cVVVV-X>
  $ X_(1; gh^a A^(b_1)_(mu_1) dots A^(b_4)_(mu_4)) = 0 . $ <eq:g1_cVVV_cVVVV-X4>
  As a local functional (bars omitted, four-dimensional):
  $ X_1 supset lf sum_f integral dif^4 x thin { & -1/3 Tr_f [ (partial^mu wm) (Am_mu Am_nu Am^nu - Am_nu Am_mu Am^nu + Am_nu Am^nu Am_mu) ] \
    & + i/3 thin s_f thin epsilon^(mu nu rho sigma) Tr_f [ (partial_mu wm) Am_nu Am_rho Am_sigma ] } $ <eq:g1_cVVV_cVVVV-Xop>
  $ X_1 supset & lf integral dif^4 x thin [ - 1/3 thin d_4^(a b c d) (partial_mu gh^a A^(b mu)) (A^c dot A^d)
     + 1/12 thin D^(a e d) fg^(b c e) epsilon^(mu nu rho sigma) (partial_mu gh^a) A^b_nu A^c_rho A^d_sigma ] \
     & - lf 2/9 sum_(r = G, W) S^((r)) N^((r))_(f f) . $ <eq:g1_cVVV_cVVVV-Xop2>
]
Here $d_4^(a b c d) = 1/6 sum_f sum_("perm" (b' c' d') "of" (b c d)) Tr_f (Tg^a Tg^(b') Tg^(c') Tg^(d'))$ is totally symmetric (a symbol local to this subsection). The term $-2/9 S^((r)) N^((r))_(f f)$ is not the whole coefficient of $N^((r))_(f f)$ for $r = W$: the $d_4$ term contributes $+1/9 gw^2 N^((W))_(f f)$ (below @tab:g1_cVVV_cVVVV-sm), so the coefficient of $N^((r))_(f f)$ in $16 pi^2 X_1$ is $-2/9 S^((G)) = -4/9 gs^2$ and $-4/9 gw^2 + 1/9 gw^2 = -1/3 gw^2$. Translation: $partial_mu gh -> -i q_mu$ and the three gauge legs are distributed over the three positions in all ways, so
$ integral Tr_f [(partial^mu wm) Am_mu Am_nu Am^nu] -> -i sum_((x,y,z)) Tr_f (Tg^a Tg^(b_x) Tg^(b_y) Tg^(b_z)) thin bar(q)^(mu_x) gmb^(mu_y mu_z), \
  integral epsilon^(mu nu rho sigma) Tr_f [(partial_mu wm) Am_nu Am_rho Am_sigma] -> -i sum_((x,y,z)) Tr_f (Tg^a Tg^(b_x) Tg^(b_y) Tg^(b_z)) thin epsilon^(alpha mu_x mu_y mu_z) bar(q)_alpha , $
and matching to @eq:g1_cVVV_cVVVV-X gives the coefficients $-1/3$ and $i s_f \/ 3$ of @eq:g1_cVVV_cVVVV-Xop. For the second form, $Am_mu Am dot Am + Am dot Am Am_mu - Am_nu Am_mu Am^nu = 1/3 (Am_mu Am dot Am + Am_nu Am_mu Am^nu + Am dot Am Am_mu) + 2/3 [Am_nu, [Am^nu, Am_mu]]$; the symmetric part gives $d_4$, and $Tr_f [(partial^mu wm) [Am_nu, [Am^nu, Am_mu]]] = sum_r S^((r))_f N^((r))_(f f)$ (density). In the Levi-Civita term, $epsilon^(mu nu rho sigma) Am_nu Am_rho Am_sigma = 1/4 epsilon^(mu nu rho sigma) {[Am_nu, Am_rho], Am_sigma}$ and $sum_f s_f Tr_f (Tg^a {Tg^e, Tg^d}) = - D^(a e d)$.

In the Standard Model, the Levi-Civita term of @eq:g1_cVVV_cVVVV-Xop2 needs a non-abelian $fg^(b c e)$ and is proportional to $D^(a e d)$ (@eq:Dabc) with $e$ in the same factor. With $D_(G G G)$ and $D_(B r r)$ as in @sec:g1_cVV ($D^(A C E) = D_(G G G) gs^3 d^(A C E)$, $D^(a b B) = D_(B r r) c_r^2 gy delta^(a b)$ for $a, b$ in the factor $r$, $B$ the hypercharge index), $D_(G G G) = 1/2 [2 - 1 - 1]$, $D_(B G G) = 2 dot 1/6 - 2/3 + 1/3$ and $D_(B W W) = 3 dot 1/6 - 1/2$ (multiplets in the order $Q, U, D$ resp. $Q, L$, each weighted with its number of weak or colour components and with the sign of $sum_"lh" - sum_"rh"$); all vanish. Per multiplet the Levi-Civita term is non-zero for $Q, U, D, L$ and zero for $E$; only the sum vanishes. In the basis of @tab:Nbasis it is
$ sum_(r = G, W) D_(B r r; f) / 12 (N^((r))_(epsilon r B f) - N^((r))_(epsilon 4 B)) + D_(G G G; f) / 12 thin N_(epsilon d f) , $ <eq:g1_cVVV_cVVVV-eps-basis>
with $N_(epsilon d f)$ of @tab:Nbasis and $D_(B r r; f)$, $D_(G G G; f)$ the contributions of the multiplet $f$ to $D_(B r r)$, $D_(G G G)$: $D_(B G G; f) = (1/3, -2/3, 1/3)$ for $f = Q, U, D$, $D_(B W W; f) = (1/2, -1/2)$ for $f = Q, L$, $D_(G G G; f) = (1, -1/2, -1/2)$ for $f = Q, U, D$. Summed over $f$, these are the terms with $N^((r))_(epsilon r B f)$, $N^((r))_(epsilon 4 B)$ and $N_(epsilon d f)$ of the anomalous part $X_"anom"$ in @sec:st-solution (below @eq:sts-anom). $N_(epsilon d f)$ is the cVVV monomial of @tab:Nbasis that no $bop Q_alpha$ contains; it multiplies $D_(G G G)$ and has zero coefficient in the Standard Model. The parity-even part in the basis of @tab:Nbasis:

#figure(
  table(
    columns: 7,
    align: (left, center, center, center, center, center, center),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([monomial], [$Q$], [$U$], [$D$], [$L$], [$E$], [sum]),
    [$N_(G 4 a)$], [$-1\/36$], [$-1\/72$], [$-1\/72$], [$0$], [$0$], [$-1\/18$],
    [$N_(G 4 b)$], [$-1\/18$], [$-1\/36$], [$-1\/36$], [$0$], [$0$], [$-1\/9$],
    [$N^((G))_(f f)$], [$-2 gs^2\/9$], [$-gs^2\/9$], [$-gs^2\/9$], [$0$], [$0$], [$-4 gs^2\/9$],
    [$N_(W 4 a)$], [$-1\/8$], [$0$], [$0$], [$-1\/24$], [$0$], [$-1\/6$],
    [$N^((W))_(f f)$], [$-gw^2\/4$], [$0$], [$0$], [$-gw^2\/12$], [$0$], [$-gw^2\/3$],
    [$N_(B 4)$], [$-1\/648$], [$-16\/81$], [$-1\/81$], [$-1\/24$], [$-1\/3$], [$-95\/162$],
    [$N_(G G G B 1)$], [$-1\/18$], [$-1\/9$], [$1\/18$], [$0$], [$0$], [$-1\/9$],
    [$N_(G G G B 2)$, $N_(G G G B 3)$], [$-1\/36$], [$-1\/18$], [$1\/36$], [$0$], [$0$], [$-1\/18$],
    [$N^((G B))_a$, $N^((G B))_b$], [$-1\/108$], [$-2\/27$], [$-1\/54$], [$0$], [$0$], [$-11\/108$],
    [$N^((G B))_c$, $N^((G B))_d$], [$-1\/54$], [$-4\/27$], [$-1\/27$], [$0$], [$0$], [$-11\/54$],
    [$N^((W B))_a$, $N^((W B))_b$], [$-1\/72$], [$0$], [$0$], [$-1\/24$], [$0$], [$-1\/18$],
    [$N^((W B))_c$, $N^((W B))_d$], [$-1\/36$], [$0$], [$0$], [$-1\/12$], [$0$], [$-1\/9$],
    [$N^((G W))_a$, $N^((G W))_b$], [$-1\/12$], [$0$], [$0$], [$0$], [$0$], [$-1\/12$],
    [$N^((G W))_c$, $N^((G W))_d$], [$-1\/6$], [$0$], [$0$], [$0$], [$0$], [$-1\/6$],
  ),
  caption: [Coefficients of $16 pi^2 X_1$ in the class $gh A A A$ (parity-even part, @eq:g1_cVVV_cVVVV-Xop2) per multiplet. The couplings are inside the monomials of @tab:Nbasis except where shown. All other monomials of the class cVVV of @tab:Nbasis have coefficient zero; the Levi-Civita monomials have the per-multiplet coefficients of @eq:g1_cVVV_cVVVV-eps-basis and zero total coefficient. The class $gh A A A A$ has $X_1 = 0$.],
) <tab:g1_cVVV_cVVVV-sm>

The entries follow from $d_4$ and $S^((r))_f$ of @tab:invariants. For a colour triplet, from $Tr(T^A T^B T^C T^D) = delta^(A B) delta^(C D)\/12 + (d^(A B E) + i f^(A B E))(d^(C D E) + i f^(C D E))\/8$ and the $"SU"(3)$ identity of @sec:Qbasis,
$ 1/6 sum_"perm" Tr(T^A T^B T^C T^D) = 1/24 (delta^(A B) delta^(C D) + delta^(A C) delta^(B D) + delta^(A D) delta^(B C)) , $
and the same holds for $tau^I\/2$ on a doublet. Hence $d_4^(G G G G) = gs^4 (delta delta + delta delta + delta delta)\/6$ (four triplets), $d_4^(W W W W) = gw^4 (delta delta + delta delta + delta delta)\/6$ (four doublets), $d_4^(B B B B) = gy^4 sum_f N^c_f N^w_f Y_f^4 = 95 gy^4 \/ 54$ (with $B$ the hypercharge index and $N^c_f N^w_f$ the number of components of $f$, as in @tab:g1_cVV-sm), $d_4^(A C E B) = gs^3 gy d^(A C E) \/ 6$, $d_4^(A C B B) = 11 gs^2 gy^2 delta^(A C) \/ 36$, $d_4^(I J B B) = gw^2 gy^2 delta^(I J) \/ 6$, $d_4^(A C I J) = gs^2 gw^2 delta^(A C) delta^(I J) \/ 4$; components with one non-abelian index or three $W$ indices with one $B$ vanish. For $"SU"(2)$, $gw^4 (partial gh^I dot W^J)(W^I dot W^J) = N_(W 4 a) - gw^2 N^((W))_(f f)$ by $epsilon^(E I J) epsilon^(E K L) = delta^(I K) delta^(J L) - delta^(I L) delta^(J K)$.

#check[*Pole.* By the check of @sec:m-oneloop, $N_0 = 0$ requires $V_(bop_d delta S^((1)))[gh A A A] = V_(bop_d delta S^((1)))[gh A A A A] = 0$. The fermion-loop poles of @sec:g0_VVV and @sec:g0_VVVV carry no hatted structure at these degrees: a hatted external momentum needs two factors $hat(k) dot hat(p)$, of degree $2 > dg(G)$. Their counterterms are $lf 1/ep (-2/3 S^((r)))$ times the three- and four-field parts of the barred Yang--Mills action $-1/4 integral dif^d x thin bar(F)^a bar(F)^a$ ($Q^((r))_(3 f)$ and $cal(Q)^((r))_"YM"$ with barred fields), the same multiple as for its transverse two-field part from @sec:g0_VV; this functional is $bop_d$-closed. The bosonic loops have no $Dbrk_1$ vertex and $d$-dimensional BRST-covariant Feynman rules, so their poles are $bop_d$-closed in classes without fermion legs. The non-closed evanescent operator $E^((r))_A$ has $bop_d E^((r))_A$ only in the classes $gh A$ and $gh A A$ (@eq:b-evanescent). This is consistent with @eq:g1_cVVV_cVVVV-pole.]

#check[*Ghost momentum, $MUV$, chirality, anomaly.* Every term of @eq:g1_cVVV_cVVVV-Xop has a differentiated ghost, as required by @eq:g1_cVVV_cVVVV-q0; the monomial $N^((r))_(B f 3)$ with an undifferentiated $gh^B$ does not appear. $MUV$ does not enter: the only integral is $intk (k^2)^2 D_MUV^4$. The parity-even part is the same for left- and right-handed multiplets. The parity-odd part is proportional to $D^(a b c)$ after the sum over multiplets, as expected for the non-abelian anomaly $epsilon^(mu nu rho sigma) Tr[(partial_mu wm) Am_nu Am_rho Am_sigma]$, and vanishes for @tab:fields. The parity-even part is $bop$-exact within the quartic sector: $bop integral d_4^(a b c d) (A^a dot A^b)(A^c dot A^d) = 4 integral d_4^(a b c d) (partial gh^a dot A^b)(A^c dot A^d)$ by rigid invariance of $d_4$ (only @eq:dict-i contributes), and $bop (Q^((r))_T + Q^((r))_(3 f)) = N^((r))_(f f)$ by @eq:btable-gauge. The cubic parts of the corresponding $Q_alpha$ are fixed together with @sec:g1_cVV in @sec:st-solution.]

#check[*Numerical.* @eq:g1_cVVV_cVVVV-chain, @eq:g1_cVVV_cVVVV-q0 and the reduction @eq:g1_cVVV_cVVVV-T1 to @eq:g1_cVVV_cVVVV-Ta were checked with explicit $8 times 8$ and $16 times 16$ representations of the Clifford algebra in the integer dimensions $d = 6, 8$ (four barred, $d - 4$ hatted directions), both chiralities: the isotropic average of the full numerator of @eq:g1_cVVV_cVVVV-int-a, linear part in $q$, equals $2 (d - 4)\/(d (d + 2)) (k^2)^2 tr[slb(q) gbar^(mu_x) gbar^(mu_y) gbar^(mu_z) Pc]$. The sum over orderings and multiplets, @eq:g1_cVVV_cVVVV-Xop, @eq:g1_cVVV_cVVVV-Xop2 and @tab:g1_cVVV_cVVVV-sm were checked against the vertex functions of the monomials with explicit $"SU"(3) times "SU"(2) times "U"(1)$ generators and random kinematics; a least-squares fit on 33 of the 34 cVVV monomials of @tab:Nbasis (all except $N_(epsilon d f)$, whose coefficient $D_(G G G) \/ 12$ vanishes for this content; rank 33) reproduces the coefficients of @tab:g1_cVVV_cVVVV-sm, confirmed exactly with symbolic couplings, and the Levi-Civita part agrees with the $D^(a e d)$ form and with @eq:g1_cVVV_cVVVV-eps-basis multiplet by multiplet. Independently of the hand reduction, $N_0$ and $N_1$ of the complete box integrand @eq:g1_cVVV_cVVVV-int-a were computed with an independent computer-algebra implementation of the expansion and the pole extraction (@sec:supporting), with $d$-dimensional and with four-dimensional external momenta, both chiralities: $N_0 = 0$ and $Pphys N_1$ equal to @eq:g1_cVVV_cVVVV-Na; the complete pentagon gives $N_0 = N_1 = 0$.]

#check[*Literature.* For a single gauge group, integrating @eq:g1_cVVV_cVVVV-Xop by parts, the parity-even part is $1/6 lf integral gh^a bar(partial)_nu (bar(A)^b_mu bar(A)^(c mu) bar(A)^(d nu)) sum_f [Tr_f (Tg^a Tg^b Tg^c Tg^d) + (a c b d) + (a d b c) + (a d c b) - (a b d c) - (a c d b)]$, where $(a c b d)$ is the trace with that order of generators. This is the one-loop result of @Belusca-Maito:2020ala in our conventions; the Levi-Civita part is proportional to the anomaly trace in both.]

#general[At constant ghost, $Dbrk_1$ is the rotation of the sterile components alone. For constant real parameters $theta^a$ let $delta_chi (theta) Psi_f = i theta^a Tg^a_f Pcb Psi_f$, $delta_chi (theta) overline(Psi)_f = -i theta^a overline(Psi)_f Pc Tg^a_f$ (only $chi_f$, $overline(chi)_f$ move). In @eq:Lferm only the kinetic term varies, with $gamma^mu Pcb - Pc gamma^mu = ghat^mu (Pcb - Pc)$, so $delta_chi (theta) S_0 = - theta^a sum_f integral overline(Psi)_f Tg^a ghat^mu (Pcb - Pc) partial_mu Psi_f$; @eq:Delta1 at constant $gh$, after integration by parts, is the same with $theta -> gh$ and the opposite sign:
$ Dbrk_1 |_(partial gh = 0) = - delta_chi (gh) thin S_0 . $
In momentum space this is @eq:g1_cVVV_cVVVV-q0: $i V_(Dbrk_1)(k, k)$ is the sterile rotation of the kinetic vertex, which cancels both neighbouring propagators and leaves $i Tg^a (Pcb S(k) - S(k) Pc)$, projected out by the physical vertices. $delta_chi$ is linear in the quantum fields, and the change of variables $Psi -> Psi + delta_chi Psi$ has unit Jacobian in dimensional regularisation. If $S_"ren" - S_0$ contains no sterile component (true for $delta S^((1))$ and for the four-dimensional $Delta S$ of @sec:m-objects; the source terms couple to $sop Psi_f prop Pc Psi_f$ only), then to all orders
$ (Dbrk_1 dot Gamma) |_(partial gh = 0) = - delta_chi (gh) thin Gamma = - integral [ (delta_chi Psi_f) (delta Gamma)/(delta Psi_f) + (delta_chi overline(Psi)_f) (delta Gamma)/(delta overline(Psi)_f) ] , $
so every 1PI function $(Dbrk_1 dot Gamma)_(gh^a (q) Phi_1 dots Phi_n)$ without external $Psi$, $overline(Psi)$ vanishes at $q = 0$ at every loop order. Consequences: in classes with only gauge (and Higgs, ghost) legs the insertion of $Dbrk_1$ is $O(q)$; a component of degree $0$, such as $gh A A A A$, has no local part from the $Dbrk_1$ insertion at any order. At one loop this is the vanishing of $Top(0) G_"(b)"$; at higher loops it applies to the graphs of class (a) of @sec:m-beyond, while the other vertices of $Dbrk^([l])$ ($bop_d delta S^((L))$, $bop_d Delta S^((L))$) are not of this form and must be inserted separately.]
