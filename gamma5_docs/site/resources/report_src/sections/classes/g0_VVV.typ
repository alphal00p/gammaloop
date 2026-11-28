#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Gauge three-point functions $A A A$ <sec:g0_VVV>

External legs $A^a_mu (p_1) thin A^b_nu (p_2) thin A^c_rho (p_3)$, all momenta incoming, $p_1 + p_2 + p_3 = 0$; $a, b, c$ run independently over the twelve gauge indices, so that $G G G$, $W W W$, $B B B$ and every mixed combination are covered at once. The sum of the diagrams below is $i thin Gamma^((1))_(A^a_mu A^b_nu A^c_rho)$. Every diagram has superficial degree of divergence $dg(G) = 1$. External momenta and indices are $d$-dimensional. Three tensors are used throughout (local to this subsection):
$ cal(V)^(lambda_1 lambda_2 lambda_3)(q_1, q_2, q_3) = g^(lambda_1 lambda_2) (q_1 - q_2)^(lambda_3) + g^(lambda_2 lambda_3) (q_2 - q_3)^(lambda_1) + g^(lambda_3 lambda_1) (q_3 - q_1)^(lambda_2), \
  cal(V)^(mu nu rho) = cal(V)^(mu nu rho)(p_1, p_2, p_3), $ <eq:g0_VVV-V>
$ Xi^(mu nu rho)(w) = g^(mu nu) w^rho + g^(nu rho) w^mu + g^(rho mu) w^nu . $ <eq:g0_VVV-Xi>
$cal(V)$ is the bracket of @eq:V-vvv; written without arguments it is taken at the external momenta, @eq:g0_VVV-V. The tree vertex function is $Gamma^((0))_(A A A) = -i fg^(a b c) cal(V)^(mu nu rho)$ (the Feynman rule @eq:V-vvv divided by $i$). $overline(cal(V))^(mu nu rho)$ and $overline(Xi)^(mu nu rho)$ denote the same tensors with $gmb$ and $bar(p)_j$. $cal(V)^(mu nu rho)$ is totally antisymmetric under permutations of the pairs $(mu, p_1)$, $(nu, p_2)$, $(rho, p_3)$, and so is $fg^(a b c)$.

*Diagrams.* The one-loop 1PI topologies with three external gauge fields are drawn in @fig:g0_VVV-diagrams. There are no others: fermions and ghosts have only cubic gauge vertices; a tadpole would need a vertex with five legs; the Yukawa and quartic scalar couplings do not enter at this order.
#figure(
  stack(dir: ttb, spacing: 4mm,
  diagrow(
    oneloop(
      ((legs: (("V", $A^a_mu$, 0),)), (legs: (("V", $A^b_nu$, 0),)), (legs: (("V", $A^c_rho$, 0),))),
      (("F", none, -1), ("F", none, -1), ("F", none, -1)), scale: 0.85, caption: [(a) $Psi_f$ triangle]),
    oneloop(
      ((legs: (("V", $A^a_mu$, 0),)), (legs: (("V", $A^b_nu$, 0),)), (legs: (("V", $A^c_rho$, 0),))),
      (("S", none, 1), ("S", none, 1), ("S", none, 1)), scale: 0.85, caption: [(b) $H$ triangle]),
    oneloop(
      ((legs: (("V", $A^a_mu$, 0), ("V", $A^b_nu$, 0))), (legs: (("V", $A^c_rho$, 0),))),
      (("S", none, 1), ("S", none, 1)), scale: 0.85, caption: [(c) $H$ seagull bubble]),
  ),
  diagrow(
    oneloop(
      ((legs: (("V", $A^a_mu$, 0),)), (legs: (("V", $A^b_nu$, 0),)), (legs: (("V", $A^c_rho$, 0),))),
      (("V", none, 0), ("V", none, 0), ("V", none, 0)), scale: 0.85, caption: [(d) gauge triangle]),
    oneloop(
      ((legs: (("V", $A^a_mu$, 0),)), (legs: (("V", $A^b_nu$, 0), ("V", $A^c_rho$, 0)))),
      (("V", none, 0), ("V", none, 0)), scale: 0.85, caption: [(e) gauge bubble]),
    oneloop(
      ((legs: (("V", $A^a_mu$, 0),)), (legs: (("V", $A^b_nu$, 0),)), (legs: (("V", $A^c_rho$, 0),))),
      (("G", none, 1), ("G", none, 1), ("G", none, 1)), scale: 0.85, caption: [(f) ghost triangle]),
  )),
  caption: [Topologies of $A A A$, one orientation or channel each. (a), (b), (f) occur with both orientations of the charged line, (c) and (e) in the three channels $(a b | c)$, $(b c | a)$, $(c a | b)$, where the pair before the bar sits at the quartic vertex.],
) <fig:g0_VVV-diagrams>

#figure(
  table(
    columns: (auto, auto, auto, auto, 1fr),
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [lines], [sym.], [sign], [group factor]),
    [(a1), (a2)], [$Psi_f$, all $f$], [$1$], [$-1$ (fermion loop)], [$Tr_f [Tg^a Tg^b Tg^c]$, $thin Tr_f [Tg^a Tg^c Tg^b]$],
    [(b1), (b2)], [$H$], [$1$], [$+1$], [$Tr_H [Tg^a Tg^c Tg^b]$, $thin Tr_H [Tg^a Tg^b Tg^c]$],
    [(c), 3 channels], [$H$], [$1$], [$+1$], [$Tr_H [{Tg^a, Tg^b} Tg^c]$ and cyclic],
    [(d)], [$A$ (all twelve)], [$1$], [$+1$], [$fg^(a x y) fg^(b y z) fg^(c z x)$],
    [(e), 3 channels], [$A$], [$1\/2$], [$+1$], [$fg^(a x y) fg^(b c e) fg^(x y e)$, $thin fg^(a x y) fg^(b x e) fg^(c y e)$, $thin fg^(a x y) fg^(b y e) fg^(c x e)$; cyclic],
    [(f1), (f2)], [$gh$, $agh$], [$1$], [$-1$ (ghost loop)], [$fg^(y a x) fg^(z b y) fg^(x c z)$, and $b <-> c$],
  ),
  caption: [Field assignments, symmetry factors, signs and group factors. Repeated adjoint indices $x, y, z, e$ are summed over the twelve gauge indices.],
) <tab:g0_VVV-diagrams>

The group factors reduce with @eq:fg, @eq:invariants and the Jacobi identity. For $a, b, c$ in one non-abelian factor $r$, with $C_"ad" = C^((r))_"ad"$,
$ fg^(a x y) fg^(b y z) fg^(c z x) = fg^(y a x) fg^(z b y) fg^(x c z) = fg^(a x y) fg^(b x e) fg^(c y e) = - fg^(a x y) fg^(b y e) fg^(c x e) = C_"ad"/2 thin fg^(a b c), \
  fg^(a x y) fg^(b c e) fg^(x y e) = C_"ad" fg^(a b c), $ <eq:g0_VVV-group>
and all of them vanish unless $a, b, c$ belong to the same non-abelian factor. For the matter loops, $Tr_f [Tg^a [Tg^b, Tg^c]] = i fg^(a b c) S_f^((a))$ for every multiplet $f$ (including $H$).

*Integrands.* Bosonic triangles (b), (d), (f): the line from the vertex of $c$ to the vertex of $a$ carries $k$, the line from $a$ to $b$ carries $k + p_1$, the line from $b$ to $c$ carries $k - p_3$ (in the direction of the arrow for $H$ and $gh$). With @eq:V-ssv, @eq:V-vvv, @eq:V-ghost and @eq:propagators,
$ G_"(b1)" = - Tr_H [Tg^a Tg^c Tg^b] intk ((2k + p_1)^mu (2k + p_1 - p_3)^nu (2k - p_3)^rho)/([k^2 - m^2][(k+p_1)^2 - m^2][(k-p_3)^2 - m^2]) , $ <eq:g0_VVV-int-b>
$ G_"(d)" = i thin fg^(a x y) fg^(b y z) fg^(c z x) intk (cal(N)_"(d)"^(mu nu rho))/(k^2 (k+p_1)^2 (k-p_3)^2) , \
  cal(N)_"(d)"^(mu nu rho) = cal(V)^(mu alpha beta)(p_1, k, -k-p_1) thin g_(beta beta') thin cal(V)^(nu beta' gamma)(p_2, k+p_1, -k+p_3) thin g_(gamma gamma') thin cal(V)^(rho gamma' alpha')(p_3, k-p_3, -k) thin g_(alpha' alpha) , $ <eq:g0_VVV-int-d>
$ G_"(f1)" = - i thin fg^(y a x) fg^(z b y) fg^(x c z) intk ((k + p_1)^mu (k - p_3)^nu k^rho)/(k^2 (k+p_1)^2 (k-p_3)^2) . $ <eq:g0_VVV-int-f>
The prefactors are $i^3 i^3 = -1$ for (b1), $(-i)^3 = i$ for (d) and $(-1)(-1)^3 i^3 = -i$ for (f1). In (f1) the ghost line leaving the vertex of $a$ has colour $y$ and momentum $k + p_1$, which is the $p'$ of @eq:V-ghost. Fermion triangle (a1): the fermion line from the vertex of $b$ to that of $a$ carries $k$, from $c$ to $b$ carries $k - p_2$, from $a$ to $c$ carries $k + p_1$. With @eq:V-ffv and @eq:phys-prop all three numerators are four-dimensional and the projectors collect into one $Pcb$:
$ G_"(a1)" = (-1) thin i^6 thin Tr_f [Tg^a Tg^b Tg^c] intk (tr[gbar^mu slb(k) gbar^nu (slb(k) - slb(p)_2) gbar^rho (slb(k) + slb(p)_1) Pcb])/(k^2 (k - p_2)^2 (k + p_1)^2), quad (-1) thin i^6 = +1 . $ <eq:g0_VVV-int-a>
The second orientation of each triangle is obtained by $(b, nu, p_2) <-> (c, rho, p_3)$: $G_"(a2)"$, $G_"(b2)"$, $G_"(f2)"$. Seagull bubble (c), channel $(a b|c)$, $H$ momentum $k$ into the vertex of $c$:
$ G_"(c)" = i^4 thin Tr_H [{Tg^a, Tg^b} Tg^c] thin g^(mu nu) intk ((2k + p_3)^rho)/([k^2 - m^2][(k+p_3)^2 - m^2]) . $ <eq:g0_VVV-int-c>
Gauge bubble (e), channel $(b c|a)$ ($a$ at the cubic vertex): the line of colour $x$ carries $k$ from the quartic into the cubic vertex, the line of colour $y$ carries $k + p_1$ back. With @eq:V-vvvv for the legs $(b, nu), (c, rho), (x, alpha'), (y, beta')$, and the factor $1\/2$ for the two identical lines,
$ G_"(e)" & = i/2 thin fg^(a x y) intk (cal(V)^(mu alpha beta)(p_1, k, -k-p_1) thin g_(alpha alpha') g_(beta beta') thin Pi^(nu rho alpha' beta'))/(k^2 (k+p_1)^2), \
  Pi^(nu rho alpha' beta') & = fg^(b c e) fg^(x y e) (g^(nu alpha') g^(rho beta') - g^(nu beta') g^(rho alpha')) + fg^(b x e) fg^(c y e) (g^(nu rho) g^(alpha' beta') - g^(nu beta') g^(rho alpha')) \
  & quad + fg^(b y e) fg^(c x e) (g^(nu rho) g^(alpha' beta') - g^(nu alpha') g^(rho beta')) , $ <eq:g0_VVV-int-e>
with prefactor $1/2 (-i)^2 (-i) = i\/2$ ($-i$ from @eq:V-vvvv, $(-i)^2$ from the propagators).

*$Rop$-operation.* At one loop $Rbar G = G$. For the triangles $dg(G) = 4 + 3 - 6 = 1$, for the bubbles $dg(G) = 4 + 1 - 4 = 1$. For each integrand, multiply the numerator by its numerical prefactor ($+1$ for (a1), $-1$ for (b1), $1$ for (c), $i$ for (d), $i\/2$ for (e), $-i$ for (f1)) and keep the group factor apart. Write the result as $n_3 + n_2 + dots$, with $n_j$ of degree $j$ in $k$ (triangles; $n_1 + n_0$ for the bubbles), and let $w$ be the sum of the momentum shifts of the denominators. By @eq:uv-expansion,
$ Top(1) G = intk [ n_2/(k^2 - MUV^2)^ell - (2 k dot w) n_3/(k^2 - MUV^2)^(ell+1) ] , $ <eq:g0_VVV-T1>
with $ell$ the number of propagators (for the bubbles $n_3 -> n_1$, $n_2 -> n_0$). The term $n_3\/(k^2 - MUV^2)^ell$ is odd in $k$ and vanishes; $m^2$, $MUV^2$ and $p_i p_j$ have degree $2$ and are dropped, so $Kop Top(1) G$ depends neither on $m^2$ nor on $MUV$. All surviving integrals are logarithmic, and by @eq:tensor-red and @eq:logpole their residues at $ep = 0$ are
$ Kop intk (k^alpha k^beta)/(k^2 - MUV^2)^3 = ilf 1/ep thin (g^(alpha beta))/4, quad
  Kop intk (k^alpha k^beta k^gamma k^delta)/(k^2 - MUV^2)^4 = ilf 1/ep thin (g^(alpha beta) g^(gamma delta) + g^(alpha gamma) g^(beta delta) + g^(alpha delta) g^(beta gamma))/24 , $ <eq:g0_VVV-res>
with $g -> gmb$ for barred loop momenta (@eq:tensor-red-split). With $Kop Top(1) G = ilf thin N_0 \/ ep$, the contributions of the $n_2$ term and of the $n_3$ term are listed for each topology. No contraction $gmh^mu_(space mu)$ occurs, so $N_1 = 0$ for every diagram of this class (@sec:m-K).

_Fermion triangle._ Here $w = p_1 - p_2$, and $Pcb = (1 plus.minus g5)\/2$ with the upper sign for $f = Q, L$. The $n_2$ term with @eq:g0_VVV-res and $gbar^alpha gbar^nu gbar_alpha = -2 gbar^nu$, $gbar^alpha gbar^nu gbar^sigma gbar^rho gbar_alpha = -2 gbar^rho gbar^sigma gbar^nu$ (@eq:contr2) gives
$ 1/8 { -2 tr[gbar^mu gbar^nu gbar^rho slb(p)_1 (1 plus.minus g5)] + 2 tr[gbar^mu gbar^rho slb(p)_2 gbar^nu (1 plus.minus g5)] }, $
and the $n_3$ term, with the four-index rule and the same identities, gives
$ 1/12 tr[gbar^mu (gbar^nu gbar^rho slb(w) + gbar^rho slb(w) gbar^nu + slb(w) gbar^nu gbar^rho)(1 plus.minus g5)] . $
With the four-dimensional traces and @eq:tr5:
$ N_0^"(a1)" = Tr_f [Tg^a Tg^b Tg^c] { & underbrace(- gmb^(mu nu) (bar(p)_1 - bar(p)_2)^rho + gmb^(nu rho) bar(p)_3^mu - gmb^(rho mu) bar(p)_3^nu, n_2 "term") + underbrace(1/3 overline(Xi)^(mu nu rho)(bar(p)_1 - bar(p)_2), n_3 "term") \
  & plus.minus underbrace(i epsilon^(mu nu rho sigma) (bar(p)_1 - bar(p)_2)_sigma, n_2 "term") minus.plus underbrace(i epsilon^(mu nu rho sigma) (bar(p)_1 - bar(p)_2)_sigma, n_3 "term") } \
  = & - 2/3 thin Tr_f [Tg^a Tg^b Tg^c] thin overline(cal(V))^(mu nu rho) . $ <eq:g0_VVV-N-a>
The $g5$-even part is the tree tensor with barred metric and momenta. It contains no evanescent structure: the external indices sit on $gbar$, and a hatted external momentum could only enter through $hat(k) dot hat(p)_j$ in the denominators, which at degree one multiplies an odd number of hatted loop components and vanishes under @eq:tensor-red-split. The $g5$-odd part has no pole: the two terms cancel. This holds for each orientation and each multiplet separately, without the anomaly cancellation $D^(a b c) = 0$, and it is forced by symmetry: a single orientation is invariant under the cyclic relabelling $(mu, p_1) -> (nu, p_2) -> (rho, p_3)$ (a shift of $k$), and $epsilon^(mu nu rho sigma) v_sigma$ with $v$ linear and cyclically invariant has $v prop p_1 + p_2 + p_3 = 0$. Since $overline(cal(V))$ is antisymmetric under $(nu, p_2) <-> (rho, p_3)$, the two orientations combine into a commutator:
$ N_0^"(a1)" + N_0^"(a2)" = - 2/3 thin Tr_f [Tg^a [Tg^b, Tg^c]] thin overline(cal(V))^(mu nu rho) = - 2/3 thin i fg^(a b c) S_f^((a)) thin overline(cal(V))^(mu nu rho) . $ <eq:g0_VVV-N-a-sum>

_Higgs triangle._ $w = p_1 - p_3$. The $n_2$ term gives $-[g^(nu rho) p_1^mu + g^(rho mu) (p_1 - p_3)^nu - g^(mu nu) p_3^rho]$; with $n_3 = -8 k^mu k^nu k^rho$, the $n_3$ term gives $-8 k^mu k^nu k^rho (-2 k dot w) -> 2/3 thin Xi^(mu nu rho)(p_1 - p_3)$. Their sum is $1/3 cal(V)^(mu nu rho)$, and with the antisymmetry of $cal(V)$,
$ N_0^"(b1)" & = 1/3 Tr_H [Tg^a Tg^c Tg^b] cal(V)^(mu nu rho), \
  N_0^"(b1)" + N_0^"(b2)" & = - 1/3 Tr_H [Tg^a [Tg^b, Tg^c]] cal(V)^(mu nu rho) = - 1/3 thin i fg^(a b c) S_H^((a)) thin cal(V)^(mu nu rho) . $ <eq:g0_VVV-N-b>

_Seagull bubbles._ $w = p_3$, $n_1 = 2 g^(mu nu) k^rho$, $n_0 = g^(mu nu) p_3^rho$: the $n_0$ term gives $g^(mu nu) p_3^rho$, the $n_1$ term $2 g^(mu nu) k^rho (-2 k dot p_3) -> - g^(mu nu) p_3^rho$. Hence $N_0^"(c)" = 0$ in each channel; the integral @eq:g0_VVV-int-c vanishes identically by $k -> -k - p_3$.

_Gauge triangle._ $w = p_1 - p_3$, and $n_j$ is the degree-$j$ part of $i thin cal(N)_"(d)"$. At $p_j = 0$ every vertex of @eq:g0_VVV-int-d becomes
$ cal(V)^(lambda alpha beta)(0, k, -k) = - g^(lambda alpha) k^beta + 2 g^(alpha beta) k^lambda - g^(beta lambda) k^alpha , $ <eq:g0_VVV-V0>
and the contraction of the three, at $d = 4$ (only the residue is needed), is
$ n_3 = i [18 thin k^mu k^nu k^rho + 2 k^2 Xi^(mu nu rho)(k)] . $ <eq:g0_VVV-n3d>
With $k^mu k^nu k^rho k^sigma w_sigma -> (k^2)^2 Xi^(mu nu rho)(w) \/ 24$ and $k^2 Xi^(mu nu rho)(k) thin k dot w -> (k^2)^2 Xi^(mu nu rho)(w) \/ 4$ (@eq:g0_VVV-res),
$ -(2 k dot w) n_3 -> -2 i (18/24 + 2/4) Xi^(mu nu rho)(w) = - (5 i)/2 thin Xi^(mu nu rho)(p_1 - p_3) . $
The part $n_2$ is linear in the $p_j$. It comes from one vertex taken at its part linear in the external momenta, the other two at $p_j = 0$. These linear parts are, in the index names of @eq:g0_VVV-int-d,
$ & g^(mu alpha) p_1^beta + g^(alpha beta) p_1^mu - 2 g^(beta mu) p_1^alpha quad & "(vertex of" a")", \
  & g^(nu beta') (p_2 - p_1)^gamma + g^(beta' gamma) (p_1 - p_3)^nu + g^(gamma nu) (p_3 - p_2)^(beta') quad & "(vertex of" b")", \
  & 2 g^(rho gamma') p_3^(alpha') - g^(gamma' alpha') p_3^rho - g^(alpha' rho) p_3^(gamma') quad & "(vertex of" c")". $ <eq:g0_VVV-dV>
After the contraction and $k^alpha k^beta -> k^2 g^(alpha beta) \/ 4$ the three insertions give
$ & i [ -2 g^(mu nu) p_1^rho + 13/4 g^(nu rho) p_1^mu + 5/2 g^(rho mu) p_1^nu ] & "(vertex of" a")", \
  & i [ g^(mu nu) (1/2 p_1 + 5/2 p_2)^rho + g^(nu rho) (1/2 p_1 - 2 p_2)^mu + 13/4 g^(rho mu) (2 p_1 + p_2)^nu ] quad & "(vertex of" b")", \
  & i [ -13/4 g^(mu nu) p_3^rho + 2 g^(nu rho) p_3^mu - 5/2 g^(rho mu) p_3^nu ] & "(vertex of" c")", $ <eq:g0_VVV-n2d-parts>
and their sum, with $p_3 = -p_1 - p_2$, is
$ n_2 -> i [ g^(mu nu) (7/4 p_1 + 23/4 p_2)^rho + g^(nu rho) (7/4 p_1 - 4 p_2)^mu + 23/4 g^(rho mu) (p_1 - p_3)^nu ] . $ <eq:g0_VVV-n2d>
Adding the $n_3$ term, $-5/2 i thin Xi^(mu nu rho)(2 p_1 + p_2)$, gives $-(13 i)/4 cal(V)^(mu nu rho)$ (e.g. the coefficient of $g^(mu nu)$ is $(7/4 - 5) p_1 + (23/4 - 5/2) p_2 = -13/4 (p_1 - p_2)$), so with @eq:g0_VVV-group
$ N_0^"(d)" = - 13/8 thin i C_"ad" fg^(a b c) cal(V)^(mu nu rho) . $ <eq:g0_VVV-N-d>

_Gauge bubbles._ $w = p_1$. The cubic vertex is $cal(V)^(mu alpha beta)(p_1, k, -k-p_1) = g^(mu alpha)(p_1 - k)^beta + g^(alpha beta)(2k + p_1)^mu - g^(beta mu)(k + 2 p_1)^alpha$, the rest of the numerator is $k$-independent. Its $n_0$ term plus its $n_1$ term, $k^lambda (-2 k dot p_1) -> -p_1^lambda\/2$, give the effective vertex $3/2 (g^(mu alpha) p_1^beta - g^(mu beta) p_1^alpha)$, antisymmetric in $alpha beta$. Contracted with the three terms of $Pi$ it gives $2 Lambda^(nu rho)$, $Lambda^(nu rho)$ and $-Lambda^(nu rho)$ with $Lambda^(nu rho) = 3/2 (g^(mu nu) p_1^rho - g^(mu rho) p_1^nu)$, and with @eq:g0_VVV-group
$ N_0^"(e)"_((b c|a)) = i/2 (2 C_"ad" + 1/2 C_"ad" + 1/2 C_"ad") fg^(a b c) Lambda^(nu rho) = 9/4 thin i C_"ad" fg^(a b c) (g^(mu nu) p_1^rho - g^(mu rho) p_1^nu) . $
The three channels add up to the tree tensor:
$ sum_"channels" N_0^"(e)" = 9/4 thin i C_"ad" fg^(a b c) cal(V)^(mu nu rho) . $ <eq:g0_VVV-N-e>

_Ghost triangles._ $w = p_1 - p_3$. The $n_2$ term gives $-i/4 [g^(nu rho) p_1^mu - g^(rho mu) p_3^nu]$; with $n_3 = -i k^mu k^nu k^rho$ the $n_3$ term gives $-i k^mu k^nu k^rho (-2 k dot w) -> i/12 thin Xi^(mu nu rho)(p_1 - p_3)$. The group factor of (f2) is that of (f1) with $b <-> c$, i.e. $-C_"ad" fg^(a b c)\/2$, so the sum is $C_"ad" fg^(a b c)\/2$ times the antisymmetrisation in $(nu, p_2) <-> (rho, p_3)$, which is $i/12 thin cal(V)^(mu nu rho)$:
$ N_0^"(f1)" + N_0^"(f2)" = 1/24 thin i C_"ad" fg^(a b c) cal(V)^(mu nu rho) . $ <eq:g0_VVV-N-f>

*Result.* Collecting @eq:g0_VVV-N-a-sum to @eq:g0_VVV-N-f: the gauge and ghost loops give $(-13/8 + 9/4 + 1/24) C_"ad" = 2/3 thin C_"ad"$. Only $fg^(a b c)$ and $Tr_f [Tg^a [Tg^b, Tg^c]]$ occur, so every combination of gauge factors other than $G G G$ and $W W W$ has no pole: $fg^(a b c) = 0$ unless $a, b, c$ lie in one non-abelian factor; the hypercharge commutes with everything, so $Tr_f [Tg^a [Tg^b, Tg^c]] = 0$ as soon as one leg is $B$; and the traces with one colour or $"SU"(2)_L$ generator and the other two from other factors vanish because the generators are traceless. The symmetric traces $Tr_f [Tg^a {Tg^b, Tg^c}]$, non-zero per multiplet for $G G G$ ($d^(A B C)$), $G G B$, $W W B$ and $B B B$, multiply the $g5$-odd structure only, and that structure has no pole in either orientation (@eq:g0_VVV-N-a): the vanishing of these poles does not use $D^(a b c) = 0$. In the $g5$-even part and in the Higgs triangles the two orientations cancel for these traces (Furry).
#result[
  $ Gamma^((1))_(A^a_mu A^b_nu A^c_rho) |_"pole" = lf 1/ep N_0, \
    N_0 = [ (4 i)/3 thin fg^(a x y) fg^(b y z) fg^(c z x) - 1/3 Tr_H [Tg^a [Tg^b, Tg^c]] ] cal(V)^(mu nu rho) - 2/3 sum_f Tr_f [Tg^a [Tg^b, Tg^c]] thin overline(cal(V))^(mu nu rho) , $ <eq:g0_VVV-pole>
  i.e. for $a, b, c$ in the factor $r = G, W$, with $S^((r)) = S^((r))_"lh" + S^((r))_"rh"$ of @eq:SLR,
  $ N_0 = -i fg^(a b c) [ z^((r))_"bos" cal(V)^(mu nu rho) + z^((r))_"ferm" overline(cal(V))^(mu nu rho) ] = z^((r))_"bos" Gamma^((0))_(A A A) + z^((r))_"ferm" overline(Gamma)^((0))_(A A A), \
    z^((r))_"bos" = - 2/3 C^((r))_"ad" + 1/3 S^((r))_H , quad z^((r))_"ferm" = 2/3 S^((r)) , $ <eq:g0_VVV-pole-r>
  and $N_0 = 0$ for all other combinations. As in @sec:g0_VV, $z^((r))_"bos"$ is the coefficient of the bosonic loops, which multiply the $d$-dimensional tree vertex function $Gamma^((0))_(A A A) = -i fg^(a b c) cal(V)^(mu nu rho)$, and $z^((r))_"ferm"$ that of the fermion loops, which multiply its barred form $overline(Gamma)^((0))_(A A A) = -i fg^(a b c) overline(cal(V))^(mu nu rho)$. The counterterm is
  $ delta S^((1)) supset & lf 1/ep integral dif^d x thin fg^(a b c) [ z^((r))_"bos" (partial_mu A^a_nu) A^(b mu) A^(c nu) + z^((r))_"ferm" (bar(partial)_mu bar(A)^a_nu) bar(A)^(b mu) bar(A)^(c nu) ] \
    = & delta Z^((r))_(3 A) thin Q^((r))_(3 f) - lf 1/ep thin z^((r))_"ferm" thin tilde(E)^((r))_F |_(A A A) , $ <eq:g0_VVV-ct>
  $ 16 pi^2 ep thin delta Z^((r))_(3 A) = -z^((r))_"bos" - z^((r))_"ferm" = 2/3 C^((r))_"ad" - 1/3 S^((r))_H - 2/3 S^((r)) . $ <eq:g0_VVV-z1>
  Here $Q^((r))_(3 f)$ of @tab:Qbasis is continued with unsplit $d$-dimensional indices (@sec:evanescent), and $tilde(E)^((r))_F |_(A A A)$ is the part of $tilde(E)^((r))_F$ (@eq:evanescent) cubic in the gauge field: the monomial $Q^((r))_(3 f)$ built from $bar(A)$ and $bar(partial)$, minus $Q^((r))_(3 f)$. In terms of the constants of the table in @sec:st-solution, $delta Z^((r))_(3 A) = delta Z_(c_r) + 3/2 delta Z^((r))_A$ (a symbol of this subsection) is the coefficient of $Q^((r))_(3 f)$ in the assembled counterterm. The coefficient $-2/3 S^((r))$ of $tilde(E)^((r))_F$ is the one of @sec:g0_VV and @sec:g0_VVVV.
]

#figure(
  table(
    columns: 9,
    align: center,
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([factor], [(d)], [(e)], [(f)], [(b)], [(a) per $f$], [$-z^((r))_"bos"$], [$-z^((r))_"ferm"$], [$16 pi^2 ep thin delta Z^((r))_(3 A)$]),
    [$G G G$], [$-39/8 gs^2$], [$27/4 gs^2$], [$1/8 gs^2$], [$0$], [$Q: -2/3, thin U, D: -1/3$ $(gs^2)$], [$2 gs^2$], [$-4/3 gs^2$], [$2/3 gs^2$],
    [$W W W$], [$-13/4 gw^2$], [$9/2 gw^2$], [$1/12 gw^2$], [$-1/6 gw^2$], [$Q: -1, thin L: -1/3$ $(gw^2)$], [$7/6 gw^2$], [$-4/3 gw^2$], [$-1/6 gw^2$],
    [other], [$0$], [$0$], [$0$], [$0$], [$0$], [$0$], [$0$], [$0$],
  ),
  caption: [Standard-Model values. Each entry is the coefficient of $i fg^(a b c) cal(V)^(mu nu rho)$ in $N_0$, except column (a), which multiplies $i fg^(a b c) overline(cal(V))^(mu nu rho)$; $fg^(A B C) = gs f^(A B C)$, $fg^(I J K) = gw epsilon^(I J K)$. Inputs from @tab:invariants: $C^((G))_"ad" = 3 gs^2$, $C^((W))_"ad" = 2 gw^2$, $S^((W))_H = gw^2\/2$, $S^((G)) = 2 gs^2$, $S^((W)) = 2 gw^2$. Hence $N_0^(G G G) = i gs^3 f^(A B C) [2 cal(V) - 4/3 overline(cal(V))]$ and $N_0^(W W W) = i gw^3 epsilon^(I J K) [7/6 cal(V) - 4/3 overline(cal(V))]$.],
) <tab:g0_VVV-sm>

#check[*Symmetry.* Each of @eq:g0_VVV-N-a-sum to @eq:g0_VVV-N-f is a product of two totally antisymmetric objects, $fg^(a b c)$ and $cal(V)$ or $overline(cal(V))$, hence Bose symmetric; the gauge bubbles are so only after the sum over the three channels, and the single ghost, Higgs and fermion orientations only after the sum over both orientations. Each single orientation of (a), (b), (f) and the gauge triangle (d) is invariant under the cyclic relabelling of the legs, as it must be by a shift of $k$: their residues are proportional to the cyclically invariant $cal(V)$ (for (a), (b), (d)) before the sum over orientations.]

#check[*Auxiliary mass and masses.* At $dg(G) = 1$ no term with $MUV^2$ or $m^2$ survives in $Top(1) G$, so $N_0$ is independent of $MUV$ and of the Higgs mass parameter, and no $lambda$ or Yukawa coupling enters.]

#check[*Known limit.* In the four-dimensional part, $Pphys N_0 = i fg^(a b c) (16 pi^2 ep thin delta Z^((r))_(3 A)) overline(cal(V))$. For $"SU"(N)$ with $n_D$ Dirac fermions in the fundamental, Feynman gauge and an anticommuting $g5$ ($C_A$ the adjoint Casimir, $T_F$ and $T_S$ the indices of the fermion and scalar representations), $Z_1 = 1 + g^2\/(16 pi^2 ep) [2/3 C_A - 4/3 T_F n_D]$, and a complex scalar contributes $-1/3 T_S$. A chiral multiplet counts as half a Dirac fermion, so $-2/3 S^((r))$ is that fermion term, and the Higgs term $-1/3 S_H$ is the scalar one. The gauge-sector split $-13/8 + 9/4 + 1/24 = 2/3$ (triangle, bubbles, ghosts) reproduces the Feynman-gauge value.]

#check[*Relations to other classes.* With $delta Z^((r))_A$ of @sec:g0_VV, the coupling constant of @sec:st-solution is $16 pi^2 ep thin delta Z_(c_r) = 16 pi^2 ep (delta Z^((r))_(3 A) - 3/2 delta Z^((r))_A) = -11/6 C^((r))_"ad" + 1/6 S^((r))_H + 1/3 S^((r))$: $-29/6 gs^2$ for $"SU"(3)_c$ and $-35/12 gw^2$ for $"SU"(2)_L$ with one family and one Higgs doublet, so that $16 pi^2 beta_(c_r) = 2 c_r (16 pi^2 ep thin delta Z_(c_r))$ is the beta function of @sec:st-solution ($c_r$ the coupling of the factor $r$, @sec:bases). The matter loops enter $delta Z^((r))_(3 A)$ and $delta Z^((r))_A$ with the same coefficients ($-2/3 S^((r))$, $-1/3 S^((r))_H$), so $16 pi^2 ep (delta Z^((r))_(3 A) - delta Z^((r))_A) = -C^((r))_"ad"$ is matter independent. It equals $16 pi^2 ep (tilde(Z)^((r))_1 - tilde(Z)^((r))_3) = -C^((r))_"ad"\/2 - C^((r))_"ad"\/2$, the difference of the renormalisation constants of $tilde(srcA)^mu_a fg^(a b c) A^b_mu gh^c$ and $tilde(srcA)^mu_a partial_mu gh^a$ (@sec:g0_ghost_sources), as required by the Slavnov--Taylor identity for the four-dimensional parts. The complete list of the cross-class conditions is in @sec:st-solution.]

#check[*Evanescent part, with @sec:g0_VV and @sec:g0_VVVV.* The bosonic loops multiply the $d$-dimensional two-, three- and four-point monomials with unsplit indices. The fermion loops multiply the barred monomials, with the common coefficient $-2/3 S^((r))$ in $16 pi^2 ep thin delta S^((1))$. With the closed directions continued with unsplit indices (@sec:evanescent), the evanescent pole of the gauge sector is therefore $-lf ep^(-1) thin 2/3 S^((r)) thin tilde(E)^((r))_F$, with one coefficient for the two-, three- and four-point parts, plus the two-point operator $-lf ep^(-1) thin 1/3 S^((r)) thin E^((r))_A$ of @sec:g0_VV. This is the gauge part of the evanescent counterterm of @sec:st-solution.]

#general[At one loop the forest is empty and the fermion numerator is four-dimensional (@eq:phys-prop), so the fermion part of the pole is exactly the barred tree structure and the $g5$-odd part has no pole by the cyclic symmetry of one fermion loop. From two loops on, $-Kop Top(dg(G)) Rbar G$ contains the insertions of the one-loop counterterms of @eq:g0_VVV-ct in their $d$-dimensional form, including the evanescent term $-lf ep^(-1) thin 2/3 S^((r)) thin tilde(E)^((r))_F |_(A A A)$. Inside a loop their hatted components are contracted with loop momenta and give $eta\/ep^2$ terms, i.e. poles and finite local terms. Gauge lines then connect barred fermion vertices through the $d$-dimensional propagator, and $epsilon$ structures inside subdiagrams can multiply poles; the statement that the $g5$-odd part of $A A A$ has no pole must be re-derived at each order.]
