#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Ghost and BRST-source functions <sec:g0_ghost_sources>

External legs: $srcA^mu_a (-q) thin gh^b (q)$; $srcA^mu_a (-p-q) thin A^b_nu (p) thin gh^c (q)$; $srcC_a (-q_1 - q_2) thin gh^b (q_1) thin gh^c (q_2)$; $overline(srcF)_f thin gh^a (q) thin Psi_f (p)$ and $gh^a (q) thin overline(Psi)_f (p_2) thin srcF_f$ for every multiplet $f$; $srcH^dagger gh^a (q) H(p)$ and $gh^a (q) H^dagger (p_2) srcH$. All momenta are incoming; $p_2$ is the incoming momentum of $overline(Psi)_f$ or $H^dagger$ (the outgoing momentum $p'$ of @sec:feynman-rules is $p' = -p_2$), and the momentum of the source leg is fixed by conservation. The sum of the diagrams below is $i thin Gamma^((1))$ for these legs. The superficial degree of divergence is $dg(G) = 1$ for $srcA gh$ and $dg(G) = 0$ for all other functions. By the ghost equation @eq:rhotilde the functions with $agh^a$ in place of $srcA^mu_a$ follow by $srcA -> tilde(srcA)$ and are not computed separately. The gauge factor $r$ is the one of the external ghost; all three factors $G$, $W$, $B$ are covered.

*Diagrams.* The ghost couples only through the ghost--gauge vertex @eq:V-ghost and the source vertices @eq:V-src-gauge and @eq:V-src-matter. Each source appears once and is external. A 1PI graph therefore contains the source vertex with two internal lines, one of which is a ghost line that ends on a ghost--gauge vertex; the other line closes the loop through one gauge exchange. Graphs in which the external ghost sits on the source vertex are one-particle reducible. There is no closed fermion or ghost loop, no Yukawa or quartic vertex can enter, and the topologies are those of @fig:g0_ghost_sources-diagrams.

#figure(
  stack(dir: ttb, spacing: 5mm,
  diagrow(
    oneloop(
      ((legs: (("K", $srcA^mu_a$, 0),), mark: "src"), (legs: (("G", $gh^b$, 1),))),
      (("G", none, -1), ("V", none, 0)), caption: [(a) $srcA gh$]),
    oneloop(
      ((legs: (("K", $srcA^mu_a$, 0),), mark: "src"), (legs: (("V", $A^b_nu$, 0),)), (legs: (("G", $gh^c$, 1),))),
      (("G", none, -1), ("G", none, -1), ("V", none, 0)), caption: [(b) $srcA A gh$, ghost line]),
    oneloop(
      ((legs: (("K", $srcA^mu_a$, 0),), mark: "src"), (legs: (("G", $gh^c$, 1),)), (legs: (("V", $A^b_nu$, 0),))),
      (("G", none, -1), ("V", none, 0), ("V", none, 0)), caption: [(c) $srcA A gh$, gauge line]),
  ),
  diagrow(
    oneloop(
      ((legs: (("K", $srcC_a$, 0),), mark: "src"), (legs: (("G", $gh^b$, 1),)), (legs: (("G", $gh^c$, 1),))),
      (("G", none, -1), ("V", none, 0), ("G", none, 1)), caption: [(d) $srcC gh gh$]),
    oneloop(
      ((legs: (("K", $overline(srcF)_f$, 0),), mark: "src"), (legs: (("G", $gh^a$, 1),)), (legs: (("F", $Psi_f$, 1),))),
      (("G", none, -1), ("V", none, 0), ("F", none, 1)), caption: [(e) $overline(srcF) gh Psi$]),
    oneloop(
      ((legs: (("K", $srcH^dagger$, 0),), mark: "src"), (legs: (("G", $gh^a$, 1),)), (legs: (("S", $H$, 1),))),
      (("G", none, -1), ("V", none, 0), ("S", none, 1)), caption: [(f) $srcH^dagger gh H$]),
  )),
  caption: [One-loop topologies of the class (square: source vertex). The conjugate functions $gh overline(Psi)_f srcF_f$ (e$'$) and $gh H^dagger srcH$ (f$'$) have the drawings (e), (f) with all matter arrows reversed.],
) <fig:g0_ghost_sources-diagrams>

#figure(
  table(
    columns: 5,
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [symmetry factor], [Wick sign], [group factor]),
    [(a)], [$gh$, $A$ (factor $r$)], [$1$], [$+1$], [$fg^(a e x) fg^(x e b) = - C_"ad"^((r)) delta^(a b)$],
    [(b)], [$gh$, $gh$, $A$], [$1$], [$+1$], [$fg^(a e x) fg^(x b y) fg^(y e c) = -1/2 C_"ad"^((r)) fg^(a b c)$],
    [(c)], [$gh$, $A$, $A$], [$1$], [$+1$], [$fg^(a e x) fg^(b e h) fg^(x h c) = 1/2 C_"ad"^((r)) fg^(a b c)$],
    [(d)], [$gh$, $gh$, $A$], [$1$], [$+1$], [$fg^(a x e) fg^(x h b) fg^(e h c) = 1/2 C_"ad"^((r)) fg^(a b c)$],
    [(e)], [$gh$, $A$, $Psi_f$; $f$ with $r in f$], [$1$], [$+1$], [$fg^(b e a) Tg^b_f Tg^e_f = i/2 C_"ad"^((r)) Tg^a_f$],
    [(e$'$)], [$gh$, $A$, $Psi_f$; $f$ with $r in f$], [$1$], [$+1$], [$fg^(b e a) Tg^e_f Tg^b_f = -i/2 C_"ad"^((r)) Tg^a_f$],
    [(f), (f$'$)], [$gh$, $A$, $H$; $r = W$], [$1$], [$+1$], [as (e), (e$'$) with $Tg_H$],
  ),
  caption: [Field assignments. Repeated adjoint indices $e, h, x, y$ are summed. The internal gauge line belongs to the factor $r$ of the external ghost, since $fg^(a b c)$ vanishes between different factors. The group factors follow from @eq:fg and @eq:invariants; the commutator gives $fg^(b e a) Tg^b Tg^e = 1/2 fg^(b e a) [Tg^b, Tg^e] = i/2 C_"ad"^((r)) Tg^a$. For the abelian ghost every group factor vanishes.],
) <tab:g0_ghost_sources-diagrams>

The Grassmann signs follow from Wick's theorem applied to $i^3 chevron.l cal(L)_1 cal(L)_2 cal(L)_3 chevron.r$ (graph (a): $i^2 chevron.l cal(L)_1 cal(L)_2 chevron.r$), where $cal(L)_1, cal(L)_2, cal(L)_3$ are the interaction terms of the three vertices, written as in @eq:Lgf and @eq:Lext-explicit and the external fields kept in the order of the tree-level term. In (a)--(c) and (f) the contracted pairs $chevron.l gh thin agh chevron.r$ already stand next to each other. In (d) the odd factors are $gh^x (z_1) gh^e (z_1) | agh (z_2) gh^b (z_2) | agh (z_3) gh^c (z_3)$, with $z_1, z_2, z_3$ the vertex positions; two transpositions give $chevron.l gh^x agh(z_2) chevron.r chevron.l gh^e agh(z_3) chevron.r gh^b gh^c$, sign $+1$. The $1\/2$ of $-1/2 srcC fg gh gh$ is cancelled by the two contractions of the ghosts at the $srcC$ vertex, and the $1\/2!$ of the two identical vertices $agh A gh^"ext"$ by their exchange; graph (d) with labelled legs has the rule $-i fg^(a x e)$ of @eq:V-src-gauge and symmetry factor $1$. In (e) the odd factors $gh^b (z_1) Psi(z_1) | agh(z_2) gh^a (z_2) | overline(Psi)(z_3) Psi^"ext" (z_3)$ need two transpositions to reach $chevron.l gh^b agh chevron.r chevron.l Psi(z_1) overline(Psi)(z_3) chevron.r gh^a Psi^"ext"$, sign $+1$; in (e$'$) $overline(Psi)(z_1)$ moves past four odd factors to the right of $Psi(z_3)$, sign $+1$; in (f$'$) two transpositions, sign $+1$. The value of each graph is therefore the product of the Feynman rules of @sec:feynman-rules, read along the open ghost line.

*Integrands.* The ghost line carries $k$ (or $k_1$, $k_2$) from the ghost--gauge vertex to the source vertex, so that the outgoing antighost momentum of @eq:V-ghost is $k$. In (a) the gauge line carries $q - k$; in (b) the ghost momentum between the $gh^c$ and $A^b$ vertices is $k - p$ and the gauge line carries $q + p - k$; in (c) the gauge lines carry $q - k$ (from the $gh^c$ vertex to the three-gauge vertex) and $q + p - k$; in (d)--(f) the gauge line carries $l$ from the external-matter (or $gh^c$) vertex to the $gh$ vertex, so that the ghost carries $q + l$ and the matter line $p - l$ (in (e$'$), (f$'$): $l - p_2$, from the source vertex to the gauge vertex). With the rules of @sec:feynman-rules (@eq:propagators to @eq:V-src-matter):
$ G_"(a)" = intk (i fg^(a e x) g_mu^(space alpha)) (i)/(k^2) (- fg^(x e b) k^beta) (-i g_(alpha beta))/((k-q)^2)
  = i C_"ad"^((r)) delta^(a b) intk k_mu/(k^2 (k-q)^2) , $ <eq:g0_ghost_sources-int-a>
$ G_"(b)" = - fg^(a e x) fg^(x b y) fg^(y e c) intk (k_nu (k - p)_mu)/(k^2 (k-p)^2 (k - p - q)^2)
  = 1/2 C_"ad"^((r)) fg^(a b c) intk (k_nu (k - p)_mu)/(k^2 (k-p)^2 (k - p - q)^2) , $ <eq:g0_ghost_sources-int-b>
$ G_"(c)" = - 1/2 C_"ad"^((r)) fg^(a b c) intk (cal(V)_(nu mu beta) k^beta)/(k^2 (k-q)^2 (k-q-p)^2), \
  cal(V)_(nu mu beta) = g_(nu mu) (2p + q - k)_beta + g_(mu beta) (2k - 2q - p)_nu + g_(beta nu) (q - k - p)_mu , $ <eq:g0_ghost_sources-int-c>
where $cal(V)_(nu mu beta)$ (a symbol of this subsection) is the bracket of @eq:V-vvv with the legs $(A^b_nu (p), A_mu (k - q - p), A_beta (q - k))$ and the numerical prefactors are $i dot i dot (-1) dot (-i) = -1$ in (a), $i dot i dot (-1) dot i dot (-1) dot (-i) = -1$ in (b) and $i dot i dot (-1) dot (-i) dot (-i) = -1$ in (c);
$ G_"(d)" = integral_l (-i fg^(a x e)) (i)/(k_1^2) (-fg^(x h b) k_1^nu) (i)/(k_2^2) (-fg^(e h c) k_2^lambda) (-i g_(nu lambda))/(l^2) \
  = 1/2 C_"ad"^((r)) fg^(a b c) integral_l (k_1 dot k_2)/(k_1^2 k_2^2 l^2), quad k_1 = q_1 + l, thin k_2 = q_2 - l , $ <eq:g0_ghost_sources-int-d>
with $l$ the loop momentum and $integral_l$ the measure of $intk$ with $k -> l$;
$ G_"(e)" = integral_l (-Tg^b_f P_(c_f)) (i)/((q+l)^2) (-fg^(b e a) (q+l)^nu) (i (sl(p) - sl(l)))/((p - l)^2) (i gbar^mu P_(c_f) Tg^e_f) (-i g_(mu nu))/(l^2) \
  = -i/2 C_"ad"^((r)) Tg^a_f integral_l (P_(c_f) (sl(p) - sl(l)) gbar^nu P_(c_f) (q + l)_nu)/((q+l)^2 (p-l)^2 l^2) , $ <eq:g0_ghost_sources-int-e>
$ G_"(e"'")" = i/2 C_"ad"^((r)) Tg^a_f integral_l (gbar^nu P_(c_f) (sl(l) - sl(p)_2) P_(overline(c)_f) (q + l)_nu)/((q+l)^2 (l - p_2)^2 l^2) , $ <eq:g0_ghost_sources-int-e2>
$ G_"(f)" & = -i/2 C_"ad"^((r)) Tg^a_H integral_l ((q + l) dot (2 p - l))/((q+l)^2 ((p - l)^2 - m^2) l^2), \
  G_"(f"'")" & = -i/2 C_"ad"^((r)) Tg^a_H integral_l ((q + l) dot (l - 2 p_2))/((q+l)^2 ((l - p_2)^2 - m^2) l^2) . $ <eq:g0_ghost_sources-int-f>
The prefactors of the rules are $(-1) dot i dot (-1) dot i dot i dot (-i) = -1$ in (e), (e$'$), (f) and $i dot i dot 1 dot i dot (-1) dot (-i) = +1$ in (f$'$), times the group factors of @tab:g0_ghost_sources-diagrams. The vertex $H^dagger H A$ of @eq:V-ssv gives $i (2p - l)^mu Tg^e_H$ in (f) and $i (l - 2 p_2)^mu Tg^e_H$ in (f$'$).

*$Rop$-operation.* At one loop $Rbar G = G$, and $Top(dg(G))$ acts on the integrands above.

_(a), $dg(G) = 1$._ By @eq:uv-expansion, $1\/k^2 = 1\/(k^2 - MUV^2) + O(MUV^2)$ and $1\/(k-q)^2 = 1\/(k^2 - MUV^2) + 2 k dot q\/(k^2 - MUV^2)^2 + O(2)$. The degree-0 term $k_mu \/ (k^2 - MUV^2)^2$ is odd; $MUV^2$ has degree two. Hence
$ Top(1) G_"(a)" = i C_"ad"^((r)) delta^(a b) intk (2 k_mu thin k dot q)/((k^2 - MUV^2)^3) = i C_"ad"^((r)) delta^(a b) 2/d thin q_mu intk (k^2)/((k^2 - MUV^2)^3) \
  = ilf 1/ep (i/2 C_"ad"^((r)) delta^(a b) q_mu) + O(ep^0) . $ <eq:g0_ghost_sources-T-a>

_(b)--(f), $dg(G) = 0$._ $Top(0)$ replaces every denominator by $1\/(k^2 - MUV^2)$ (or $1\/(l^2 - MUV^2)$) and keeps the part of the numerator quadratic in the loop momentum; masses, $MUV$ and external momenta do not enter. With @eq:tensor-red and @eq:logpole the numerators reduce as follows.
- (b): $k_nu k_mu -> g_(mu nu) k^2 \/ d$.
- (c): $cal(V)_(nu mu beta) k^beta -> -g_(mu nu) k^2 + k_mu k_nu -> (-1 + 1\/d) g_(mu nu) k^2 = -3/4 g_(mu nu) k^2 + O(ep)$.
- (d): $k_1 dot k_2 -> -l^2$.
- (e): $P_(c_f) sl(l) gbar^nu P_(c_f) l_nu = P_(c_f) sl(l) P_(overline(c)_f) slb(l) = slb(l) slb(l) P_(c_f) = bar(l)^2 P_(c_f)$ by @eq:phys-prop, so the numerator is $-bar(l)^2 P_(c_f) -> -(4\/d) l^2 P_(c_f)$ by @eq:tensor-red-split.
- (e$'$): $gbar^nu P_(c_f) sl(l) P_(overline(c)_f) l_nu = slb(l) slb(l) P_(overline(c)_f) -> (4\/d) l^2 P_(overline(c)_f)$.
- (f): $(q+l) dot (2p - l) -> -l^2$; (f$'$): $(q+l) dot (l - 2 p_2) -> + l^2$.
The Dirac algebra of (e), (e$'$) contains only $gbar$ and $bar(l)$: no hatted object and no trace $gmh^mu_(space mu)$ occurs, and no hatted external momentum survives at degree zero. The factors $1\/d$ and $4\/d$ come from unsplit or barred pairings and are $1\/4$, $1$ at $ep = 0$. With $intk (k^2)\/(k^2 - MUV^2)^3 -> ilf \/ ep$, $Top(dg(G)) G = ilf thin N_0 \/ ep + dots$ with
$ N_0^"(a)" = i/2 C_"ad"^((r)) delta^(a b) q_mu, quad
  N_0^"(b)" = 1/8 C_"ad"^((r)) fg^(a b c) g_mu^(space nu), quad
  N_0^"(c)" = 3/8 C_"ad"^((r)) fg^(a b c) g_mu^(space nu), quad
  N_0^"(d)" = -1/2 C_"ad"^((r)) fg^(a b c), $ <eq:g0_ghost_sources-N0>
$ N_0^"(e)" = i/2 C_"ad"^((r)) Tg^a_f P_(c_f), quad
  N_0^"(e"'")" = i/2 C_"ad"^((r)) P_(overline(c)_f) Tg^a_f, quad
  N_0^"(f)" = i/2 C_"ad"^((r)) Tg^a_H, quad
  N_0^"(f"'")" = -i/2 C_"ad"^((r)) Tg^a_H . $ <eq:g0_ghost_sources-N0b>
$N_1 = 0$ for every graph: there is no $Dbrk_1$ vertex and no closed evanescent contraction. The pole parts do not depend on $MUV$ or $m^2$.

*Result.* The tree-level vertex functions in the leg order written are the rules of @eq:V-src-gauge and @eq:V-src-matter divided by $i$: $Gamma^((0))_(srcA^mu_a gh^b (q)) = -i q_mu delta^(a b)$, $Gamma^((0))_(srcA^mu_a A^b_nu gh^c) = fg^(a b c) g_mu^(space nu)$, $Gamma^((0))_(srcC_a gh^b gh^c) = -fg^(a b c)$, $Gamma^((0))_(overline(srcF)_f gh^a Psi_f) = i Tg^a_f P_(c_f)$, $Gamma^((0))_(gh^a overline(Psi)_f srcF_f) = i P_(overline(c)_f) Tg^a_f$, $Gamma^((0))_(srcH^dagger gh^a H) = i Tg^a_H$, $Gamma^((0))_(gh^a H^dagger srcH) = -i Tg^a_H$. Reordering two odd legs changes the sign of both the tree and the one-loop function, so the ratios below are independent of the order. $C_"ad"^((r))$ is the adjoint invariant of the factor of the external ghost, and $fg^(a b c) Tg^b_f Tg^c_f = i/2 C_"ad"^((r)) Tg^a_f$.
#result[
  $ Gamma^((1))_(srcA^mu_a gh^b (q)) |_"pole" & = lf 1/ep thin i/2 C_"ad"^((r)) delta^(a b) q_mu = - 1/2 C_"ad"^((r)) lf 1/ep thin Gamma^((0))_(srcA^mu_a gh^b), \
    Gamma^((1))_(srcA^mu_a A^b_nu gh^c) |_"pole" & = lf 1/ep thin 1/2 C_"ad"^((r)) fg^(a b c) g_mu^(space nu), quad
    Gamma^((1))_(srcC_a gh^b gh^c) |_"pole" = lf 1/ep (-1/2 C_"ad"^((r)) fg^(a b c)), \
    Gamma^((1))_(overline(srcF)_f gh^a Psi_f) |_"pole" & = lf 1/ep thin fg^(a b c) Tg^b_f Tg^c_f P_(c_f), quad
    Gamma^((1))_(gh^a overline(Psi)_f srcF_f) |_"pole" = lf 1/ep thin P_(overline(c)_f) fg^(a b c) Tg^b_f Tg^c_f, \
    Gamma^((1))_(srcH^dagger gh^a H) |_"pole" & = lf 1/ep thin fg^(a b c) Tg^b_H Tg^c_H, quad
    Gamma^((1))_(gh^a H^dagger srcH) |_"pole" = - lf 1/ep thin fg^(a b c) Tg^b_H Tg^c_H . $ <eq:g0_ghost_sources-pole>
  All functions except $srcA gh$ are $+1/2 C_"ad"^((r)) lf \/ ep$ times their tree. With $tilde(srcA)$ of @eq:rhotilde and the factor-$r$ restrictions $wm^((r))$, $Tg^a_((r))$ of @sec:bases,
  $ delta S^((1)) supset - lf 1/ep sum_(r = G, W, B) C_"ad"^((r))/2 integral dif^d x thin [ - tilde(srcA)^mu_a partial_mu gh^a_((r)) + fg^(a b c)_((r)) tilde(srcA)^mu_a A^b_mu gh^c - 1/2 fg^(a b c)_((r)) srcC_a gh^b gh^c \
    + sum_(f: thin r in f) ( i thin overline(srcF)_f wm^((r)) P_(c_f) Psi_f + i thin gh^a_((r)) overline(Psi)_f P_(overline(c)_f) Tg^a_((r)) srcF_f )
    + delta_(r W) ( i thin srcH^dagger wm^((W)) H - i thin gh^a_((W)) H^dagger Tg^a_((W)) srcH ) ] , $ <eq:g0_ghost_sources-ct>
  i.e. in the basis of @tab:Qbasis
  $ delta S^((1)) supset - lf 1/ep sum_(r) C_"ad"^((r))/2 [ - Q^((r))_(rho partial) + Q_rho [fg_((r))] + Q_zeta [-1/2 fg_((r))] + sum_(f: thin r in f) Q^(f r)_R + delta_(r W) Q^((W))_sigma ] \
    = - lf 1/ep sum_r C_"ad"^((r))/2 [ cal(O)^((r))_c - 2 Q^((r))_(rho partial) ] . $ <eq:g0_ghost_sources-ct-basis>
  In @eq:g0_ghost_sources-ct-basis the $Q_alpha$ stand for their $d$-dimensional continuations defined in @sec:evanescent (unsplit indices in the bosonic monomials, the physical fermion components in the fermion source terms), with every index and derivative $d$-dimensional as written in @eq:g0_ghost_sources-ct: the poles of @eq:g0_ghost_sources-pole contain $q_mu = bar(q)_mu + hat(q)_mu$ and $g_mu^(space nu) = gmb_mu^(space nu) + gmh_mu^(space nu)$ with equal weights, so the hatted parts $integral hat(tilde(srcA))^mu_a hat(partial)_mu gh^a$ and $integral fg^(a b c) hat(tilde(srcA))^mu_a hat(A)^b_mu gh^c$ are included in $Q^((r))_(rho partial)$ and $Q_rho [fg_((r))]$ and are not separate counterterms. The fermion and Higgs source terms contain no hatted object. All other ghost-sector monomials ($Q_rho [gs d]$, the seven Kronecker $Q_rho [t]$, $Q_zeta [(srcC^A gh^A gh^B)]$, $Q_zeta [(srcC^I gh^I gh^B)]$, $Q^(f B)_R$, $Q^((B))_sigma$) have coefficient zero. There is no evanescent pole in this class.
]

#figure(
  table(
    columns: 4,
    align: (left, center, center, center),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([operator in @eq:g0_ghost_sources-ct-basis], [$r = G$], [$r = W$], [$r = B$]),
    [$Q^((r))_(rho partial)$], [$-3/2 gs^2$], [$-gw^2$], [$0$],
    [$Q_rho [fg_((r))]$, $Q_zeta [-1/2 fg_((r))]$], [$3/2 gs^2$], [$gw^2$], [--],
    [$Q^(f r)_R$], [$3/2 gs^2$ ($f = Q, U, D$)], [$gw^2$ ($f = Q, L$)], [$0$ (all $f$)],
    [$Q^((r))_sigma$], [--], [$gw^2$], [$0$],
    [$tilde(Z)^((r))_3 - 1$], [$3 gs^2 \/ (32 pi^2 ep)$], [$gw^2 \/ (16 pi^2 ep)$], [$0$],
    [$tilde(Z)^((r))_1 - 1$], [$-3 gs^2 \/ (32 pi^2 ep)$], [$-gw^2 \/ (16 pi^2 ep)$], [$0$],
  ),
  caption: [Standard Model values. Rows 1--4: coefficient $c$ in $delta S^((1)) supset -lf (1\/ep) thin c thin Q$, from $C_"ad"^((G)) = 3 gs^2$, $C_"ad"^((W)) = 2 gw^2$, $C_"ad"^((B)) = 0$. Rows 5--6: the constants of @eq:g0_ghost_sources-Z. "--" means that the monomial does not exist.],
) <tab:g0_ghost_sources-sm>

Adding @eq:g0_ghost_sources-ct to @eq:Lgf and @eq:Lext gives, per factor $r$,
$ tilde(Z)^((r))_3 tilde(srcA)^mu_a partial_mu gh^a_((r)) + tilde(Z)^((r))_1 [ fg^(a b c)_((r)) tilde(srcA)^mu_a A^b_mu gh^c - 1/2 fg^(a b c)_((r)) srcC_a gh^b gh^c + (sop Psi, sop H "terms of factor" r) ], \
  tilde(Z)^((r))_3 = 1 + C_"ad"^((r))/2 lf 1/ep, quad tilde(Z)^((r))_1 = 1 - C_"ad"^((r))/2 lf 1/ep . $ <eq:g0_ghost_sources-Z>
$tilde(Z)_3$ multiplies the ghost kinetic term $partial^mu agh^a partial_mu gh^a$ (through $tilde(srcA)$): it is the ghost field renormalisation. These are the constants $tilde(Z)^((r))_3 - 1 = 1/2 (delta Z^((r))_gh - delta Z^((r))_A)$ and $tilde(Z)^((r))_1 - 1 = delta Z_(c_r) + 1/2 delta Z^((r))_gh$ of the table of constants in @sec:st-solution. In the assembled counterterm the coefficient of $Q^((r))_(rho partial)$ is $tilde(Z)^((r))_3 - 1$ and that of every other source monomial is $tilde(Z)^((r))_1 - 1$. In terms of @eq:b-action the source terms of $delta S^((1))$ change the BRST transformations to
$ sop A^a_mu -> tilde(Z)_3 thin partial_mu gh^a + tilde(Z)_1 fg^(a b c) A^b_mu gh^c, quad
  sop gh^a -> tilde(Z)_1 thin sop gh^a, \
  sop Psi_f -> i (tilde(Z)_1 wm) P_(c_f) Psi_f, quad
  sop H -> i (tilde(Z)_1 wm) H , $
with $tilde(Z)_1 wm = sum_r tilde(Z)^((r))_1 wm^((r))$: the linear term of $sop A$ is renormalised by $tilde(Z)_3$, all non-linear terms by the common factor $tilde(Z)_1$. For the abelian ghost, $C_"ad"^((B)) = 0$ and no 1PI graph exists. Every three-leg vertex of $gh^B$, $agh^B$, $srcA_B$ or $srcC_B$ in @eq:V-ghost and @eq:V-src-gauge carries $fg$ and vanishes, and the two-point vertex $srcA_B gh^B$ is tree level only; the only remaining vertices of $gh^B$ are the matter source vertices of @eq:V-src-matter, where an external $gh^B$ next to the external source leaves a single internal line, so the graph cannot be 1PI. $tilde(Z)^((B))_3 = tilde(Z)^((B))_1 = 1$: $gh^B$ stays free, the hypercharge rotations in $sop Psi_f$, $sop H$ are not renormalised at one loop, and $delta S^((1))$ contains no $gh^B$ from this class, in agreement with the abelian ghost equation and with (N4) of @sec:normalisation.

#check[*Slavnov--Taylor structure.* In four dimensions the ghost-sector monomials enter the closed directions of @eq:closed only through $cal(O)^((r))_A$ ($-Q^((r))_(rho partial)$), $cal(O)^((r))_c$ and $cal(O)^((r))_g$ (both with $Q_rho [fg]$, $Q_zeta [-fg\/2]$, $Q^(f r)_R$, $Q^((W))_sigma$ with coefficient $1$). An invariant one-loop pole $1/2 delta Z^((r))_A cal(O)^((r))_A + 1/2 delta Z^((r))_gh cal(O)^((r))_c + delta Z_(c_r) cal(O)^((r))_g$ (as in the assembled counterterm of @sec:st-solution) therefore needs equal coefficients of $Q_rho [fg_((r))]$, $Q_zeta [-fg_((r)) \/ 2]$, every $Q^(f r)_R$ and $Q^((W))_sigma$. The independent graphs (b)+(c), (d), (e), (e$'$), (f), (f$'$) give exactly this, $delta Z_(c_r) + 1/2 delta Z^((r))_gh = -lf C_"ad"^((r)) \/ (2 ep)$, with the Grassmann signs derived above; a wrong sign in any Wick reordering would break the equality. No ghost-number-1 graph of $Dbrk_1 dot Gamma$ has these external legs at one loop (@sec:normalisation), so the four-dimensional poles of this class must be $bop$-closed by themselves, as they are.]

#check[*Relations to other classes.* With $tilde(Z)^((r))_3 - 1 = lf C_"ad"^((r)) \/ (2 ep)$ and $tilde(Z)^((r))_1 - 1 = -lf C_"ad"^((r)) \/ (2 ep)$ from this class, the closed directions predict for the four-dimensional invariant parts of the poles (coefficients $[Q]$ in $16 pi^2 ep thin delta S^((1))$, $r = G, W$): $[Q^((r))_(3 f)] - [Q^((r))_T] = 16 pi^2 ep (tilde(Z)^((r))_1 - tilde(Z)^((r))_3) = -C_"ad"^((r))$ (@sec:g0_VV, @sec:g0_VVV), i.e. $delta Z^((r))_(3 A) - delta Z^((r))_A = -lf C_"ad"^((r)) \/ ep$ for the three-gauge-vertex and gauge-field constants; $[Q^(f r)_V] - [Q^f_K] = -C_"ad"^((r))$ (@sec:g0_FF, @sec:g0_FFV); $[Q^((W))_J] + [Q_(S S)] = C_"ad"^((W))$ (@sec:g0_SS, @sec:g0_SSV_SSVV). For $r = B$ the ghost sector vanishes, which forces $[Q^(f B)_V] = [Q^f_K]$ and $[Q^((B))_J] + [Q_(S S)] = 0$ (abelian Ward identities), while $[Q^((B))_T]$ is unconstrained. All these conditions hold; they are listed with the assembled counterterm in @sec:st-solution.]

#check[*Ghost equation.* Replacing the external $srcA^mu_a$ by the antighost replaces the vertex $srcA A gh$ by $agh A gh$ with the factor $-fg^(a b c) p_0^mu$, with $p_0$ the outgoing momentum of the external antighost; the integrands of (a)--(c) are unchanged otherwise, so $Gamma_(agh gh)$ and $Gamma_(agh A gh)$ are $i p_(0 mu)$ times $Gamma_(srcA^mu gh)$ and $Gamma_(srcA^mu A gh)$, as @eq:rhotilde requires. *Hermiticity.* (e$'$) and (f$'$) carry the same coefficient as (e) and (f), so the counterterm is a multiple of the hermitian $Q^(f r)_R$, $Q^((r))_sigma$. *Mass independence.* $m^2$ and $MUV$ do not enter any pole.]

#general[At one loop the ghost and source functions contain no fermion loop and no $Dbrk_1$ vertex, so their poles are $d$-dimensional covariant, and their four-dimensional parts are $bop$-closed (in $d$ dimensions $bop_d cal(O)^((r))_c = Dbrk^((r))_1$, @eq:bd-Oc). From two loops on the forest of @eq:m-Rbar subtracts subgraphs with the counterterms of @eq:g0_ghost_sources-ct and of the gauge self-energy, including its evanescent parts $E^((r))_A$ of @eq:evanescent inserted on the internal gauge line; fermion loops enter the gauge line, and source-dependent counterterms make the source--field products $cal(S)_"q" (C)$ of @eq:m-Delta-l non-zero. The ghost-sector poles are then no longer closed by themselves and the source classes of $X_l$ need not vanish.]
