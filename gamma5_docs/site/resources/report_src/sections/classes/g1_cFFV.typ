#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Breaking of $omega overline(Psi) Psi A$ <sec:g1_cFFV>

External legs $gh^a (q)$, $Psi_f (p)$ incoming, $overline(Psi)_f (p')$ outgoing, $A^b_mu (p_A)$ incoming, $p' = p + q + p_A$, for every multiplet $f$ and every pair of gauge indices $a$, $b$ (same or different gauge factors; $r, s in {G, W, B}$ label gauge factors). The sum of the diagrams below is $i thin (Dbrk_1 dot Gamma)^((1))_(gh^a overline(Psi)_f Psi_f A^b_mu)$. Every graph has four propagators and a numerator of degree four, so its superficial degree of divergence is $dg(G) = 4 + 4 - 2 dot 4 = 0$; $N_0$ and $N_1$ are constants, and $Top(0)$ sets $p, p', q, p_A$ and $m^2$ to zero and replaces every denominator by $D = 1\/(k^2 - MUV^2)$.

*Diagrams.* A one-loop graph has exactly two internal lines at the $Dbrk_1$ vertex. Either both fermion legs of $Dbrk_1$ are internal and the ghost is external (graphs (a)--(d) of @fig:g1_cFFV-diagrams), or the ghost leg and one fermion leg are internal and the other fermion leg is external (graphs (e)--(g)). In (a)--(d) the fermion line runs from the external $Psi_f$ through $Dbrk_1$ to the external $overline(Psi)_f$; the loop closes through a gauge field or a Higgs field, and $A^b$ sits on one of the two internal fermion segments or on the exchanged line. In (e)--(g) the internal ghost $gh^c$ of $Dbrk_1$ ends on a ghost--gauge vertex @eq:V-ghost; the loop closes through one gauge line to the vertex @eq:V-ffv of the other external fermion; $A^b$ sits on the internal fermion line, on the gauge line, or on the ghost line. No Higgs line reaches the ghost, so (e)--(g) contain no Yukawa coupling. Each of (e)--(g) has two versions: $Dbrk_1$ at the outgoing leg $overline(Psi)_f$ (L, drawn) or at the incoming leg $Psi_f$ (R, all fermion arrows reversed). In (e)--(g) the internal fermion leg of $Dbrk_1$ is the sterile one; the term of @eq:V-delta with the sterile leg external carries the external momentum and vanishes under $Top(0)$.

#figure(
  stack(dir: ttb, spacing: 5mm,
  diagrow(
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("F", $overline(Psi)_f$, -1),)), (legs: (("F", $Psi_f$, 1),)), (legs: (("V", $A^b_mu$, 0),))),
      (("F", none, 1), ("V", $A^e$, 0), ("F", none, 1), ("F", none, 1)), scale: 0.85, caption: [(a)]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("F", $overline(Psi)_f$, -1),)), (legs: (("V", $A^b_mu$, 0),)), (legs: (("F", $Psi_f$, 1),))),
      (("F", none, 1), ("V", none, 0), ("V", none, 0), ("F", none, 1)), scale: 0.85, caption: [(b)]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("F", $overline(Psi)_f$, -1),)), (legs: (("F", $Psi_f$, 1),)), (legs: (("V", $A^b_mu$, 0),))),
      (("F", $Psi_(f')$, 1), ("S", $H$, 0), ("F", none, 1), ("F", none, 1)), scale: 0.85, caption: [(c)]),
    oneloop(
      ((legs: (("G", $gh^a$, 1),), mark: "delta"), (legs: (("F", $overline(Psi)_f$, -1),)), (legs: (("V", $A^b_mu$, 0),)), (legs: (("F", $Psi_f$, 1),))),
      (("F", $Psi_(f')$, 1), ("S", none, 0), ("S", none, 0), ("F", none, 1)), scale: 0.85, caption: [(d)]),
  ),
  diagrow(
    oneloop(
      ((legs: (("F", $overline(Psi)_f$, -1),), mark: "delta"), (legs: (("V", $A^b_mu$, 0),)), (legs: (("F", $Psi_f$, 1),)), (legs: (("G", $gh^a$, 1),))),
      (("F", $chi$, -1), ("F", none, -1), ("V", $A^e$, 0), ("G", $gh^c$, 1)), scale: 0.85, caption: [(e) L]),
    oneloop(
      ((legs: (("F", $overline(Psi)_f$, -1),), mark: "delta"), (legs: (("F", $Psi_f$, 1),)), (legs: (("V", $A^b_mu$, 0),)), (legs: (("G", $gh^a$, 1),))),
      (("F", $chi$, -1), ("V", $A^e$, 0), ("V", $A^(e')$, 0), ("G", $gh^c$, 1)), scale: 0.85, caption: [(f) L]),
    oneloop(
      ((legs: (("F", $overline(Psi)_f$, -1),), mark: "delta"), (legs: (("F", $Psi_f$, 1),)), (legs: (("G", $gh^a$, 1),)), (legs: (("V", $A^b_mu$, 0),))),
      (("F", $chi$, -1), ("V", $A^e$, 0), ("G", $gh^x$, 1), ("G", $gh^c$, 1)), scale: 0.85, caption: [(g) L]),
  )),
  caption: [Topologies of $(Dbrk_1 dot Gamma)^((1))_(gh overline(Psi) Psi A)$ (crossed circle: @eq:V-delta). (a), (b): gauge exchange; (c), (d): Higgs exchange with the Yukawa partner $f'$ of $f$ on the internal fermion line. In (a) and (c) $A^b$ may also sit on the segment between $Dbrk_1$ and $overline(Psi)_f$. (e)--(g): internal ghost; the line marked $chi$ is the sterile component next to $Dbrk_1$. The R versions of (e)--(g) have $Dbrk_1$ at the leg $Psi_f$.],
) <fig:g1_cFFV-diagrams>

#figure(
  table(
    columns: (auto, 1fr, auto, auto, 1.15fr),
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [sym.], [sign], [group factor (string from $overline(Psi)_f$ to $Psi_f$)]),
    [(a)], [$Psi_f$ (3 segments), $A^e$; $e$ in any factor carried by $f$; 2 positions of $A^b$], [$1$], [$+1$], [$Tg^e Tg^a Tg^b Tg^e$, $Tg^e Tg^b Tg^a Tg^e$],
    [(b)], [$Psi_f$ (2), $A^e$, $A^(e')$; $e, e'$ in the factor of $b$], [$1$], [$+1$], [$fg^(e b e') Tg^e Tg^a Tg^(e')$],
    [(c)], [$Psi_(f')$ (3), $H$; $f'$ the Yukawa partner of $f$; 2 positions], [$1$], [$+1$], [$y_(f f')^2 Tg^a_(f') Tg^b_(f')$ or $y_(f f')^2 Tg^b_(f') Tg^a_(f')$, contracted with the Yukawa tensors],
    [(d)], [$Psi_(f')$ (2), $H$ (2)], [$1$], [$+1$], [$y_(f f')^2 thin Tg^a_(f') Tg^b_H$, contracted],
    [(e) L, R], [$Psi_f$ (2), $A^e$, $gh^c$; $a, c, e in r$ non-abelian], [$1$], [$+1$], [L: $fg^(c e a) Tg^c Tg^b Tg^e$; R: $fg^(c e a) Tg^e Tg^b Tg^c$],
    [(f) L, R], [$Psi_f$, $A^e$, $A^(e')$, $gh^c$; $a, b, c, e, e' in r$], [$1$], [$+1$], [L: $fg^(e b e') fg^(c e' a) Tg^c Tg^e$; R: $fg^(e b e') fg^(c e' a) Tg^e Tg^c$],
    [(g) L, R], [$Psi_f$, $A^e$, $gh^x$, $gh^c$; $a, b, c, e, x in r$], [$1$], [$+1$], [L: $fg^(c b x) fg^(x e a) Tg^c Tg^e$; R: $fg^(c b x) fg^(x e a) Tg^e Tg^c$],
  ),
  caption: [Field assignments of @fig:g1_cFFV-diagrams. All generators act on the multiplet of the line they sit on. There is no closed loop; the sign $+1$ is the Wick sign of the open fermion and ghost lines relative to the leg order $gh, overline(Psi), Psi$ of @eq:V-delta (derived below @eq:g1_cFFV-int-e). $(f, f') in {(Q,U),(Q,D),(U,Q),(D,Q),(L,E),(E,L)}$ and $y_(f f')$ is $y_u$, $y_d$ or $y_e$. Per external field assignment the numbers of graphs agree with the automatic graph enumeration of @sec:automation (e.g. 17 for $gh^G overline(Q) Q G$: 7 of type (a), (b), 4 of type (c), 6 of type (e)--(g); 220 labelled multiplet-level graphs over the 27 assignments with a non-vanishing group factor). Ten further assignments ($gh^G overline(U) U W$, $gh^W overline(U) U G$, $gh^W overline(U) U B$, $gh^B overline(U) U W$, the same four for $D$, $gh^W overline(E) E B$, $gh^B overline(E) E W$) have graphs of type (a)--(d) whose group factors vanish identically. Graphs (e)--(g) exist only for $a$ in $"SU"(3)_c$ or $"SU"(2)_L$ and $f$ charged under that factor; (f), (g) need $b$ in the same factor.],
) <tab:g1_cFFV-diagrams>

*Integrands.* With all external momenta zero, every line of the loop carries the loop momentum $k$. In (a)--(d), $k$ flows along the fermion arrow through $Dbrk_1$, so $p = p' = k$ in @eq:V-delta. A graph of type (a)--(d) has the form
$ G_"(a)-(d)" = intk thin dots thin (i V_"out") thin (i sl(k))/(k^2) thin Tg^a (slh(k) Pcb - slh(k) Pc) thin (i sl(k))/(k^2) thin (i V_"in") thin dots , $ <eq:g1_cFFV-int-ad>
where $i V_"out"$ and $i V_"in"$ are the Feynman rules of the gauge (@eq:V-ffv) or Yukawa (@eq:V-yuk) vertices next to $Dbrk_1$ ($V_"out"$, $V_"in"$ their vertex functions), with $V_"out" = V_"out" Pc$ and $V_"in" = Pcb V_"in"$ ($Pc$ the projector of the line through $Dbrk_1$). In (e)--(g) L, $k$ flows along the fermion arrow into $Dbrk_1$ ($p = k$, $p' = 0$, vertex $Tg^c slh(k) Pcb$); the ghost line carries $-k$ along its arrow, so every ghost--gauge vertex @eq:V-ghost gives $+fg k$. In the R versions $k$ flows along the fermion arrow out of $Dbrk_1$ ($p = 0$, $p' = k$, vertex $-Tg^c slh(k) Pc$), the ghost carries $+k$ along its arrow, and every ghost--gauge vertex gives $-fg k$. The L graphs are
$ G_"(e,L)" &= intk Tg^c slh(k) Pcb (i sl(k))/(k^2) (i gbar^mu Pc Tg^b) (i sl(k))/(k^2) (i gbar^nu Pc Tg^e) (-i g_(nu rho))/(k^2) fg^(c e a) k^rho thin i/k^2 , $ <eq:g1_cFFV-int-e>
$ G_"(f,L)" &= intk Tg^c slh(k) Pcb (i sl(k))/(k^2) (i gbar^nu Pc Tg^e) (-i g_(nu nu'))/(k^2) fg^(e b e') Lambda^(nu' mu rho') (-k, 0, k) (-i g_(rho' rho))/(k^2) fg^(c e' a) k^rho thin i/k^2 , $ <eq:g1_cFFV-int-f>
$ G_"(g,L)" &= intk Tg^c slh(k) Pcb (i sl(k))/(k^2) (i gbar^nu Pc Tg^e) (-i g_(nu rho))/(k^2) fg^(x e a) k^rho thin i/k^2 thin fg^(c b x) k^mu thin i/k^2 , $ <eq:g1_cFFV-int-g>
with $Lambda^(nu' mu rho') (p_1, p_2, p_3) = g^(nu' mu) (p_1 - p_2)^(rho') + g^(mu rho') (p_2 - p_3)^(nu') + g^(rho' nu') (p_3 - p_1)^mu$ the bracket of @eq:V-vvv for the legs $(A^e_(nu'), A^b_mu, A^(e')_(rho'))$ with incoming momenta $(-k, 0, k)$. In the R versions the Dirac and group string is reversed, e.g. $G_"(e,R)" = intk (i gbar^nu Pc Tg^e) (i sl(k))/(k^2) (i gbar^mu Pc Tg^b) (i sl(k))/(k^2) (-Tg^c slh(k) Pc) (-i g_(nu rho))/(k^2) (-fg^(c e a) k^rho) i/k^2$.

_Wick sign._ In $(i Dbrk_1)(i S_"int")^3$ the fermion bilinears are contracted in the order $Psi thin overline(Psi)$ of neighbouring bilinears, which gives the string $overline(Psi) [dots] Psi$ without sign. The ghost $gh^c$ of $Dbrk_1$ stands to the left of the even string; moving it next to $partial agh^c$ of the ghost vertex and contracting in the order $gh thin agh$ gives no sign, and the remaining external $gh^a$ is moved back to the left of the even string, again without sign. The result is $gh^a overline(Psi) [dots] Psi thin A^b$ with the plain product of the rules, in the leg order of @eq:V-delta. The same steps applied to the graph (e) of @sec:g0_ghost_sources ($overline(srcF) gh Psi$, source vertex in place of $Dbrk_1$) reproduce its pole $fg^(a b c) Tg^b Tg^c Pc$, which fixes the conventions of the ghost and gauge propagators and of @eq:V-ghost used here.

*$Rop$-operation.* At one loop $Rbar G = G$, $dg(G) = 0$, and $Top(0) G$ is the integrand at zero external momenta with $1\/k^2 -> D$.

_Graphs (a)--(d) vanish._ With @eq:sandwich, $Pc sl(k) slh(k) Pcb sl(k) Pcb = hat(k)^2 slb(k) Pcb$ and $Pc sl(k) slh(k) Pc sl(k) Pcb = hat(k)^2 slb(k) Pcb$, hence
$ Pc thin (i sl(k))/(k^2) thin Tg^a (slh(k) Pcb - slh(k) Pc) thin (i sl(k))/(k^2) thin Pcb = - Tg^a (hat(k)^2 slb(k) - hat(k)^2 slb(k)) Pcb \/ (k^2)^2 = 0 . $ <eq:g1_cFFV-ad-zero>
The integrand of every graph (a)--(d) vanishes identically at $q = 0$ (for any $p$, $p_A$), so $Top(0) G = 0$: $N_0 = N_1 = 0$ for all graphs (a)--(d), gauge and Yukawa alike. This is @eq:m-tree-WI at $q = 0$: there the Feynman rule @eq:V-delta is $i V_(Dbrk_1) = Tg^a (slh(k) Pcb - slh(k) Pc) = Tg^a (Pcb thin sl(k) - sl(k) thin Pc)$, the difference of two rigid rotations of the inverse propagator. Each term cancels one adjacent propagator, and between physical neighbours the two remaining rotations are equal.

_Dirac algebra of (e)--(g)._ By @eq:m-sterile, $Tg^c slh(k) Pcb thin i sl(k) thin Pcb = i Tg^c hat(k)^2 Pcb$ (L) and $Pc thin i sl(k) (-slh(k) Pc) = -i hat(k)^2 Pc$ (R); the rest of the string is four-dimensional. The bracket of (f) contracted with the ghost momentum is $Lambda^(nu' mu rho') (-k,0,k) k_(rho') = -k^2 g^(nu' mu) + k^(nu') k^mu$. Dropping the constants, the numerators are
$ "(e)": thin hat(k)^2 bar(k)^2 thin gbar^mu Pc, quad
  "(f)": thin hat(k)^2 (-k^2 gbar^mu + slb(k) k^mu) Pc, quad
  "(g)": thin hat(k)^2 thin slb(k) k^mu thin Pc , $ <eq:g1_cFFV-num>
the same for L and R and for both chiralities ($g5$ enters only through $Pc$). The constants (products of the factors $i$, $-i$ and of the signs of the $Dbrk_1$ and ghost vertices) are $1, i, -i$ for (e, f, g) L and $1, -i, i$ for (e, f, g) R. All four propagators give $D^4$.

_Tensor reduction with the tag._ With @eq:tensor-red, $gmh^rho_(space rho) = -2 eta$ and $d(d+2) -> 24$,
$ hat(k)^2 bar(k)^2 -> (-2 eta dot 4)/24 (k^2)^2 = -eta/3 (k^2)^2, quad
  hat(k)^2 k^2 -> (-2 eta (d + 2))/(d(d+2)) (k^2)^2 = -eta/2 (k^2)^2, \
  hat(k)^2 slb(k) k^mu -> ((-2 eta) gbar^mu + 2 gmh^(alpha mu) gbar_alpha)/24 (k^2)^2 = -eta/12 gbar^mu (k^2)^2 , $ <eq:g1_cFFV-red>
and $Kop intk (k^2)^2 D^4 = ilf thin 1\/ep$. Every term is proportional to $eta$: $N_0 = 0$ graph by graph, also for an evanescent index $mu$ ($gmh^(alpha mu) gbar_alpha = 0$). The numerator of (f) gives $eta (1/2 - 1/12) = 5/12 eta$. Writing $Top(0) G = ilf (N_0 \/ ep + eta N_1 \/ ep) + dots$,
$ N_1^"(e,L)" = -1/3 fg^(c e a) Tg^c Tg^b Tg^e thin gbar^mu Pc, & quad N_1^"(e,R)" = -1/3 fg^(c e a) Tg^e Tg^b Tg^c thin gbar^mu Pc, \
  N_1^"(f,L)" = (5 i)/12 fg^(e b e') fg^(c e' a) Tg^c Tg^e thin gbar^mu Pc, & quad N_1^"(f,R)" = -(5 i)/12 fg^(e b e') fg^(c e' a) Tg^e Tg^c thin gbar^mu Pc, \
  N_1^"(g,L)" = i/12 fg^(c b x) fg^(x e a) Tg^c Tg^e thin gbar^mu Pc, & quad N_1^"(g,R)" = -i/12 fg^(c b x) fg^(x e a) Tg^e Tg^c thin gbar^mu Pc . $ <eq:g1_cFFV-N>

_Group algebra._ In (e) the sum L + R is $-1/3 fg^(c e a) (Tg^c Tg^b Tg^e + Tg^e Tg^b Tg^c) = 0$, by the antisymmetry of $fg^(c e a)$ in $c, e$. In (f) and (g) the sums contain $[Tg^c, Tg^e] = i fg^(c e g) Tg^g$ and a product of three structure constants, reduced with the Jacobi identity in the form $fg^(a i j) fg^(b j k) fg^(c k i) = 1/2 C_"ad" fg^(a b c)$ (@eq:invariants):
$ fg^(e b e') fg^(c e' a) fg^(c e g) = fg^(c b x) fg^(x e a) fg^(c e g) = -1/2 C_"ad" fg^(a b g) . $ <eq:g1_cFFV-fff>
Hence
$ N_1^"(e)" = 0, quad
  N_1^"(f)" = 5/24 C_"ad" fg^(a b c) Tg^c gbar^mu Pc, quad
  N_1^"(g)" = 1/24 C_"ad" fg^(a b c) Tg^c gbar^mu Pc, quad
  N_1^"(a)-(d)" = 0 . $ <eq:g1_cFFV-N-sum>
The symmetric products $fg fg {Tg, Tg}$, which would contain $d^(A B C)$ for quarks and $delta^(I J)$ for doublets, cancel between L and R in each topology. For $b$ in a factor other than that of $a$ only (e) exists, and it vanishes.

*Result.* $C_"ad" = C_"ad"^((r))$ is the adjoint invariant of the factor $r$ of $a$ (@tab:invariants); $Pphys N_1 = N_1$.
#result[
  $ (Dbrk_1 dot Gamma)^((1))_(gh^a overline(Psi)_f Psi_f A^b_mu) |_"pole" = lf 1/ep N_0, quad quad N_0 = 0 , $ <eq:g1_cFFV-pole>
  $ X_(1; gh^a overline(Psi)_f Psi_f A^b_mu) = lf thin Pphys N_1, quad
    Pphys N_1 = 1/4 C_"ad"^((r)) fg^(a b c) Tg^c_f thin gbar^mu P_(c_f) = -i/4 C_"ad"^((r)) [Tg^a_f, Tg^b_f] thin gbar^mu P_(c_f) , $ <eq:g1_cFFV-X>
  for $a, b in r$, $r in {G, W}$ and $f$ charged under $r$,
  and $Pphys N_1 = 0$ if $a$ is abelian, if $a$ and $b$ are in different factors, or if $f$ is a singlet of the factor of $a$. As a local functional,
  $ X_1 & supset lf thin i/4 sum_f sum_(r in {G, W}, thin r in f) C_"ad"^((r)) integral dif^4 x thin overline(psi)_f gamma^mu [wm^((r)), Am^((r))_mu] psi_f \
    & = lf sum_f sum_(r in {G,W}, thin r in f) C_"ad"^((r))/4 thin i (N^(f r r)_(F V) - N^(f r)_(V F)) . $ <eq:g1_cFFV-X-op>
]
Translation: in $N^(f r r)_(F V) = integral overline(psi)_f gamma^mu wm Am_mu psi_f = -integral gh^a A^b_mu thin overline(psi)_f gamma^mu Tg^a Tg^b psi_f$ the ghost is moved to the left of $overline(psi)_f$ (one sign). The coefficient of $gh^a A^b_mu overline(Psi) [dots] Psi$ in the leg order of @eq:V-delta is the vertex function, so $i (N^(f r r)_(F V) - N^(f r)_(V F))$ has the vertex function $-i [Tg^a, Tg^b] gbar^mu Pc = fg^(a b c) Tg^c gbar^mu Pc$, and @eq:g1_cFFV-X-op reproduces @eq:g1_cFFV-X. The structure is exactly the $gh overline(psi) psi A$ part of $bop Q^(f r)_R$ in @eq:btable-fermion; no $N^(f r s)_(F V)$ with $r != s$ and no $N_cal(J)$ structure arises.

#figure(
  table(
    columns: 7,
    align: (left, center, center, center, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$r$ of $gh^a$], [$Q$], [$U$], [$D$], [$L$], [$E$], [vertex function $16 pi^2 X_1$]),
    [$G$], [$3 gs^2 \/ 4$], [$3 gs^2 \/ 4$], [$3 gs^2 \/ 4$], [$0$], [$0$], [$3/4 gs^4 f^(A B C) T^C gbar^mu P_(c_f)$],
    [$W$], [$gw^2 \/ 2$], [$0$], [$0$], [$gw^2 \/ 2$], [$0$], [$1/4 gw^4 epsilon^(I J K) tau^K gbar^mu PL$],
    [$B$], [$0$], [$0$], [$0$], [$0$], [$0$], [$0$],
    [mixed $a$, $b$], [$0$], [$0$], [$0$], [$0$], [$0$], [$0$],
  ),
  caption: [Standard-Model values: coefficient of $i (N^(f r r)_(F V) - N^(f r)_(V F))$ in $16 pi^2 X_1$, i.e. $C_"ad"^((r))\/4$ with $C_"ad"^((G)) = 3 gs^2$, $C_"ad"^((W)) = 2 gw^2$, for each multiplet charged under $r$. Per topology: (f) gives $5/6$ and (g) $1/6$ of each entry; (a)--(e) give $0$. The pole @eq:g1_cFFV-pole vanishes for every entry.],
) <tab:g1_cFFV-sm>

#check[*Cancellation of the pole.* By @eq:m-b-rules at $q = 0$ the source-free counterterms of @sec:g0_FF and @sec:g0_FFV give no $gh overline(Psi) Psi A$ vertex under $bop_d$: the gauge-leg contraction is proportional to $q$, the rotations of the legs cancel by rigid invariance, and $overline(psi) Psi A A$ vertices do not exist. The source terms of @sec:g0_ghost_sources, $delta S^((1)) supset -lf 1/ep C_"ad"^((r))/2 (Q_rho [fg_((r))] + sum_f Q^(f r)_R + dots)$, contribute through $bop_d tilde(srcA)^mu_a supset overline(psi) gbar^mu Tg^a psi$ and $bop_d overline(srcF) supset -overline(psi) gbar^mu Am_mu$, $bop_d srcF supset -gbar^mu Am_mu psi$: $bop_d Q_rho [fg] |_(gh overline(psi) psi A) = integral fg^(a b c) (overline(psi) gbar^mu Tg^a psi) A^b_mu gh^c$ and $bop_d Q^(f r)_R |_(gh overline(psi) psi A) = i (N^(f r r)_(F V) - N^(f r)_(V F)) = - integral fg^(a b c) (overline(psi) gbar^mu Tg^a psi) A^b_mu gh^c$. These cancel, with no evanescent remainder (only barred vertices of $S_0$ enter). So $bop_d delta S^((1)) |_(gh overline(Psi) Psi A) = 0$, in agreement with $N_0 = 0$ in @eq:g1_cFFV-pole, as required by the check of @sec:m-oneloop.]

#check[*Position in the image of #bop.* $X_1$ of @eq:g1_cFFV-X-op lies in $bop("span"{Q_alpha})$: $bop Q^(f r)_R supset i (N^(f r r)_(F V) - N^(f r)_(V F))$ (@eq:btable-fermion). The only other monomial with a commutator $gh overline(psi) psi A$ part is $Q_rho [fg_((r))]$, through $bop tilde(srcA)$ (@eq:b-rho), and it is excluded by (N3). Hence $bop Delta S^((1)) = -X_1$ in this class requires $Delta S^((1)) supset -lf C_"ad"^((r))\/4 sum_(f: thin r in f) Q^(f r)_R$ for $r = G, W$. The same monomial then predicts, through the other terms of $bop Q^(f r)_R$, the contribution $lf C_"ad"^((r))\/4 (N^(f r r)_(R R) + overline(N)^(f r r)_(R R))$ to $X_1$ in @sec:g1_ccRF, and enters the $N^(f r)_F$ balance of @sec:g1_cFF and the Yukawa structures of @sec:g1_cFFS. The final solution is in @sec:st-solution. No mixed $N^(f r s)_(F V)$ ($r != s$) and no $N_cal(J)$ component arises: the mixed and symmetric parts cancel between the L and R graphs.]

#check[*Structure and auxiliary mass.* $MUV$ enters only through $Kop intk (k^2)^2 D^4 = ilf \/ ep$, which has no $MUV$ dependence; no $MUV^2$ term exists at degree $0$. $X_1$ is the same for left- and right-handed multiplets, involves no $epsilon$ tensor, and is universal: one coefficient $C_"ad"^((r))\/4$ for every multiplet charged under $r$, independent of the Yukawa couplings (graphs (c), (d) vanish). The Dirac strings @eq:g1_cFFV-num and the vanishing @eq:g1_cFFV-ad-zero were checked numerically with an explicit $8 times 8$ representation of the Clifford algebra with four barred and two hatted directions, for both chiralities; the group sums @eq:g1_cFFV-fff and @eq:g1_cFFV-N-sum were checked with explicit $"SU"(3)$ and $"SU"(2)$ generators. An independent implementation of all twelve integrands from the Feynman rules (symbolic BMHV algebra, and explicit Clifford matrices in $d = 4, dots, 12$ with $eta$ interpolated in $d - 4$; @sec:supporting) reproduces every $N_1$ of @eq:g1_cFFV-N, $N_0 = 0$, and the totals for all multiplets and all nine pairs of gauge factors.]

#check[*Literature.* For a simple group, @eq:g1_cFFV-X-op is eq. (6.25b) of #cite(<Belusca-Maito:2020ala>) at $xi = 1$, $(C_"ad"\/4) integral dif^4 x thin fg^(a b c) gh^a thin overline(psi) gamma^mu A^b_mu Tg^c psi$, including the sign; there too only the four graphs (f), (g) contribute, the two graphs (e) cancel and the graphs with scalar exchange vanish.]

#general[The vanishing of the graphs with the external ghost on $Dbrk_1$ uses only $q = 0$, which the degree-$0$ local part imposes, and the physical neighbours of $Dbrk_1$: by @eq:m-tree-WI the vertex is then the difference of two rigid rotations at the ends of a physical-to-physical propagator, which cancel inside the integrand. It holds in any loop graph whose $Dbrk_1$ vertex carries an external ghost of zero momentum and has two physical neighbours, including graphs with insertions of $delta S^((1))$ and $Delta S^((1))$, whose fermion vertices are physical (@sec:evanescent). At two loops the local part of this class therefore comes from graphs with an internal $Dbrk_1$ ghost and from the insertions (b)--(d) of @sec:m-beyond, among them $bop_d Delta S^((1))$ with the $Q^(f r)_R$ term required here.]
