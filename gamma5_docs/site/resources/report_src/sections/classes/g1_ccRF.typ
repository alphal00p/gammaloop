#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Breaking with BRST sources: $omega omega overline(R) Psi$ <sec:g1_ccRF>

External legs $gh^a (q_a)$, $gh^c (q_c)$, $overline(srcF)_f (p_R)$, $Psi_f (p)$, and the conjugate set $gh^a (q_a)$, $gh^c (q_c)$, $overline(Psi)_f (p_2)$, $srcF_f (p_R)$; all momenta are incoming with zero sum, $a, c$ are any two of the twelve ghost indices and $f in {Q, U, D, L, E}$. $p_2$ is the incoming momentum of $overline(Psi)_f$; the outgoing momentum $p'$ of @sec:feynman-rules is $p' = -p_2$. The sum of the diagrams below is the corresponding part of $i thin (Dbrk_1 dot Gamma)^((1))$. The legs have total dimension $5$, so $dg(G) = 0$ and $Top(0)$ evaluates the integrand at zero external momenta. The same subsection shows that the other ghost-number-1 functions with a BRST source, $gh gh tilde(srcA)$, $gh gh tilde(srcA) A$, $gh gh gh srcC$, $gh gh srcH^dagger H$ and $gh gh H^dagger srcH$, have no one-loop graph with a $Dbrk_1$ vertex.

*Diagrams.* A source never propagates, so $overline(srcF)_f$ sits on the vertex $overline(srcF)_f gh^b Psi_f$ of @eq:V-src-matter. In a one-loop 1PI graph every vertex lies on the loop, so $gh^b$ and $Psi_f$ of the source vertex are internal. The internal ghost ends on an antighost, i.e. on a vertex @eq:V-ghost. An external ghost cannot sit on the source vertex (that vertex would have a single internal line). It therefore sits on $Dbrk_1$ or on a vertex @eq:V-ghost. This gives two topologies:
- (a) The $Dbrk_1$ ghost is external ($gh^a$). The vertex @eq:V-ghost next to the source carries $gh^c$, and its gauge line ends on the gauge vertex @eq:V-ffv of the external $Psi_f$. Along the fermion line the order is $Psi_f ->$ @eq:V-ffv $-> Dbrk_1 ->$ source. The other order, $Psi_f -> Dbrk_1 ->$ @eq:V-ffv $->$ source, is one-particle reducible: $Dbrk_1$ is joined to the loop by a single fermion propagator.
- (b) The $Dbrk_1$ ghost is internal. It ends on a second vertex @eq:V-ghost, and the two vertices @eq:V-ghost carry the external ghosts: $gh^c$ next to the source and $gh^a$ next to $Dbrk_1$. They are joined by their gauge line. The external $Psi_f$ enters at $Dbrk_1$, and the loop is source $<-$ ghost $<-$ @eq:V-ghost $-$ gauge $-$ @eq:V-ghost $->$ ghost $-> Dbrk_1 ->$ fermion $->$ source, with $(I, V) = (4, 4)$.

Every other assignment has $L = I - V + 1 >= 2$. Yukawa vertices on the fermion line come in pairs joined by an internal Higgs line, and in (b) the gauge lines could instead both end on the fermion line; in both cases $(I, V) = (7, 6)$. The conjugate legs give (a$'$) and (b$'$). In (b) the ghost chain fixes one non-abelian factor $r$ for $gh^a$, $gh^c$ and all internal ghost and gauge lines. In (a) only $gh^c$ must be non-abelian. The external ghosts enter as fields of $i (Dbrk_1 dot Gamma)$. The exchange $a <-> c$ is contained in the antisymmetry of $gh^a gh^c$ and is not a separate graph. In (b) the two placements of the identical vertices @eq:V-ghost cancel the $1\/2!$. In (a) the two terms of @eq:V-delta put the sterile component on the line that enters $Dbrk_1$ or on the line that leaves it. They are the rows (a1), (a2) of @tab:g1_ccRF-diagrams, which are parts of one graph.

#figure(
  grid(columns: 2, column-gutter: 14mm, row-gutter: 6mm,
    oneloop(
      ((legs: (("K", $overline(srcF)_f$, 0),), mark: "src"), (legs: (("G", $gh^c$, 1),)), (legs: (("F", [#h(3mm)$Psi_f$], 1),)), (legs: (("G", $gh^a$, 1),), mark: "delta")),
      (("G", $gh^b$, -1), ("V", $A^d$, 0), ("F", $k$, 1), ("F", [$k + q_a$#h(5mm)], 1)), r: 0.75, leg: 0.8, caption: [(a)]),
    oneloop(
      ((legs: (("K", $overline(srcF)_f$, 0),), mark: "src"), (legs: (("G", $gh^c$, 1),)), (legs: (("G", $gh^a$, 1),)), (legs: (("F", $Psi_f$, 1),), mark: "delta")),
      (("G", $ell$, -1), ("V", $A^d$, 0), ("G", $ell'$, 1), ("F", $k$, 1)), r: 0.75, leg: 0.8, caption: [(b)]),
    oneloop(
      ((legs: (("K", $srcF_f$, 0),), mark: "src"), (legs: (("G", $gh^c$, 1),)), (legs: (("F", [#h(3mm)$overline(Psi)_f$], -1),)), (legs: (("G", $gh^a$, 1),), mark: "delta")),
      (("G", $gh^b$, -1), ("V", $A^d$, 0), ("F", $k$, -1), ("F", [$k - q_a$#h(5mm)], -1)), r: 0.75, leg: 0.8, caption: [(a$'$)]),
    oneloop(
      ((legs: (("K", $srcF_f$, 0),), mark: "src"), (legs: (("G", $gh^c$, 1),)), (legs: (("G", $gh^a$, 1),)), (legs: (("F", $overline(Psi)_f$, -1),), mark: "delta")),
      (("G", $ell$, -1), ("V", $A^d$, 0), ("G", $ell'$, 1), ("F", $k$, -1)), r: 0.75, leg: 0.8, caption: [(b$'$)]),
  ),
  caption: [The one-loop graphs of $(Dbrk_1 dot Gamma)^((1))_(gh gh overline(srcF) Psi)$, (a) and (b), and of the conjugate function, (a$'$) and (b$'$). Square: source vertex @eq:V-src-matter; crossed circle: @eq:V-delta; the vertices with a single external ghost and no mark are @eq:V-ghost, the remaining one in (a), (a$'$) is @eq:V-ffv.],
) <fig:g1_ccRF-diagrams>

#figure(
  {
  set par(justify: false)
  table(
    columns: (auto, 1fr, auto, auto, 1fr),
    align: (center, left, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([diagram], [internal lines], [sym.], [sign], [group factor]),
    [(a1)], [$gh^b$, $A^d$ ($b, d, c$ in one non-abelian factor), $Psi_f$; sterile part on the line $k$ (entering $Dbrk_1$)], [$1$], [$+1$], [$Upsilon^"(a)"_f = fg^(b d c) gh^a gh^c Tg^b Tg^a Tg^d$],
    [(a2)], [same; sterile part on the line $k + q_a$ (leaving $Dbrk_1$)], [$1$], [$+1$], [$Upsilon^"(a)"_f$],
    [(a$'$1), (a$'$2)], [same, fermion arrows reversed], [$1$], [$+1$], [$overline(Upsilon)^"(a)"_f = fg^(b d c) gh^a gh^c Tg^d Tg^a Tg^b$],
    [(b)], [$gh^b$, $A^d$, $gh^e$ (all indices, $a$ and $c$ included, in one non-abelian factor $r$), $Psi_f$ (sterile part)], [$1$], [$+1$], [$Upsilon^"(b)"_f = fg^(b d c) fg^(e d a) gh^a gh^c Tg^b Tg^e$],
    [(b$'$)], [same, fermion arrow reversed], [$1$], [$+1$], [$overline(Upsilon)^"(b)"_f = fg^(b d c) fg^(e d a) gh^a gh^c Tg^e Tg^b$],
  )},
  caption: [Field assignments of @fig:g1_ccRF-diagrams. There is no closed loop; the sign is the Wick sign derived below, for the ghosts in the order $gh^a gh^c$ ($gh^a$ on the $Dbrk_1$ side). The generators act on the multiplet $f$. In (a) $gh^c$ belongs to $"SU"(3)_c$ or $"SU"(2)_L$ and $gh^a$ to any factor carried by $f$, $gh^B$ included; in (b) both belong to the same non-abelian factor. The multiplet $E$ has no graph.],
) <tab:g1_ccRF-diagrams>

_Other source functions._ The legs $gh gh tilde(srcA)$ (one derivative), $gh gh tilde(srcA) A$, $gh gh gh srcC$, $gh gh srcH^dagger H$, $gh gh H^dagger srcH$ contain no fermion and no fermion source. Without $overline(srcF)$, $srcF$, every vertex of $S_0$ that contains $Psi_f$ contains exactly one $Psi_f$ and one $overline(Psi)_(f')$ (@eq:V-ffv, @eq:V-yuk, @eq:V-delta), so fermion lines close into loops. The two fermion legs of $Dbrk_1$ are internal, and the loop of a one-loop graph is the closed fermion loop through $Dbrk_1$: every vertex on it carries two fermion lines of the loop. The source vertices @eq:V-src-gauge, $srcH^dagger gh H$ and $gh H^dagger srcH$ contain no fermion and are not on the loop. Joined to it, they are attached through a bridge, and the graph is not 1PI. Hence
$ (Dbrk_1 dot Gamma)^((1))_(gh gh tilde(srcA)) = (Dbrk_1 dot Gamma)^((1))_(gh gh tilde(srcA) A) = (Dbrk_1 dot Gamma)^((1))_(gh gh gh srcC) = (Dbrk_1 dot Gamma)^((1))_(gh gh srcH^dagger H) = (Dbrk_1 dot Gamma)^((1))_(gh gh H^dagger srcH) = 0 . $ <eq:g1_ccRF-others>
All five first appear at two loops:
- $gh gh tilde(srcA)$: the vertex $srcA A gh$ with internal $A$ and $gh$ (the vertex $srcA gh$ has a single internal line, always a bridge). The ghost goes to @eq:V-ghost with the second external ghost, and both gauge lines go to a fermion loop with $Dbrk_1$: $(I, V) = (6, 5)$.
- $gh gh tilde(srcA) A$: the same with the external gauge field on the fermion loop, $(7, 6)$.
- $gh gh gh srcC$: both ghosts of $srcC gh gh$ are internal, each goes to a vertex @eq:V-ghost with an external ghost, and both gauge lines go to the fermion loop with $Dbrk_1$: $(7, 6)$.
- $gh gh srcH^dagger H$ and the conjugate: $gh$ and $H$ of the source vertex are internal. The ghost goes through @eq:V-ghost and a gauge line to the fermion loop, and the Higgs line and the external $H$ sit on two Yukawa vertices of the loop: $(7, 6)$.
Each has $L = 2$.

*Integrands.* Each graph is a term $gh^a gh^c thin overline(srcF)_f (p_R) [G] Psi_f (p)$ of $i (Dbrk_1 dot Gamma)^((1))$, with the external fields as Grassmann-odd background fields. _Wick signs._ The vertices of @eq:Delta1, @eq:Lext-explicit, @eq:Lgf and @eq:Lferm are written in the order $Dbrk_1$, source, ghost vertices, gauge vertex ($overline(srcF)$ is even). In (a) the odd factors are
$ gh^a overline(Psi)(x_1) Psi(x_1) | gh^b (x_2) Psi(x_2) | agh^e (x_3) gh^c | overline(Psi)(x_4) Psi^"ext" . $
The permutation to $gh^a gh^c chevron.l gh^b agh^e chevron.r chevron.l Psi(x_2) overline(Psi)(x_1) chevron.r chevron.l Psi(x_1) overline(Psi)(x_4) chevron.r Psi^"ext"$ has $12$ inversions. In (b) they are
$ gh^e overline(Psi)(x_1) Psi^"ext" | gh^b (x_2) Psi(x_2) | agh^(b') (x_3) gh^c | agh^(e') (x_4) gh^a , $
and the permutation to $gh^a gh^c chevron.l gh^b agh^(b') chevron.r chevron.l gh^e agh^(e') chevron.r chevron.l Psi(x_2) overline(Psi)(x_1) chevron.r Psi^"ext"$ has $26$ inversions. For (a$'$), (b$'$) the analogous reorderings, ending in $overline(Psi)^"ext" chevron.l Psi overline(Psi) chevron.r dots$, are even as well. All four signs are $+1$, and each value is the product of the Feynman rules of @sec:feynman-rules read along the fermion line from the source vertex.

_(a)._ The routing is $k_1 = k$ on the fermion line from @eq:V-ffv to $Dbrk_1$, $k_2 = k + q_a$ from $Dbrk_1$ to the source vertex, $ell = -(k + q_a + p_R)$ on the ghost line along its arrow and $k - p$ on the gauge line. In @eq:V-delta $p -> k_1$, $p' -> k_2$, and in @eq:V-ghost $p' -> ell$:
$ G_"(a)" = intk (-Tg^b Pc) (i sl(k)_2)/(k_2^2) thin Tg^a (slh(k)_1 Pcb - slh(k)_2 Pc) thin (i sl(k)_1)/(k_1^2) (i gbar^mu Tg^d Pc) (-i g_(mu nu))/((k - p)^2) (-fg^(b d c) ell^nu) i/(ell^2) , $ <eq:g1_ccRF-int>
with $Pc = P_(c_f)$. The sterile propagator sits between $Dbrk_1$ and the physical projectors of its neighbours. By @eq:m-sterile and @eq:phys-prop, for any two momenta,
$ Pc (i sl(k)_2)/(k_2^2) (slh(k)_1 Pcb - slh(k)_2 Pc) (i sl(k)_1)/(k_1^2) Pcb
  = (hat(k_2)^2 thin slb(k)_1 - hat(k_1)^2 thin slb(k)_2)/(k_1^2 k_2^2) thin Pcb . $ <eq:g1_ccRF-sandwich>
The term $-hat(k_1)^2 slb(k)_2$ is (a1) and $hat(k_2)^2 slb(k)_1$ is (a2). With $Pcb gbar^mu Pc = gbar^mu Pc$,
$ G_"(a)" = i thin Upsilon^"(a)"_f intk ((hat(k_2)^2 slb(k)_1 - hat(k_1)^2 slb(k)_2) slb(ell) thin Pc)/(k_1^2 k_2^2 (k - p)^2 ell^2) . $ <eq:g1_ccRF-G>
For (a$'$), with $k_1 = k$ from $Dbrk_1$ to @eq:V-ffv, $k_2 = k - q_a$ from the source vertex to $Dbrk_1$ and $ell' = k - q_a - p_R$, the same steps give
$ G_"(a')" = i thin overline(Upsilon)^"(a)"_f intk (slb(ell)' (hat(k_1)^2 slb(k)_2 - hat(k_2)^2 slb(k)_1) Pcb)/(k_1^2 k_2^2 (k + p_2)^2 ell'^2) . $ <eq:g1_ccRF-Gbar>

_(b)._ The routing is $k$ on the fermion line from $Dbrk_1$ to the source vertex, $ell = -(k + p_R)$ on the ghost line into the source vertex, $ell' = k - p$ on the ghost line into $Dbrk_1$, and $k + p_R + q_c$ on the gauge line. In @eq:V-delta $p' -> k$, and the antighost momenta of @eq:V-ghost are $ell$ and $ell'$:
$ G_"(b)" = intk (-Tg^b Pc) (i sl(k))/(k^2) thin Tg^e (slh(p) Pcb - slh(k) Pc) (-fg^(b d c) ell^mu) (-fg^(e d a) ell'^nu) (-i g_(mu nu))/((k + p_R + q_c)^2) i/(ell^2) i/(ell'^2) . $
The constants multiply to $1$, and $Pc sl(k) slh(p) Pcb = slb(k) slh(p) Pcb$, $Pc sl(k) slh(k) Pc = hat(k)^2 Pc$, so
$ G_"(b)" = Upsilon^"(b)"_f intk ((slb(k) slh(p) Pcb - hat(k)^2 Pc) thin ell dot ell')/(k^2 (k + p_R + q_c)^2 ell^2 ell'^2) . $ <eq:g1_ccRF-Gb>
For (b$'$), with $k$ on the fermion line from the source vertex to $Dbrk_1$, $ell = k - p_R$, $ell' = -(k + p_2)$ and the outgoing fermion momentum of @eq:V-delta equal to $p' = -p_2$,
$ G_"(b')" = overline(Upsilon)^"(b)"_f intk ((hat(k)^2 + slh(p)_2 slb(k)) Pcb thin ell dot ell')/(k^2 (k - p_R - q_c)^2 ell^2 ell'^2) . $ <eq:g1_ccRF-Gbp>
In @eq:g1_ccRF-Gb and @eq:g1_ccRF-Gbp, $ell dot ell'$ is the unsplit $d$-dimensional product of the gauge propagator.

_Group factors._ From $[Tg^d, wm] = i fg^(d c b) gh^c Tg^b$ (@eq:fg) and $gh^a gh^c = -gh^c gh^a$, summed over $d$:
$ Upsilon^"(a)"_f = i [Tg^d, wm] wm Tg^d, quad overline(Upsilon)^"(a)"_f = -i Tg^d wm [Tg^d, wm], quad
  Upsilon^"(b)"_f = [Tg^d, wm] [Tg^d, wm] = - overline(Upsilon)^"(b)"_f . $ <eq:g1_ccRF-Theta1>
Only $d$ in a non-abelian factor $r$ contributes, and then $[Tg^d, wm] = [Tg^d, wm^((r))]$. Write $wm = wm^((r)) + wm'$ and define, for this subsection only, $C^((r))_f$ by $sum_(d in r) Tg^d_f Tg^d_f = C^((r))_f bb(1)$, so that $sum_r C^((r))_f = C_f$ of @eq:invariants. Then $sum_(d in r) Tg^d X Tg^d = (C^((r))_f - C^((r))_"ad"\/2) X$ for $X$ in the Lie algebra of $r$ (in particular $X = wm^((r))$ and $X = wm^((r)) wm^((r)) = i/2 thin fg^(x y z) gh^x gh^y Tg^z$), and $C^((r))_f X$ for $X$ commuting with the factor $r$. Hence
$ Upsilon^"(a)"_f = -i/2 sum_(r = G, W) C^((r))_"ad" thin wm^((r)) (wm - wm^((r))) = - overline(Upsilon)^"(a)"_f, quad quad
  Upsilon^"(b)"_f = - 1/2 sum_(r = G, W) C^((r))_"ad" thin wm^((r)) wm^((r)) = - overline(Upsilon)^"(b)"_f . $ <eq:g1_ccRF-Theta>
$Upsilon^"(a)"$ contains only ghosts of two different factors, $Upsilon^"(b)"$ only ghosts of one non-abelian factor.

*$Rop$-operation.* At one loop $Rbar G = G$. $Top(0)$ sets $q_a = q_c = p = p_2 = p_R = 0$, every denominator becomes $D_MUV = 1\/(k^2 - MUV^2)$, and the numerators of @eq:uv-expansion, as well as the terms with $slh(p)$, $slh(p)_2$, have degree $>= 1$ and do not contribute. By @eq:tensor-red with split metrics (cross pairings vanish, $gmh dot gmb = 0$),
$ hat(k)^2 bar(k)^2 -> -(8 eta)/(d (d+2)) (k^2)^2, quad quad hat(k)^2 k^2 -> ((k^2)^2)/(d (d+2)) (gmh^alpha_(space alpha) g^beta_(space beta) + 2 gmh^alpha_(space alpha)) = -(2 eta)/d (k^2)^2 , $ <eq:g1_ccRF-tens>
and $intk (k^2)^2 D_MUV^4 -> ilf 1\/ep$ (@eq:logpole).

_(a)._ $k_1 = k_2 = k$, $ell = -k$. With $slb(ell) = -slb(k)$ the two terms of the numerator of @eq:g1_ccRF-G are $-hat(k)^2 slb(k) (-slb(k)) = hat(k)^2 bar(k)^2$ (a1) and $hat(k)^2 slb(k) (-slb(k)) = -hat(k)^2 bar(k)^2$ (a2), so
$ N_0^"(a1)" = N_0^"(a2)" = 0, quad N_1^"(a1)" = -i/3 Upsilon^"(a)"_f Pc = - N_1^"(a2)", quad Top(0) G_"(a)" = 0 . $ <eq:g1_ccRF-N>
The sum vanishes before integration: at zero external momenta @eq:g1_ccRF-sandwich is zero, identically in $MUV$, $eta$ and $ep$. The same holds for (a$'$), with $slb(k) (hat(k)^2 slb(k) - hat(k)^2 slb(k)) = 0$. The higher terms in the external momenta of @eq:g1_ccRF-G have degree $< 0$: under $Pphys$ they are $eta$ times a finite integral and vanish at $ep -> 0$.

_(b)._ $ell dot ell' = -k^2$, and the numerator of @eq:g1_ccRF-Gb becomes $hat(k)^2 k^2 Pc$; in @eq:g1_ccRF-Gbp it becomes $-hat(k)^2 k^2 Pcb$. With @eq:g1_ccRF-tens and $2\/d -> 1\/2$:
$ Top(0) G_"(b)" = ilf thin eta/ep (-1/2 Upsilon^"(b)"_f Pc) + dots, quad N_0^"(b)" = 0, \
  N_1^"(b)" = -1/2 Upsilon^"(b)"_f Pc = sum_(r = G, W) C^((r))_"ad"/4 thin wm^((r)) wm^((r)) Pc , $ <eq:g1_ccRF-Nb>
$ N_0^"(b')" = 0, quad N_1^"(b')" = 1/2 overline(Upsilon)^"(b)"_f Pcb = sum_(r = G, W) C^((r))_"ad"/4 thin wm^((r)) wm^((r)) Pcb . $ <eq:g1_ccRF-Nbp>
$N_1^"(b)"$ is the $hat(k)^2$ of the sterile insertion @eq:m-sterile, contracted with the loop momenta of the ghost vertices. In (a) the same insertion is cancelled by the second term of @eq:V-delta.

*Result.*
#result[
  $ (Dbrk_1 dot Gamma)^((1))_(gh^a gh^c overline(srcF)_f Psi_f) |_"pole" = lf 1/ep N_0, quad N_0 = N_0^"(a)" + N_0^"(b)" = 0 , $ <eq:g1_ccRF-pole>
  and the same for $gh^a gh^c overline(Psi)_f srcF_f$: no pole, evanescent or not. The finite breaking follows from @eq:g1_ccRF-Nb and the left derivatives of @eq:m-Gamma ($gh^a gh^c$ in front, $overline(srcF)$ even): in $wm^((r)) wm^((r)) = sum_(x, y in r) gh^x gh^y Tg^x Tg^y$ the coefficient of $gh^a gh^c$ is $Tg^a Tg^c - Tg^c Tg^a = i fg^(a c e) Tg^e$ (@eq:fg), so
  $ X_(1; gh^a gh^c overline(srcF)_f Psi_f) & = lf thin Pphys N_1 = i lf thin C^((r))_"ad"/4 thin fg^(a c e) Tg^e_f thin P_(c_f), \
    X_(1; gh^a gh^c overline(Psi)_f srcF_f) & = i lf thin C^((r))_"ad"/4 thin fg^(a c e) Tg^e_f thin P_(overline(c)_f) $ <eq:g1_ccRF-X>
  for $a, c$ in the same non-abelian factor $r$ carried by $f$, and zero otherwise ($a$, $c$ in different factors, or $gh^B$). Topology (a) contributes nothing; the whole result is (b), (b$'$). As a local functional,
  $ X_1 supset lf sum_f sum_(r = G, W; thin r in f) C^((r))_"ad"/4 [ N^(f r r)_(R R) + overline(N)^(f r r)_(R R) ] , $ <eq:g1_ccRF-X-op>
  and the mixed $N^(f r s)_(R R)$, $overline(N)^(f r s)_(R R)$ ($r != s$) have coefficient zero. The class contributes no singular counterterm. By @eq:g1_ccRF-others, $X_1$ vanishes in the classes $gh gh tilde(srcA)$, $gh gh tilde(srcA) A$, $gh gh gh srcC$, $gh gh srcH^dagger H$, $gh gh H^dagger srcH$, which first receive graphs with $Dbrk_1$ at two loops.
]

#figure(
  table(
    columns: 5,
    align: (left, left, center, center, center),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([$(r, s)$], [multiplets], [$N_0$], [coefficient of $N^(f r s)_(R R)$ and of $overline(N)^(f r s)_(R R)$ in $16 pi^2 X_1$], [from]),
    [$(G, G)$], [$Q, U, D$], [$0$], [$3/4 thin gs^2$], [(b), (b$'$)],
    [$(W, W)$], [$Q, L$], [$0$], [$1/2 thin gw^2$], [(b), (b$'$)],
    [$(G, W)$], [$Q$], [$0$], [$0$], [(a) vanishes],
    [$(G, B)$], [$Q, U, D$], [$0$], [$0$], [(a) vanishes],
    [$(W, B)$], [$Q, L$], [$0$], [$0$], [(a) vanishes],
  ),
  caption: [Standard-Model values, from $C^((G))_"ad" = 3 gs^2$, $C^((W))_"ad" = 2 gw^2$ (@tab:invariants). In the mixed rows the group factor of (a) is non-zero, $Upsilon^"(a)"_f = -i/2 (3 gs^2 - 2 gw^2) wm^((G)) wm^((W))$, $-3 i\/2 thin gs^2 wm^((G)) wm^((B))$, $-i gw^2 wm^((W)) wm^((B))$ (@eq:g1_ccRF-Theta), and it multiplies the vanishing Dirac integral. The lepton $E$ carries only $"U"(1)_Y$ and has no graph.],
) <tab:g1_ccRF-sm>

#check[*Mechanism of (a).* At zero ghost momentum the Feynman rule @eq:V-delta is $i V_(Dbrk_1) (k, k) = Tg^a (slh(k) Pcb - slh(k) Pc) = Tg^a (sl(k) Pcb - Pc sl(k))$ (@eq:m-tree-WI), a difference of two kinetic vertices. Between two propagators it gives $(i sl(k)\/k^2) thin i V_(Dbrk_1)(k,k) thin (i sl(k)\/k^2) = i Tg^a (Pcb (i sl(k))\/k^2 - (i sl(k))\/k^2 Pc)$, and the projectors of the neighbouring physical vertices annihilate both terms, $Pc Pcb = 0$. In (b) the $Dbrk_1$ vertex has only one internal fermion neighbour and carries the loop momentum on its ghost leg, so the mechanism does not apply.]

#check[*Independent evaluations.* The Wick signs were recomputed by explicit permutation parity. The Dirac parts of (a), (a$'$), (b), (b$'$) were recomputed for both chiralities with all external momenta $d$-dimensional: $N_0 = 0$, $Pphys N_1 = 0, 0, -1/2 Pc, +1/2 Pcb$ per unit group factor; $MUV$-independence holds and a shift of the loop momentum changes nothing. The same values follow from an explicit representation of the Clifford algebra at $d = 6, 8, 10, 12$, in which the reduced numerator of (b) is $(d-4)\/d$ times $Pc$. The group factors @eq:g1_ccRF-Theta were fitted with the explicit generators of every multiplet, with zero residual. The computer-algebra evaluation splits $ell dot ell'$ into barred and hatted parts; this changes only an $eta^2$ term and leaves $Pphys N_1$ unchanged.]

#check[*Cancellation of the pole.* The check of @sec:m-oneloop requires $bop_d delta S^((1)) = 0$ in this class, since $N_0 = 0$. By @eq:b-R and @eq:b-zeta, the $gh gh overline(srcF) psi$ components of #bop come from $bop overline(srcF)_f supset -i overline(srcF)_f wm$, $sop wm^((r)) = i wm^((r)) wm^((r))$, $sop psi_f = i wm psi_f$ and $bop srcC_a supset i overline(srcF)_f Tg^a_f psi_f$:
$ bop Q^(f r)_R |_(gh gh overline(srcF) psi) = N^(f r r)_(R R), quad quad bop Q_zeta [-1/2 fg_((r))] |_(gh gh overline(srcF) psi) = - sum_(f: thin r in f) N^(f r r)_(R R) , $ <eq:g1_ccRF-bct>
the first as in @eq:btable-fermion; the conjugate components give $overline(N)^(f r r)_(R R)$ in the same way. @sec:g0_ghost_sources gives $delta S^((1)) supset -lf (1\/ep) sum_r C^((r))_"ad"\/2 thin [Q_zeta [-1/2 fg_((r))] + sum_f Q^(f r)_R + dots]$ with equal coefficients, so $bop_d delta S^((1))$ vanishes in this class and in the conjugate class.]

#check[*Group structure and limits.* For a simple gauge group $Upsilon^"(a)"_f = 0$, and (b) alone gives $X_1 supset lf (C_"ad"\/4) integral [overline(srcF) wm wm psi + overline(psi) wm wm srcF]$. In the abelian factor there is no graph (b) ($fg^(a b c)_((B)) = 0$ and $wm^((B)) wm^((B)) = 0$), and $gh^B$ appears only on the $Dbrk_1$ vertex of (a), which vanishes. $N^(f r s)_(R R)$ and $overline(N)^(f r s)_(R R)$ appear with equal coefficients, the combination that $bop Q^(f r)_R$ produces in @eq:btable-fermion.]

#check[*Consequence for @sec:st-solution.* With @eq:g1_ccRF-bct, $bop Delta S^((1)) = -X_1$ restricted to this class reads $r[Q^(f r)_R] - r[Q_zeta [-1/2 fg_((r))]] = -C^((r))_"ad"\/4$ for every $f$ and $r = G, W$ with $r in f$. The condition is the same for every $f$, so it can be met. The mixed monomials $N^(f G B)_(R R)$, $N^(f W B)_(R R)$ come only from $bop srcC_a supset i overline(srcF)_f Tg^a_f psi_f$ in $bop Q_zeta [(srcC^A gh^A gh^B)]$ and $bop Q_zeta [(srcC^I gh^I gh^B)]$, whose coefficients must therefore vanish, as (N4) of @sec:normalisation imposes anyway. No monomial gives $N^(Q G W)_(R R)$, and $X_1$ has no component along it.]

#general[The finite part comes from the graph in which the $Dbrk_1$ ghost is internal: its momentum is a loop momentum and the sterile insertion $hat(k)^2$ is not cancelled. The graphs with an external $Dbrk_1$ ghost vanish at degree $0$ for the reason given in the first check. This holds at every order whenever both fermion neighbours of $Dbrk_1$ are physical vertices. At two loops the class receives graphs with one $Dbrk_1$ vertex and the one-loop graphs with one vertex of $bop_d delta S^((1))$ (class (c) of @sec:m-beyond). The five other source functions also start there. From two loops on, the glued products $cal(S)_"q" (C)$ are non-zero for source-dependent $C$ (@eq:m-Delta-l), so these classes must be computed at that order.]
