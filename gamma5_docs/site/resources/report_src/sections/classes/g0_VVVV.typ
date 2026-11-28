#import "../../macros.typ": *
#import "../../diagrams.typ": *

=== Gauge four-point functions $A A A A$ <sec:g0_VVVV>

External legs $A^a_mu A^b_nu A^c_rho A^d_sigma$, all at zero momentum (the pole of a function of dimension zero is momentum independent). The sum of the diagrams below is $i thin Gamma^((1))_(A^a_mu A^b_nu A^c_rho A^d_sigma)$. Every diagram has superficial degree of divergence $dg(G) = 0$. The indices $a, b, c, d$ run over all twelve gauge fields, so all gauge factors are covered at once: $G^4$, $W^4$, $B^4$ and all mixed combinations.

*Diagrams.* No internal line can mix species: a vertex with one external gauge leg and two internal lines has both internal lines of the same field (@sec:feynman-rules), and no vertex has two external legs and two internal lines of different species. The topologies are the eight of @fig:g0_VVVV-diagrams.

#figure(
  stack(dir: ttb, spacing: 5mm,
  diagrow(
    oneloop(
      ((legs: (("V", $a$, 0),)), (legs: (("V", $b$, 0),)), (legs: (("V", $c$, 0),)), (legs: (("V", $d$, 0),))),
      (("F", none, 1), ("F", none, 1), ("F", none, 1), ("F", $Psi_f$, 1)), caption: [(a) fermion box]),
    oneloop(
      ((legs: (("V", $a$, 0),)), (legs: (("V", $b$, 0),)), (legs: (("V", $c$, 0),)), (legs: (("V", $d$, 0),))),
      (("S", none, 1), ("S", none, 1), ("S", none, 1), ("S", $H$, 1)), caption: [(b) Higgs box]),
    oneloop(
      ((legs: (("V", $a$, 0), ("V", $b$, 0))), (legs: (("V", $c$, 0),)), (legs: (("V", $d$, 0),))),
      (("S", none, 1), ("S", none, 1), ("S", $H$, 1)), caption: [(c) Higgs triangle]),
    oneloop(
      ((legs: (("V", $a$, 0), ("V", $b$, 0))), (legs: (("V", $c$, 0), ("V", $d$, 0)))),
      (("S", none, 1), ("S", $H$, 1)), caption: [(d) Higgs bubble]),
  ),
  diagrow(
    oneloop(
      ((legs: (("V", $a$, 0),)), (legs: (("V", $b$, 0),)), (legs: (("V", $c$, 0),)), (legs: (("V", $d$, 0),))),
      (("V", none, 0), ("V", none, 0), ("V", none, 0), ("V", none, 0)), caption: [(e) gauge box]),
    oneloop(
      ((legs: (("V", $a$, 0), ("V", $b$, 0))), (legs: (("V", $c$, 0),)), (legs: (("V", $d$, 0),))),
      (("V", none, 0), ("V", none, 0), ("V", none, 0)), caption: [(f) gauge triangle]),
    oneloop(
      ((legs: (("V", $a$, 0), ("V", $b$, 0))), (legs: (("V", $c$, 0), ("V", $d$, 0)))),
      (("V", none, 0), ("V", none, 0)), caption: [(g) gauge bubble]),
    oneloop(
      ((legs: (("V", $a$, 0),)), (legs: (("V", $b$, 0),)), (legs: (("V", $c$, 0),)), (legs: (("V", $d$, 0),))),
      (("G", none, 1), ("G", none, 1), ("G", none, 1), ("G", $gh$, 1)), caption: [(h) ghost box]),
  )),
  caption: [Topologies of $Gamma^((1))_(A A A A)$, each drawn once. External legs are labelled by their gauge index; the legs $(a, mu)$, $(b, nu)$, $(c, rho)$, $(d, sigma)$ are distributed over the vertices in all inequivalent ways (@tab:g0_VVVV-diagrams).],
) <fig:g0_VVVV-diagrams>

The group factors are written with ordered traces. A _word_ $w = (w_1 w_2 w_3 w_4)$ is a cyclic ordering of the four legs; up to cyclic rotation there are six,
$ cal(W) = {(a b c d), (a b d c), (a c b d), (a c d b), (a d b c), (a d c b)}, quad
  Tr_cal(R) [Tg^((w))] = Tr_cal(R) [Tg^(w_1) Tg^(w_2) Tg^(w_3) Tg^(w_4)] , $ <eq:g0_VVVV-words>
where the representation $cal(R)$ is a fermion multiplet $f$, the Higgs doublet $H$, or the adjoint representation, with the coupling-weighted generators $(Tg^a_"adj")_(x y) = -i fg^(a x y)$. These obey $[Tg^a_"adj", Tg^b_"adj"] = i fg^(a b c) Tg^c_"adj"$ and $Tr_"adj" (Tg^a_"adj" Tg^b_"adj") = C_"ad"^((a)) delta^(a b)$ (@eq:invariants); the symbol $Tg_"adj"$ and the trace $Tr_"adj"$ are used only in this subsection. Inside a gauge or ghost loop the structure constants of the vertices are rewritten as
$ fg^(a x y) = i (Tg^a_"adj")_(x y), quad fg^(a x e) fg^(b y e) = (Tg^a_"adj" Tg^b_"adj")_(x y), quad fg^(a b e) fg^(x y e) = ([Tg^a_"adj", Tg^b_"adj"])_(x y) , $ <eq:g0_VVVV-fadj>
with $x$ the internal line entering the vertex and $y$ the one leaving it in the direction of the loop momentum.

#figure(
  table(
    columns: (auto, 1.25fr, auto, auto, auto, 1.9fr),
    align: (center, left, center, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([], [internal lines], [number], [sym.], [sign], [group factor]),
    [(a)], [$Psi_f$, $f in {Q, U, D, L, E}$], [$6$ per $f$], [$1$], [$-1$], [$Tr_f [Tg^((w))]$, $w in cal(W)$],
    [(b)], [$H$], [$6$], [$1$], [$+1$], [$Tr_H [Tg^((w))]$, $w in cal(W)$],
    [(c)], [$H$, seagull on a pair $(i j)$], [$6 times 2$], [$1$], [$+1$], [$Tr_H [{Tg^(i), Tg^(j)} Tg^(k) Tg^(l)]$, both orders of $k l$],
    [(d)], [$H$, seagulls on $(i j)$, $(k l)$], [$3$], [$1$], [$+1$], [$Tr_H [{Tg^(i), Tg^(j)} {Tg^(k), Tg^(l)}]$],
    [(e)], [$A^x$; four $A A A$], [$3$], [$1$], [$+1$], [$Tr_"adj" [Tg^((w))]$, $w$ and its reverse],
    [(f)], [$A^x$; $A A A A$ on $(i j)$], [$6$], [$1$], [$+1$], [$Tr_"adj"$ of $Tg^(i) Tg^(j) Tg^(k) Tg^(l)$, $Tg^(j) Tg^(i) Tg^(k) Tg^(l)$],
    [(g)], [$A^x$; $A A A A$ on $(i j)$, $(k l)$], [$3$], [$1\/2$], [$+1$], [$Tr_"adj"$ of $(Tg^(i) Tg^(j) "or" Tg^(j) Tg^(i))(Tg^(k) Tg^(l) "or" Tg^(l) Tg^(k))$],
    [(h)], [$gh^x$, $agh^x$], [$6$], [$1$], [$-1$], [$Tr_"adj" [Tg^((w))]$, $w in cal(W)$],
  ),
  caption: [Field assignments of @fig:g0_VVVV-diagrams for four distinct external legs; the sign is that of a closed fermion or ghost loop. In (a), (b), (h) the two orientations of the complex line are distinct diagrams, so all six words occur. The real gauge loops (e), (f) are written below as $1\/2$ times the sum over both orientations, which counts each diagram once; the bubble (g) has the symmetry factor $1\/2$. Internal gauge and ghost lines belong to the non-abelian factor of the external legs; the Higgs and fermion loops couple to all three factors.],
) <tab:g0_VVVV-diagrams>

*Integrands.* All external momenta vanish, so every internal line carries the loop momentum $k$, oriented along the loop (along the arrow for $Psi_f$, $H$, $gh$). For the word $w$ the vertex of leg $w_i$ is the $i$-th vertex along the loop; $mu_(w_i)$ denotes its Lorentz index ($mu_a = mu$, $mu_b = nu$, $mu_c = rho$, $mu_d = sigma$). From @eq:propagators, @eq:V-ffv, @eq:V-ssv, @eq:V-vvv, @eq:V-vvvv and @eq:V-ghost, with the Dirac and group matrices ordered against the fermion flow:
$ G^((w))_"(a)" = - intk tr[ gbar^(mu_(w_1)) P_(c_f) sl(k) gbar^(mu_(w_2)) P_(c_f) sl(k) gbar^(mu_(w_3)) P_(c_f) sl(k) gbar^(mu_(w_4)) P_(c_f) sl(k) ] / (k^2)^4 thin Tr_f [Tg^((w))] , $ <eq:g0_VVVV-int-a>
where the factors $i^4$ of the vertices and $i^4$ of the propagators give $1$. For the Higgs loop the cubic vertex is $i (2 k)^mu Tg^a_H$ ($p = p' = k$) and the propagator $i\/(k^2 - m^2)$:
$ G^((w))_"(b)" & = 16 intk (k^(mu_(w_1)) k^(mu_(w_2)) k^(mu_(w_3)) k^(mu_(w_4)))/(k^2 - m^2)^4 thin Tr_H [Tg^((w))], \
  G_"(c)" & = -4 g^(mu_i mu_j) intk (k^(mu_k) k^(mu_l))/(k^2 - m^2)^3 thin Tr_H [{Tg^(i), Tg^(j)} Tg^(k) Tg^(l)] , \
  G_"(d)" & = g^(mu_i mu_j) g^(mu_k mu_l) intk 1/(k^2 - m^2)^2 thin Tr_H [{Tg^(i), Tg^(j)} {Tg^(k), Tg^(l)}] . $ <eq:g0_VVVV-int-bcd>
In the gauge loop the Feynman rule of the $A A A$ vertex with external leg $(a, alpha)$, entering line $(x, beta)$ of momentum $k$ and leaving line $(y, gamma)$ of incoming momentum $-k$ is, by @eq:V-vvv and @eq:g0_VVVV-fadj, $i (Tg^a_"adj")_(x y) cal(V)^(alpha beta gamma)(0, k, -k)$, with $cal(V)$ the bracket of @eq:V-vvv as in @sec:g0_VVV:
$ cal(V)^(alpha beta gamma)(0, k, -k) = - g^(alpha beta) k^gamma + 2 g^(beta gamma) k^alpha - g^(gamma alpha) k^beta . $ <eq:g0_VVVV-V3>
The Feynman rule of the $A A A A$ vertex with external legs $(a, alpha)$, $(b, beta)$, entering line $(x, gamma)$ and leaving line $(y, delta)$ is, by @eq:V-vvvv and @eq:g0_VVVV-fadj,
$ -i [ (Tg^a_"adj" Tg^b_"adj")_(x y) Theta^(alpha beta gamma delta) + (Tg^b_"adj" Tg^a_"adj")_(x y) Theta^(beta alpha gamma delta) ], quad
  Theta^(alpha beta gamma delta) = g^(alpha beta) g^(gamma delta) + g^(alpha gamma) g^(beta delta) - 2 g^(alpha delta) g^(beta gamma) . $ <eq:g0_VVVV-V4>
With the propagator $-i g_(gamma beta) \/ k^2$ joining the leaving index $gamma$ of one vertex to the entering index $beta$ of the next, the factors are $i^4 (-i)^4 = 1$ in (e), $(-i) dot i^2 dot (-i)^3 = -1$ in (f) and $(-i)^2 (-i)^2 = 1$ in (g):
$ G^((w))_"(e)" = 1/2 intk (cal(C)^(mu_(w_1) mu_(w_2) mu_(w_3) mu_(w_4))_"(e)")/(k^2)^4 thin Tr_"adj" [Tg^((w))] , \
  cal(C)^(mu_1 mu_2 mu_3 mu_4)_"(e)" = cal(V)^(mu_1 beta_1 gamma_1) g_(gamma_1 beta_2) cal(V)^(mu_2 beta_2 gamma_2) g_(gamma_2 beta_3) cal(V)^(mu_3 beta_3 gamma_3) g_(gamma_3 beta_4) cal(V)^(mu_4 beta_4 gamma_4) g_(gamma_4 beta_1) , $ <eq:g0_VVVV-int-e>
with every $cal(V)$ taken at $(0, k, -k)$. In (f) the vertex @eq:g0_VVVV-V4 sits on two adjacent legs of the word, in (g) on both pairs of adjacent legs:
$ G^((w))_"(f)" = -1/2 sum_(s = 1)^4 intk (cal(C)_"(f)"^(mu_(w_s) mu_(w_(s+1)) mu_(w_(s+2)) mu_(w_(s+3))))/(k^2)^3 thin Tr_"adj" [Tg^((w))], \
  cal(C)^(mu_1 mu_2 mu_3 mu_4)_"(f)" = Theta^(mu_1 mu_2 gamma delta) g_(delta beta_3) cal(V)^(mu_3 beta_3 gamma_3) g_(gamma_3 beta_4) cal(V)^(mu_4 beta_4 gamma_4) g_(gamma_4 gamma) , $ <eq:g0_VVVV-int-f>
$ G^((w))_"(g)" = 1/2 sum_(s = 1)^2 intk (cal(C)_"(g)"^(mu_(w_s) mu_(w_(s+1)) mu_(w_(s+2)) mu_(w_(s+3))))/(k^2)^2 thin Tr_"adj" [Tg^((w))], \
  cal(C)^(mu_1 mu_2 mu_3 mu_4)_"(g)" = Theta^(mu_1 mu_2 gamma delta) g_(delta gamma') Theta^(mu_3 mu_4 gamma' delta') g_(delta' gamma) , $ <eq:g0_VVVV-int-g>
with the positions $s$ counted cyclically ($w_5 = w_1$, ...). The sum over $s$ collects, for the word $w$, the diagrams in which the quartic vertex carries the adjacent legs $(w_s, w_(s+1))$ (in (g) also $(w_(s+2), w_(s+3))$), with the term of @eq:g0_VVVV-V4 whose matrices stand in the order of $w$; with the factor $1\/2$ for the two orientations of the real loop in (e), (f) and the symmetry factor $1\/2$ in (g), every diagram of @tab:g0_VVVV-diagrams is counted once. In the ghost loop the ghost--gauge vertex @eq:V-ghost, with the ghost line $c$ entering and the antighost line $e$ leaving in the direction of $k$ and outgoing antighost momentum $k$, is $-fg^(e a c) k^mu = -i (Tg^a_"adj")_(c e) k^mu$ at every vertex; the factors $(-i)^4$ of the vertices and $i^4$ of the propagators give $1$:
$ G^((w))_"(h)" = - intk (k^(mu_(w_1)) k^(mu_(w_2)) k^(mu_(w_3)) k^(mu_(w_4)))/(k^2)^4 thin Tr_"adj" [Tg^((w))] . $ <eq:g0_VVVV-int-h>

*$Rop$-operation.* At one loop $Rbar G = G$. Since $dg(G) = 0$, $Top(0)$ keeps the leading term of @eq:uv-expansion only: every denominator becomes $k^2 - MUV^2$, and $m^2$ (degree two) drops out. All resulting integrals are $intk (k^2)^n \/ (k^2 - MUV^2)^(n+2)$ with the universal pole $ilf \/ ep$ (@eq:logpole). With @eq:tensor-red, $k^alpha k^beta -> k^2 g^(alpha beta) \/ d$ and $k^alpha k^beta k^gamma k^delta -> (k^2)^2 (g^(alpha beta) g^(gamma delta) + g^(alpha gamma) g^(beta delta) + g^(alpha delta) g^(beta gamma)) \/ (d (d+2))$; in $N_0$ these factors, and the traces $g^alpha_(space alpha) = d$ of the bosonic loops, are taken at $d = 4$, since their $O(ep)$ parts multiply the pole and give finite terms. No hatted object appears: the fermion numerators are four-dimensional by @eq:phys-prop, the bosonic ones are unsplit, and the external momenta vanish. Hence $N_1 = 0$ for every diagram, and $Kop = Kop_eta$.

The results are expressed through two metric structures per word, with $g^(a b) equiv g^(mu nu)$, $g^(a c) equiv g^(mu rho)$, and so on (local notation):
$ cal(G)_"adj" [w] = g^(w_1 w_2) g^(w_3 w_4) + g^(w_1 w_4) g^(w_2 w_3), quad cal(G)_"opp" [w] = g^(w_1 w_3) g^(w_2 w_4) , $ <eq:g0_VVVV-G>
and $overline(cal(G))_"adj"$, $overline(cal(G))_"opp"$ for the same products of $gmb$.

_Fermion box._ $P_(c_f) sl(k) gbar^nu P_(c_f) = P_(c_f) sl(k) P_(overline(c)_f) gbar^nu = slb(k) P_(overline(c)_f) gbar^nu$, so the trace of @eq:g0_VVVV-int-a is $tr[gbar^(mu_1) slb(k) gbar^(mu_2) slb(k) gbar^(mu_3) slb(k) gbar^(mu_4) slb(k) P_(overline(c)_f)]$ (indices $mu_i = mu_(w_i)$), with only barred matrices. The reduction $bar(k)^4 -> (k^2)^2 (gmb gmb + gmb gmb + gmb gmb) \/ 24$ and @eq:contr2 ($gbar^alpha gbar^nu gbar_alpha = -2 gbar^nu$, $gbar^alpha gbar^nu gbar^beta gbar^rho gbar_alpha = -2 gbar^rho gbar^beta gbar^nu$, $gbar^beta gbar^nu gbar^sigma gbar_beta = 4 gmb^(nu sigma)$) give for the three pairings
$ gbar^(mu_1) gbar^alpha gbar^(mu_2) gbar^beta gbar^(mu_3) gbar^gamma gbar^(mu_4) gbar^delta (gmb_(alpha beta) gmb_(gamma delta) + gmb_(alpha gamma) gmb_(beta delta) + gmb_(alpha delta) gmb_(beta gamma))
  \ = 4 gbar^(mu_1) gbar^(mu_2) gbar^(mu_3) gbar^(mu_4) - 8 gmb^(mu_2 mu_4) gbar^(mu_1) gbar^(mu_3) + 4 gbar^(mu_1) gbar^(mu_4) gbar^(mu_3) gbar^(mu_2) . $ <eq:g0_VVVV-dirac>
With $P_(overline(c)_f) = (1 plus.minus g5)\/2$ the $g5$ part is $plus.minus 1/2 (-4 i)(4 epsilon^(mu_1 mu_2 mu_3 mu_4) + 4 epsilon^(mu_1 mu_4 mu_3 mu_2)) = 0$ by @eq:tr5, separately for each word and for either chirality. The remaining trace is $16 overline(cal(G))_"adj" [w] - 32 overline(cal(G))_"opp" [w]$, and with the fermion-loop sign and the factor $1\/24$,
$ N_0^"(a)" = sum_f sum_(w in cal(W)) Tr_f [Tg^((w))] ( -2/3 overline(cal(G))_"adj" [w] + 4/3 overline(cal(G))_"opp" [w] ) . $ <eq:g0_VVVV-N-a>

_Higgs loops._ In (b) $16 dot 1\/24$ gives $2\/3 (cal(G)_"adj" + cal(G)_"opp")$ per word. In (c) $-4 dot 1\/4 = -1$ multiplies $g^(mu_i mu_j) g^(mu_k mu_l)$; expanding the anticommutators, the word $w$ receives four terms (seagull on $(w_1 w_2)$, $(w_3 w_4)$, $(w_2 w_3)$ or $(w_4 w_1)$), all with adjacent pairs: $-2 cal(G)_"adj" [w]$. In (d) the word $w$ receives the pairings $(w_1 w_2)(w_3 w_4)$ and $(w_4 w_1)(w_2 w_3)$: $+cal(G)_"adj" [w]$.

_Gauge and ghost loops._ The chains of @eq:g0_VVVV-int-e to @eq:g0_VVVV-int-g are contracted at $d = 4$, with $g_(i j) equiv g^(mu_i mu_j)$ and $k_i equiv k^(mu_i)$ (local notation):
$ cal(C)_"(e)" & = (k^2)^2 (g_(1 2) g_(3 4) + g_(1 4) g_(2 3)) + 3 k^2 (g_(1 2) k_3 k_4 + g_(2 3) k_4 k_1 + g_(3 4) k_1 k_2 + g_(4 1) k_2 k_3) + 34 thin k_1 k_2 k_3 k_4 , \
  cal(C)_"(f)" & = k^2 (2 g_(1 2) g_(3 4) - 2 g_(1 3) g_(2 4) + g_(1 4) g_(2 3)) + 6 g_(1 2) k_3 k_4 - g_(3 4) k_1 k_2 + 3 g_(1 4) k_2 k_3 + 3 g_(2 3) k_1 k_4 , \
  cal(C)_"(g)" & = 2 g_(1 2) g_(3 4) - 4 g_(1 3) g_(2 4) + 5 g_(1 4) g_(2 3) . $ <eq:g0_VVVV-chains>
The last line follows by hand from the nine products of the three terms of $Theta$: with $g^gamma_(space gamma) = 4$, $g_(1 2) g_(3 4)$ collects $4 + 1 - 2 + 1 - 2 = 2$, $g_(1 3) g_(2 4)$ collects $-2 - 2 = -4$ and $g_(1 4) g_(2 3)$ collects $1 + 4 = 5$. After $Top(0)$ the integrals are $intk (k^2)^j \/ (k^2 - MUV^2)^(j+2)$ with $j = 2, 1, 0$, and the reductions above give, per position $s$, the residues
$ cal(C)_"(e)" -> 47/12 (g_(1 2) g_(3 4) + g_(1 4) g_(2 3)) + 17/12 g_(1 3) g_(2 4), \
  cal(C)_"(f)" -> 13/4 g_(1 2) g_(3 4) - 2 g_(1 3) g_(2 4) + 5/2 g_(1 4) g_(2 3), quad
  cal(C)_"(g)" -> cal(C)_"(g)" . $ <eq:g0_VVVV-chains-red>
For (e): $1 + 3 dot 2\/4 + 34\/24 = 47\/12$ and $34\/24 = 17\/12$; times $1\/2$ this is $(x, y) = (47\/24, 17\/24)$. For (f) the sum over the four cyclic positions gives $2 (13/4 + 5/2) = 23\/2$ for each of $g_(1 2) g_(3 4)$ and $g_(1 4) g_(2 3)$ and $-8$ for $g_(1 3) g_(2 4)$; times $-1\/2$ this is $(-23\/4, 4)$. For (g) the two positions give $2 + 5 = 7$ for each adjacent pairing and $-8$ for the opposite one; times $1\/2$ this is $(7\/2, -4)$. In (h) the numerator is totally symmetric: $-1\/24 (cal(G)_"adj" + cal(G)_"opp")$ per word.

All topologies have the form $N_0^X = sum_(w in cal(W)) Tr_cal(R) [Tg^((w))] (x_X cal(G)_"adj" [w] + y_X cal(G)_"opp" [w])$, with the coefficients of @tab:g0_VVVV-N. Then $lf N_0^X \/ ep$ is the vertex function of the local functional
$ lf 1/ep integral dif^d x thin ( x_X/2 thin Tr_cal(R) [Am_mu Am^mu Am_nu Am^nu] + y_X/4 thin Tr_cal(R) [Am_mu Am_nu Am^mu Am^nu] ) , $ <eq:g0_VVVV-op>
since the functional derivative of $Tr_cal(R) [Am_mu Am^mu Am_nu Am^nu]$ gives $2 cal(G)_"adj" [w]$ and that of $Tr_cal(R) [Am_mu Am_nu Am^mu Am^nu]$ gives $4 cal(G)_"opp" [w]$ per word (four cyclic rotations each). Here $Am_mu = Tg^a A^a_mu$ in the representation $cal(R)$ (@eq:matrix-fields), with $Tg_"adj"$ for $cal(R) = "adj"$, and $bar(A)$ in place of $A$ for (a).

#figure(
  table(
    columns: 5,
    align: (left, center, center, center, center),
    stroke: (x, y) => if y == 0 or y == 1 or y == 5 { (bottom: 0.6pt) } else { none },
    table.header([topology], [$cal(R)$], [$x_X$ (of $cal(G)_"adj"$)], [$y_X$ (of $cal(G)_"opp"$)], [$x_X + y_X\/2$]),
    [(a) fermion box (barred)], [$f$], [$-2\/3$], [$4\/3$], [$0$],
    [(b) Higgs box], [$H$], [$2\/3$], [$2\/3$], [$1$],
    [(c) Higgs triangle], [$H$], [$-2$], [$0$], [$-2$],
    [(d) Higgs bubble], [$H$], [$1$], [$0$], [$1$],
    [(b)+(c)+(d)], [$H$], [$-1\/3$], [$2\/3$], [$0$],
    [(e) gauge box], [adj], [$47\/24$], [$17\/24$], [$37\/16$],
    [(f) gauge triangle], [adj], [$-23\/4$], [$4$], [$-15\/4$],
    [(g) gauge bubble], [adj], [$7\/2$], [$-4$], [$3\/2$],
    [(h) ghost box], [adj], [$-1\/24$], [$-1\/24$], [$-1\/16$],
    [(e)+(f)+(g)+(h)], [adj], [$-1\/3$], [$2\/3$], [$0$],
  ),
  caption: [Residues $N_0$ per topology in the form $sum_w Tr_cal(R) [Tg^((w))] (x_X cal(G)_"adj" [w] + y_X cal(G)_"opp" [w])$, summed over all diagrams of the topology (symmetry factors included). The last column vanishes exactly when the functional @eq:g0_VVVV-op is a pure commutator term.],
) <tab:g0_VVVV-N>

*Result.* For each species the sum has $x = -y\/2$, so @eq:g0_VVVV-op is proportional to $Tr_cal(R) [Am_mu Am^mu Am_nu Am^nu] - Tr_cal(R) [Am_mu Am_nu Am^mu Am^nu] = -1/2 Tr_cal(R) ([Am_mu, Am_nu][Am^mu, Am^nu])$. With $[Am_mu, Am_nu] = i fg^(a b e) A^a_mu A^b_nu Tg^e$ and @eq:invariants,
$ Tr_cal(R) ([Am_mu, Am_nu][Am^mu, Am^nu]) = 4 sum_(r = G, W) S^((r))_(cal(R)) thin L^((r))_4 [A], quad
  L^((r))_4 [A] = -1/4 fg^(a b e) fg^(c d e) (A^a dot A^c)(A^b dot A^d) , $ <eq:g0_VVVV-L4>
with $a, dots, e$ in the factor $r$ and $S^((r))_(cal(R))$ the index of @eq:invariants for $cal(R)$, $S^((r))_"adj" = C^((r))_"ad"$; $L^((r))_4$ (local notation) is the quartic part of $-1/4 F^a_(mu nu) F^(a mu nu)$ restricted to the factor $r$, and $L^((r))_4 [bar(A)]$ the same with $bar(A)$. Generators of different factors commute, so $[Am_mu, Am_nu]$ is a sum of commutators within one non-abelian factor, and the cross terms $Tr_cal(R) ([Am^((G))_mu, Am^((G))_nu][Am^((W) mu), Am^((W) nu)])$ are proportional to $Tr(T^E) Tr(tau^K) = 0$. Hence the functions $B B B B$, $W W B B$, $G G B B$, $G G W W$, $G G G B$ and $W W W B$ vanish after the sum over words, although single topologies do not (e.g. the Higgs loops in $B B B B$ and $W W B B$, @tab:g0_VVVV-topologies); the remaining combinations ($G G G W$, $G W W W$, $G B B B$, $W B B B$, $G G W B$, $G W W B$, $G W B B$) vanish diagram by diagram because a single generator of a factor is traceless.

#result[
  $ Gamma^((1))_(A^a_mu A^b_nu A^c_rho A^d_sigma) |_"pole" = lf 1/ep N_0, quad
    N_0 = -1/3 sum_(w in cal(W)) [ 2 sum_f Tr_f [Tg^((w))] (overline(cal(G))_"adj" - 2 overline(cal(G))_"opp") [w] \
      + (Tr_H [Tg^((w))] + Tr_"adj" [Tg^((w))]) (cal(G)_"adj" - 2 cal(G)_"opp") [w] ] $ <eq:g0_VVVV-pole-w>
  $ N_0 = sum_(r = G, W) [ z^((r))_"ferm" thin overline(Gamma)^((0, r))_(A A A A) + z^((r))_"bos" thin Gamma^((0, r))_(A A A A) ], quad
    z^((r))_"ferm" = 2/3 S^((r)), quad z^((r))_"bos" = 1/3 (S^((r))_H + C^((r))_"ad") , $ <eq:g0_VVVV-pole>
  with $S^((r)) = sum_f S^((r))_f = S^((r))_"lh" + S^((r))_"rh"$ (@eq:SLR); $z^((r))_"ferm"$ and $z^((r))_"bos"$ are the fermion-loop and bosonic-loop coefficients of the tree vertex function, as in @sec:g0_VV and @sec:g0_VVV. Here $Gamma^((0, r))_(A A A A) = -[fg^(a b e) fg^(c d e) (g^(mu rho) g^(nu sigma) - g^(mu sigma) g^(nu rho)) + fg^(a c e) fg^(b d e) (g^(mu nu) g^(rho sigma) - g^(mu sigma) g^(nu rho)) + fg^(a d e) fg^(b c e) (g^(mu nu) g^(rho sigma) - g^(mu rho) g^(nu sigma))]$ is the tree vertex function (@eq:V-vvvv divided by $i$) with all four legs in the factor $r$, and $overline(Gamma)^((0, r))_(A A A A)$ the same with $gmb$. In matrix-valued form,
  $ delta S^((1)) supset & - lf 1/ep integral dif^d x [ 1/6 sum_f Tr_f ([bar(Am)_mu, bar(Am)_nu][bar(Am)^mu, bar(Am)^nu]) \
    & + 1/12 Tr_H ([Am_mu, Am_nu][Am^mu, Am^nu]) + 1/12 Tr_"adj" ([Am_mu, Am_nu][Am^mu, Am^nu]) ] $ <eq:g0_VVVV-ct-matrix>
  $ delta S^((1)) supset - lf 1/ep sum_(r = G, W) integral dif^d x [ z^((r))_"bos" thin L^((r))_4 [A] + z^((r))_"ferm" thin L^((r))_4 [bar(A)] ] . $ <eq:g0_VVVV-ct>
  The bosonic loops multiply the $d$-dimensional quartic term $L^((r))_4 [A]$ with unsplit indices, the fermion loops the barred term $L^((r))_4 [bar(A)]$. In the basis of @sec:bases, with $cal(Q)^((r))_"YM" = integral dif^d x thin L^((r))_4 [A]$ of @eq:QYM continued with unsplit $d$-dimensional indices (@sec:evanescent) and $tilde(E)^((r))_F |_(A A A A) = integral dif^d x (L^((r))_4 [bar(A)] - L^((r))_4 [A])$ the part of $tilde(E)^((r))_F$ (@eq:evanescent) quartic in the gauge field,
  $ delta S^((1)) supset sum_(r = G, W) [ delta Z^((r))_(4 A) thin cal(Q)^((r))_"YM" - lf 1/ep z^((r))_"ferm" thin tilde(E)^((r))_F |_(A A A A) ], \
    16 pi^2 ep thin delta Z^((r))_(4 A) = - z^((r))_"bos" - z^((r))_"ferm" = -2/3 S^((r)) - 1/3 S^((r))_H - 1/3 C^((r))_"ad" . $ <eq:g0_VVVV-Z4> In terms of the constants of the table in @sec:st-solution, $delta Z^((r))_(4 A) = 2 delta Z_(c_r) + 2 delta Z^((r))_A$ (a symbol of this subsection) is the coefficient of $cal(Q)^((r))_"YM"$ in the assembled counterterm. The coefficient $-z^((r))_"ferm" = -2/3 S^((r))$ of $tilde(E)^((r))_F$ is the one of @sec:g0_VV and @sec:g0_VVV.
]

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, 1fr),
    align: (left, center, center, center, center, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([factor], [$z_"ferm"$], [$S_H\/3$], [$C_"ad"\/3$], [$z_"bos" + z_"ferm" = -16 pi^2 ep thin delta Z_(4 A)$], [$16 pi^2 ep thin delta S^((1))$ (physical)]),
    [$G$], [$4/3 gs^2$], [$0$], [$gs^2$], [$7/3 gs^2$], [$7/24 Q_(G 4 a) - 7/12 Q_(G 4 b) + 7/8 Q_(G 4 c)$],
    [$W$], [$4/3 gw^2$], [$1/6 gw^2$], [$2/3 gw^2$], [$13/6 gw^2$], [$13/24 Q_(W 4 a) - 13/24 Q_(W 4 b)$],
    [$B$, mixed], [$0$], [$0$], [$0$], [$0$], [$0$],
  ),
  caption: [Standard-Model values of @eq:g0_VVVV-pole and @eq:g0_VVVV-ct: fermion part $z^((r))_"ferm"$, Higgs and gauge-plus-ghost parts of $z^((r))_"bos"$, and their sum. The last column uses $cal(Q)^((G))_"YM" = gs^(-2) (-1/8 Q_(G 4 a) + 1/4 Q_(G 4 b) - 3/8 Q_(G 4 c))$ and $cal(Q)^((W))_"YM" = gw^(-2) (-1/4 Q_(W 4 a) + 1/4 Q_(W 4 b))$ (@eq:QYM). The coefficient of $tilde(E)^((r))_F |_(A A A A)$ in $16 pi^2 ep thin delta S^((1))$ is $-z_"ferm"$: $-4/3 gs^2$ ($G$) and $-4/3 gw^2$ ($W$). "Mixed" covers $Q^((r s))_(a,b)$, $Q_(G G G B)$, $Q_(B 4)$, $Q^((r))_(4 epsilon)$: all vanish.],
) <tab:g0_VVVV-values>

#figure(
  table(
    columns: 5,
    align: (left, left, left, left, left),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { none },
    table.header([topology], [$G^4$: $(Q_(G 4 a), Q_(G 4 b), Q_(G 4 c))$], [$W^4$: $(Q_(W 4 a), Q_(W 4 b))$], [$W^2 B^2$: $(Q^((W B))_a, Q^((W B))_b)$], [$B^4$: $Q_(B 4)$]),
    [(a), all $f$], [$(-1/6, 1/3, -1/2)$], [$(-1/3, 1/3)$], [$0$], [$0$],
    [(b)], [$0$], [$(1/48, 1/24)$], [$(1/8, 1/4)$], [$1/16$],
    [(c)], [$0$], [$(-1/8, 0)$], [$(-1/4, -1/2)$], [$-1/8$],
    [(d)], [$0$], [$(1/16, 0)$], [$(1/8, 1/4)$], [$1/16$],
    [(e)], [$(131/128, 91/64, 15/32)$], [$(47/48, 4/3)$], [$0$], [$0$],
    [(f)], [$(-21/8, -3/8, -117/32)$], [$(-23/8, -7/8)$], [$0$], [$0$],
    [(g)], [$(3/2, -3/4, 45/16)$], [$(7/4, -1/4)$], [$0$], [$0$],
    [(h)], [$(-3/128, -3/64, 0)$], [$(-1/48, -1/24)$], [$0$], [$0$],
    [sum], [$(-7/24, 7/12, -7/8)$], [$(-13/24, 13/24)$], [$0$], [$0$],
  ),
  caption: [The pole $16 pi^2 ep thin Gamma^((1))|_"pole"$ per topology as a local functional in the monomials of @tab:Qbasis (couplings inside the monomials; fields $bar(A)$ in (a), $d$-dimensional $A$ in (b)--(h)). The traces of four generators are reduced with $Tr_"adj" [Tg^A Tg^B Tg^C Tg^D] = gs^4 [delta^(A B) delta^(C D) + delta^(A D) delta^(B C) + 3/4 (d^(A B E) d^(C D E) - d^(A C E) d^(B D E) + d^(A D E) d^(B C E))]$ for $"SU"(3)$, $Tr_"adj" [Tg^I Tg^J Tg^K Tg^L] = gw^4 (delta^(I J) delta^(K L) + delta^(I L) delta^(J K))$ and $Tr_H [Tg^I Tg^J Tg^K Tg^L] = gw^4 \/ 8 thin (delta^(I J) delta^(K L) - delta^(I K) delta^(J L) + delta^(I L) delta^(J K))$ for $"SU"(2)$ (checked with explicit matrices), and the $d d$ identity of @sec:Qbasis; e.g. $Tr_"adj" [Am_mu Am^mu Am_nu Am^nu] = Q_(G 4 a) + Q_(G 4 b) + 3/4 Q_(G 4 c)$ and $Tr_"adj" [Am_mu Am_nu Am^mu Am^nu] = 1/4 Q_(G 4 a) + 5/2 Q_(G 4 b) - 3/2 Q_(G 4 c)$ for gluons. The fermion boxes are free of mixed terms separately for each multiplet; the counterterm is minus the last row.],
) <tab:g0_VVVV-topologies>

#check[*Auxiliary mass, Higgs mass, $N_1$.* Every term of $Top(0) G$ is a logarithmically divergent vacuum integral, whose pole does not depend on $MUV$; $m^2$ does not enter. No hatted object appears, so $N_1 = 0$, as required for a ghost-number-0 graph (@sec:m-K).]

#check[*Chirality and $epsilon$ terms.* The $g5$ trace of the fermion box vanishes word by word (@eq:g0_VVVV-dirac), so left- and right-handed multiplets contribute identically and the result depends only on $S^((r)) = sum_f S^((r))_f$. No $epsilon^(mu nu rho sigma)$ structure appears, in agreement with the absence of a totally antisymmetric four-gauge tensor without an abelian leg (@sec:Qbasis); the $Q^((r))_(4 epsilon)$ coefficients vanish.]

#check[*Abelian limits.* Each species sums separately to a commutator, so every function with a $B$ leg vanishes. For the fermion box this is the finiteness of the four-photon amplitude in QED; for the Higgs loop it is the cancellation $1/16 - 1/8 + 1/16 = 0$ of box, triangles and bubbles in scalar QED (@tab:g0_VVVV-topologies, last column).]

#check[*Invariance and other classes.* The physical part of @eq:g0_VVVV-ct is proportional to $cal(Q)^((r))_"YM"$; the non-invariant combinations of $Q_(G 4 a,b,c)$, $Q_(W 4 a,b)$ and all mixed monomials are absent, although single topologies contain them. In the closed operators of @eq:closed, $cal(Q)^((r))_"YM"$ enters $cal(O)^((r))_A$ with weight $4$ and $cal(O)^((r))_g$ with weight $2$; with the four-dimensional counterterms $delta Z^((r))_A thin Q^((r))_T$, $delta Z^((r))_(3 A) thin Q^((r))_(3 f)$ and $delta Z^((r))_(4 A) thin cal(Q)^((r))_"YM"$ of the two-, three- and four-point functions, invariance requires $delta Z^((r))_(4 A) = 2 delta Z^((r))_(3 A) - delta Z^((r))_A$. With @sec:g0_VV and @sec:g0_VVV, $16 pi^2 ep (2 delta Z^((r))_(3 A) - delta Z^((r))_A) = -1/3 C_"ad" - 1/3 S_H - 2/3 S^((r))$, which is $16 pi^2 ep thin delta Z^((r))_(4 A)$ of @eq:g0_VVVV-Z4, i.e. $Z_4 Z_A = Z_1^2$ with $Z_1 = 1 + delta Z^((r))_(3 A)$, $Z_4 = 1 + delta Z^((r))_(4 A)$ the conventional vertex constants: $2 dot 2/3 gs^2 - 11/3 gs^2 = -7/3 gs^2$ and $2 dot (-1/6) gw^2 - 11/6 gw^2 = -13/6 gw^2$. For QCD with two Dirac quarks this is the textbook $Z_4 - 1 = -(gs^2 \/ 16 pi^2 ep)(C_A\/3 + 4/3 T_F n_f)$ at $xi = 1$ ($C_A = 3$, $T_F = 1\/2$, $n_f = 2$). The $bop_d$ variation of the four-point counterterms together with the two- and three-point ones is compared with the poles of @sec:g1_cVVV_cVVVV.]

#check[*Evanescent part across the gauge classes.* The fermion loops give $z^((r))_"ferm" = 2/3 S^((r))$ times the barred tree monomial in the two-, three- and four-point functions (@sec:g0_VV, @sec:g0_VVV and @eq:g0_VVVV-pole), and the bosonic loops give $d$-dimensional monomials with unsplit indices. With the closed directions continued with unsplit indices (@sec:evanescent), the evanescent pole of the gauge sector is therefore $-lf ep^(-1) z^((r))_"ferm" thin tilde(E)^((r))_F$, with one coefficient for all numbers of legs, plus the two-point operator $E^((r))_A$ of @sec:g0_VV. This is the gauge part of the evanescent counterterm of @sec:st-solution.]

#check[*Literature.* The pole @eq:g0_VVVV-pole, with the gauge and ghost part $1/3 C_"ad"$ times the $d$-dimensional quartic Yang--Mills term, the Higgs part $-1/12 Tr_H (F F)$ ($d$-dimensional) and the fermion part $-1/6 sum_f Tr_f (bar(F) bar(F))$ (barred, both chiralities with equal weight), agrees with the one-loop divergent effective action of @Belusca-Maito:2020ala at $xi = 1$, translated to our conventions.]

#general[At one loop and zero external momenta, $Top(0)$ only sets all external momenta and masses to zero, and the barred structure of fermion loops versus the unsplit structure of bosonic loops follows from @eq:phys-prop alone. From two loops on, $Rbar$ inserts the one-loop counterterms of @sec:g0_VV, @sec:g0_VVV and this subsection into the subgraphs, including the evanescent monomial $tilde(E)^((r))_F |_(A A A A)$ of @eq:g0_VVVV-Z4; contracted with a loop momentum it gives $eta\/ep^2$ terms, which are kept as described in @sec:m-beyond.]
