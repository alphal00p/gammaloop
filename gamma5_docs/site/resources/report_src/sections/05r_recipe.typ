#import "../macros.typ": *

// local shorthands (this file only)
#let Dl(l) = $Dbrk^([#l])$
#let STq = $cal(S)_"q"$
#let Kf = $Kop^((eta))$

#block(breakable: false)[
== The recipe in brief <sec:recipe>

The ten steps of @sec:m-algorithm at one loop, with the place where each was carried out:
+ *Model:* propagators and vertices of $S_0$, the $Dbrk_1$ vertex @eq:V-delta, source legs (@sec:feynman-rules).
+ *Ghost-number-0 graphs:* ten classes with $sum dim <= 4$ (@sec:g0).
+ *$Rop$-operation:* at one loop the forest is empty, $Z(G) = -Kop Top(dg(G)) G$ (@eq:m-Z1).
+ *Numerator algebra:* BMHV splitting, traces with $g5$, tensor reduction, the tag $eta$ (@sec:bmhv, @sec:integrals).
+ *Pole matching:* $delta S^((1))$ in the basis $Q_alpha$ plus evanescent operators (@sec:st-dS, @eq:c-dS-inv, @eq:c-dS-ev).
+ *Breaking vertices:* $Dl(1) = Dbrk_1 + bop_d delta S^((1))$ (@eq:m-Delta-l), with $bop_d delta S^((1))$ of @eq:c-bd-dS.
+ *Insertion graphs:* eight ghost-number-1 classes with one $Dbrk_1$ vertex (@sec:g1).
+ *$X_1$:* poles cancel class by class (@tab:sts-poles); finite part $N_1$ projected with $Pphys$ (@eq:m-X1, @tab:sts-X).
+ *Solve $bop Delta S^((1))_"ST" = -X_1$* under (N1)--(N4) (@sec:st-solve, @eq:c-DS).
+ *Update the model:* $delta S^((1))$ and $Delta S^((1))_"ST"$ in their $d$-dimensional forms (@sec:Sren1).
]

*At two loops* (@sec:m-beyond) step 3 has non-empty forests, and child counterterms keep their evanescent structures and hatted momenta. $X_2$ is the finite part of the sum of four classes: (a) two-loop graphs with one $Dbrk_1$ vertex, with the ghost-number-0 subdivergences subtracted by $Rbar$; (b) one-loop graphs with $Dbrk_1$ and one vertex of $Delta S^((1))_"ST"$; (c) one-loop graphs with one vertex of $bop_d delta S^((1)) + bop_d Delta S^((1))_"ST"$ (@eq:c-bd-dS, @eq:c-bd-DS); (d) the tree vertices $bop_d delta S^((2)) + STq(C^((1)))$, with $STq$ the source--field part of $ST$ (@eq:m-S-exact). Two bookkeepings of this sum are possible (@tab:m-bookkeeping). In bookkeeping A the classes are taken as listed: the one-loop subgraphs of (a) that contain $Dbrk_1$ are not subtracted, their poles cancel against (c), and $X_2$ requires the complete finite parts of all classes, not only UV-expanded integrals. In bookkeeping B #cite(<Lazopoulos:2026two>) these subgraphs are subtracted with $Kf Top(dg(gamma))$, the pole part at fixed $eta$ with $eta$ kept until the last integration, and class (c) is replaced by the purely evanescent vertex $hat(Dbrk)^((1))$ of @eq:m-LP. Per class: integrate at fixed $eta$, take $Kf Top(dg(G)) Rbar_eta G$ per graph, sum, apply $Pphys$, set $eta = ep$; $X_2$ is the $ep^0$ term, i.e. the $eta^r ep^(-r)$ terms, plus the finite part of (d) (@eq:m-X2B). The tag follows one rule in graphs and subgraphs: $d$-dimensional tensor reduction, then $-2 eta$ per closed hatted cycle (@sec:m-K). Bookkeeping B needs UV-expanded integrals only, and @sec:m-beyond shows that it gives the whole finite breaking (checked on a toy two-loop graph). At order $l$ every proper subgraph is subtracted at fixed $eta$ and the remainders $hat(Dbrk)^((L))$, $E^((L))$ of @eq:m-LPl are inserted; this extension is argued, not tested. In both, all one-loop counterterms enter in the $d$-dimensional forms of @sec:Sren1, the evanescent parts of @eq:c-dS-ev included. The singular counterterm can contain non-evanescent non-invariant poles. The source counterterms enter quadratically through $STq$ (example in @sec:m-beyond). Since $gh^B$ can then sit on an internal line, the existence condition of (N4) must be checked again, with (N4$'$) as the fallback.
