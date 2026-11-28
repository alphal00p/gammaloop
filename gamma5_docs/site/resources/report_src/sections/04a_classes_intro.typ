#import "../macros.typ": *

A class is a group of 1PI functions with related external legs, for example $A A$ or $gh overline(Psi) Psi A$, taken for all gauge factors and multiplets at once; it is one row of @tab:classes and one subsection. Rows that combine several sets of external legs list them separated by semicolons. The ghost-number-0 classes (@sec:g0) have $sum_j dim Phi_j <= 4$ and give the singular counterterm @eq:m-Z1, $i thin delta S^((1))_(Phi_1 dots Phi_n) = -ilf sum_G N_0 (G) \/ ep$. The ghost-number-1 classes (@sec:g1) carry one vertex $Dbrk_1$, have $sum_j dim Phi_j <= 5$, and give the finite breaking @eq:m-X1, $i thin X_(1; Phi_1 dots Phi_n) = sum_G Pphys [ilf N_1 (G)]_(ep = 0)$, and the cancellation of their poles against $bop_d delta S^((1))$.

*Layout.* An opening paragraph gives the external legs (all momenta incoming), the function computed, $i Gamma^((1))$ or $i (Dbrk_1 dot Gamma)^((1))$, the superficial degree of divergence $dg(G)$ (@sec:integrals) and the gauge factors and multiplets covered. Four run-in titles follow in a fixed order. _Diagrams_: one drawing per topology; a table of field assignments with symmetry factor, sign of a closed fermion or ghost loop, and group factor. _Integrands_: each topology from the Feynman rules of @sec:feynman-rules, with the momentum routing. _$Rop$-operation_: $Top(dg(G))$ (@eq:uv-expansion), the algebra of @sec:bmhv and the reduction @eq:tensor-red-split, ending with $N_0$ and $N_1$ per topology in the form @eq:N0N1. _Result_: the framed result and a table of Standard-Model values. Check blocks close the subsection. Coefficients of tree structures in $N_0$ get names local to their subsection, stated where they are defined. Renormalisation constants ($delta Z_x$, $delta m^2$, $delta lambda$, $delta y_k$) always include the factor $lf 1\/ep$; the number without it is written $16 pi^2 ep thin delta Z_x$.

*Blocks.* A block with a blue rule is the result. At ghost number 0 it states $Gamma^((1))|_"pole" = lf N_0 \/ ep$ and $delta S^((1)) supset -lf ep^(-1) integral dif^d x [dots]$. At ghost number 1 it states $(Dbrk_1 dot Gamma)^((1))|_"pole" = lf N_0 \/ ep$, $X_(1; Phi_1 dots Phi_n) = lf Pphys N_1$ and $X_1 supset lf integral dif^4 x [dots]$ in the basis $N_beta$ (@tab:Nbasis). A block with a green rule, headed *Check.*, tests the result: independence of $MUV$, symmetries, known limits, relations to other classes, independent evaluations. A block with a purple rule, headed *Beyond one loop.*, states the general-order construction of which the step is the one-loop case; it is not used at one loop.

#figure(
  text(size: 8.6pt, table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (left, left, center, right, left, left),
    inset: (x: 3.2pt, y: 2.6pt),
    stroke: (x, y) => if y == 0 or y == 10 { (bottom: 0.6pt) } else { none },
    table.header([class], [external legs], [$dg(G)$], [diagrams], [pole], [$X_1$]),
    [@sec:g0_VV], [$A A$], [$2$], [$20$], [yes, evan.], [--],
    [@sec:g0_VVV], [$A A A$], [$1$], [$54$], [yes, evan.], [--],
    [@sec:g0_VVVV], [$A A A A$], [$0$], [$237$], [yes, evan.], [--],
    [@sec:g0_FF], [$overline(Psi) Psi$], [$1$], [$16$], [yes], [--],
    [@sec:g0_FFV], [$overline(Psi) Psi A$], [$0$], [$46$], [yes], [--],
    [@sec:g0_SS], [$H^dagger H$], [$2$], [$8$], [yes, evan.], [--],
    [@sec:g0_SSV_SSVV], [$H^dagger H A$; $H^dagger H A A$], [$1$; $0$], [$146$], [yes, evan.], [--],
    [@sec:g0_SSSS], [$H^dagger H H^dagger H$], [$0$], [$63$], [yes], [--],
    [@sec:g0_FFS], [$overline(Psi) Psi H$, $overline(Psi) Psi H^dagger$], [$0$], [$32$], [yes], [--],
    [@sec:g0_ghost_sources], [$srcA gh$; $srcA A gh$, $srcC gh gh$, $overline(srcF) gh Psi$, $srcH^dagger gh H$, conj.], [$1$; $0$], [$20$], [yes], [--],
    [@sec:g1_cV], [$gh A$], [$3$], [$10$], [evan.], [yes],
    [@sec:g1_cVV], [$gh A A$], [$2$], [$40$], [evan.], [yes],
    [@sec:g1_cVVV_cVVVV], [$gh A A A$; $gh A A A A$], [$1$; $0$], [$1320$], [$0$], [$gh A A A$ only],
    [@sec:g1_cFF], [$gh overline(Psi) Psi$], [$1$], [$42$], [evan.], [yes],
    [@sec:g1_cFFV], [$gh overline(Psi) Psi A$], [$0$], [$220$], [$0$], [$gh^G$, $gh^W$],
    [@sec:g1_cSS_cSSV], [$gh H^dagger H$; $+A$; $+A A$; $gh (H^dagger H)^2$], [$2$; $1$; $0$; $0$], [$525$], [evan. ($gh H^dagger H$)], [$gh H^dagger H$, $gh H^dagger H A$],
    [@sec:g1_cFFS], [$gh overline(Psi) Psi H$, $gh overline(Psi) Psi H^dagger$], [$0$], [$122$], [$0$], [$gh^W$ only],
    [@sec:g1_ccRF], [$gh gh overline(srcF) Psi$, conj.; source functions], [$0$], [$54$], [$0$], [$gh^G$, $gh^W$],
  )),
  caption: [The eighteen classes. "$dg(G)$": superficial degree of divergence. "Diagrams": one-loop 1PI diagrams at multiplet level, generated with the diagram generator HepKIT (@sec:automation), summed over the external assignments allowed by rigid invariance. A multiplet-level diagram has fixed multiplets and gauge factors on its external legs and its internal lines labelled by multiplet and orientation ($W$ for $W^1, W^2, W^3$, $Q$ for both $"SU"(2)_L$ components); identical external legs are distinguishable, and diagrams with identically vanishing couplings are not generated. The class subsections count fewer graphs, because there one graph has generic gauge indices and a group factor summed over the gauge factors of its external and internal lines (e.g. $29$ graphs in @tab:g0_SSSS-diagrams against $63$ here). The 11251 component diagrams of @sec:auto-generation are the same diagrams split further by the $"SU"(2)_L$ components of the external and internal lines; each multiplet-level diagram is the union of its component diagrams, and the column sums to $2975$. "Pole": "evan." marks an evanescent part (at ghost number 0: barred fermion-loop tensors next to $d$-dimensional bosonic ones, or hatted external momenta). "$X_1$": finite breaking in the Standard Model; at ghost number 0, $N_1 = 0$ (@sec:m-K). "Source functions": the five functions of @sec:g1 without one-loop diagram. In @sec:g1_cFF, @sec:g1_cFFV, @sec:g1_cFFS and @sec:g1_ccRF part of the diagrams have the ghost of $Dbrk_1$ internal, ending on a ghost--gauge vertex that carries an external ghost (@sec:g1).],
) <tab:classes>
