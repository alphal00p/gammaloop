// Notation macros for the report. Every section file starts with:  #import "../macros.typ": *
// (class files in sections/classes/ use  #import "../../macros.typ": * ).
// Keep formulas explicit: prefer writing indices out over inventing new shorthands.

// --- Dirac algebra -------------------------------------------------------------------------------
#let sl(x) = math.cancel(x, angle: 30deg, stroke: 0.45pt)      // Feynman slash:  sl(k) = gamma^mu k_mu  (d-dimensional)
#let gbar = $overline(gamma)$                                   // four-dimensional part of gamma^mu
#let ghat = $hat(gamma)$                                        // (d-4)-dimensional part of gamma^mu
#let gmb = $overline(g)$                                        // four-dimensional metric
#let gmh = $hat(g)$                                             // (d-4)-dimensional metric
#let bar(x) = $overline(#x)$                                    // four-dimensional part of a vector / index
#let slb(x) = $overline(#math.cancel(x, angle: 30deg, stroke: 0.45pt))$   // slashed four-dimensional part: bar above the slashed letter
#let slh(x) = $hat(#math.cancel(x, angle: 30deg, stroke: 0.45pt))$        // slashed evanescent part: hat above the slashed letter
#let g5 = $gamma_5$
#let PL = $P_upright(L)$
#let PR = $P_upright(R)$
#let Pc = $P_c$                                                 // projector of the physical chirality of a multiplet
#let Pcb = $P_(overline(c))$                                    // the opposite projector
#let tr = math.op("tr")                                         // Dirac trace
#let Tr = math.op("Tr")                                         // trace over gauge / multiplet indices

// --- gauge structure -----------------------------------------------------------------------------
#let Tg = $cal(T)$                                              // coupling-weighted generator  cal(T)^a = g_a T^a
#let fg = $frak(f)$                                             // coupling-weighted structure constant
#let Am = $cal(A)$                                              // matrix-valued gauge field  cal(A)_mu = cal(T)^a_f A^a_mu  (in the representation of the multiplet at hand)
#let wm = $bold(omega)$                                         // matrix-valued ghost  bold(omega) = cal(T)^a_f omega^a
#let gs = $g_s$
#let gw = $g$
#let gy = $g'$

// --- fields, ghosts, sources ---------------------------------------------------------------------
#let gh = $omega$                                               // ghost
#let agh = $overline(omega)$                                    // antighost
#let srcA = $rho$                                               // BRST source of the gauge field (gamma_a^mu in the starter pack)
#let srcC = $zeta$                                              // BRST source of the ghost
#let srcF = $R$                                                 // BRST source of a fermion
#let srcH = $sigma$                                             // BRST source of the Higgs doublet

// --- operators -----------------------------------------------------------------------------------
#let sop = $bold(s)$                                            // BRST operator
#let bop = $bold(b)$                                            // linearised Slavnov-Taylor operator of the tree-level action
#let ST = $cal(S)$                                              // Slavnov-Taylor operator
#let Kop = $K$                                                  // pole-part operator (MS-bar)
#let Rop = $cal(R)$                                             // R-operation (calligraphic: R_f is the fermion BRST source)
#let Rbar = $overline(cal(R))$                                  // R-bar: subtraction of all sub-divergences
#let Top(n) = $T^((#n))$                                        // Taylor (UV-expansion) operator of degree n
#let Pphys = $cal(P)$                                           // projector on the physical (four-dimensional) part
#let Evan = $cal(E)$                                            // projector on the evanescent part,  E = 1 - P
#let Dbrk = $Delta$                                             // breaking operator  Delta = S(S_ren)

// --- loop integrals ------------------------------------------------------------------------------
#let ep = $epsilon.alt$                                         // regulator, d = 4 - 2 ep; lunate glyph. The curly epsilon is reserved for Levi-Civita symbols
#let intk = $integral_k$                                        // = mu^(2 eps) int d^d k / (2 pi)^d
#let lf = $1/(16 pi^2)$                                         // one-loop factor
#let ilf = $i/(16 pi^2)$
#let MUV = $M$                                                  // auxiliary UV mass of the expansion operator
#let dg = math.op("deg")                                        // superficial degree of divergence of a graph: dg(G)  (omega is the ghost)

// --- typography ----------------------------------------------------------------------------------
#let eqbox(body) = align(center, box(stroke: 0.5pt + luma(120), inset: 7pt, radius: 2pt, body))
#let result(body) = block(width: 100%, stroke: (left: 1.5pt + rgb("#1f4e79")), inset: (left: 9pt, y: 5pt), body)
#let check(body) = block(width: 100%, stroke: (left: 1.5pt + rgb("#2e7d32")), inset: (left: 9pt, y: 5pt), [*Check.* #body])
#let general(body) = block(width: 100%, stroke: (left: 1.5pt + rgb("#8e44ad")), inset: (left: 9pt, y: 5pt), [*Beyond one loop.* #body])
