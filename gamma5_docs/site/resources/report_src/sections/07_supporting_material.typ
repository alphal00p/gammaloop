#import "../macros.typ": *

= Supporting material <sec:supporting>

The PDF is accompanied by the project directory `one_loop_explicit_calculation/`; all paths below are relative to it. Each entry gives a directory, its content and one check command, run inside that directory with a Python 3 that provides sympy, numpy and mpmath. Entries marked (H) also need HepKIT, i.e. the Symbolica community package, with the licence variable `SYMBOLICA_LICENSE` set in the shell. The commands followed by a time were rerun on the final result for this appendix and passed; the others were run when the results were produced (the regeneration commands overwrite the stored diagram lists).

#[
#set text(size: 9pt)
#set par(leading: 0.5em, spacing: 0.75em)
#set terms(separator: [ --- ], hanging-indent: 1em, spacing: 0.75em)
#let cmd(c, t: none) = [ _Check:_ #raw(c)#if t != none [ (#t)].]

/ `report/`: Typst sources of this report and the bibliography (`refs.bib`). #cmd("typst compile one_loop_gamma5.typ")
/ `bmhv/`: the gammaloop checkout (branch `bmhv`, commit `055fd0c1d`) at which the code locations of @sec:automation were read.
/ `SCRATCHPAD/toolkit/`: the stand-alone BMHV and one-loop toolkit `bmhvkit` (sympy, numpy). BMHV Dirac algebra, symbolic and with explicit matrices for $d = 4, 6, dots, 12$; $Top(dg(G))$ with auxiliary mass $MUV$, tensor reduction @eq:tensor-red-split, vacuum integrals, $N_0$ and $N_1$; group theory of the Standard Model. `README.md` documents the interface. #cmd("python tests/run_tests.py", t: [31 tests, 20 s])
/ `SCRATCHPAD/st/`: the Slavnov--Taylor plane-wave engine. The Standard Model in explicit components; $sop$ and #bop on random plane waves; the bases $Q_alpha$, $N_beta$, the #bop table, the kernel and the normalisation (@sec:bases). #cmd("python test_dictionary.py", t: [examples of @sec:dictionary, 1 s])
/ `SCRATCHPAD/st/verify/`: an independent re-implementation of the engine (own Grassmann algebra and bases), in $D = 4$ and $D = 6$. #cmd("python check_tree.py", t: [tree-level BRST identities, 1 min])
/ `SCRATCHPAD/classes/<key>/`: one directory per class; `<key>` is the label of the class subsection, e.g. `g0_VV` for @sec:g0_VV. It holds the final numbers of the class (`final_result.json`), the scripts of the two independent calculations and the referee notes. Example, in `g1_cFF/`: #cmd("python check_st_pole.py", t: [pole identity, 9 s])
/ `SCRATCHPAD/phase4/ct_assembly/`: the over-determined assembly of $delta S^((1))$ from the eighteen classes (exact, sympy) and its check with the engine. #cmd("python assemble.py")
/ `SCRATCHPAD/phase4/x_solve_engine/`, `x_solve_hand/`, `st_referee/`: the two independent solutions of $bop Delta S^((1)) = -X_1$, their exact comparison, the final check of $bop Delta S^((1))_"ST" + X_1 = 0$ and the literature comparison (@sec:st-solution). In `x_solve_hand/`: #cmd("python check_matrix_forms.py", t: [1 s])
/ `SCRATCHPAD/phase6/ct_rules/`: checks of the counterterm Feynman rules of @sec:ct-rules: the rules of $Delta S^((1))_"ST"$ against plane-wave functional differentiation with the engine, the evanescent rules in $d = 6$ (`engine6.py`), and the factors of $delta S^((1))$ against the class poles (sympy). #cmd("python check_finite_rules.py", t: [65 comparisons, 5 s])
/ `SCRATCHPAD/phase4/pole_consistency/`: the engine in $D = 6$ (four barred and two hatted components) and the pole identity of @sec:st-poles. #cmd("python extra_checks.py", t: [1 min])
/ `SCRATCHPAD/phase4/broken_lagrangian/`, `broken_lagrangian_verify/`: the two independent expansions of the broken phase (@sec:Sren-broken). In `broken_lagrangian/`: tree action and counterterms after the shift (`tree_and_dS.py`), the finite counterterm after the shift (`finite_st.py`), propagators (`propagators.py`). In `broken_lagrangian_verify/`: the independent expansion (`expand.py`), propagators (`propagators.py`), the massless photon, $Z$--photon and gluon at zero momentum (`massless_check.py`) and its cross-check by the Ward identity (`ward_check.py`). In `broken_lagrangian_verify/`: #cmd("python expand.py", t: [123 checks, 6 s])
/ `SCRATCHPAD/scripts/broken_phase_draft/`: the algebraic checks of @sec:broken: the shifted $H$-dependent monomials of @tab:Qbasis, the tree masses and mixings, the $R_xi$ gauge-fixing function and the ghost mass matrix (`check_shift.py`, sympy) and the fermion propagator with sterile partners and its limits (`check_propagator.py`). #cmd("python check_propagator.py", t: [1 s])
/ `SCRATCHPAD/phase5/e2e_scalar/`: the checks of the restored identity with complete finite parts of @sec:st-checks: the gauge two-point function (`partA_fermion.py`, `partA_bosonic.py`) and the Ward identity with legs $(gh^B, H^dagger, H)$ (`partB_higgs.py`). #cmd("python partA_fermion.py", t: [1 s])
/ `SCRATCHPAD/phase5/e2e_fermion/`: the Ward identity with legs $(gh, overline(Psi)_f, Psi_f)$ with complete finite parts, at three off-shell spacelike kinematic points, for the hypercharge ghost (all five multiplets, $m^2 = 0.59$ and $m^2 = 0$) and for the gluon ghost ($Q$, $U$). The $1\/ep$ poles cancel, the finite violation equals $X_1$ of @sec:g1_cFF, and the terms $Q^(f r)_V$, $Q^(f r)_R$ of $Delta S^((1))_"ST"$ cancel it to $10^(-13)$. #cmd("python abelian_wi.py", t: [3 min])
/ `SCRATCHPAD/phase5/numbers_g1/`: an independent numerical audit, $bop Delta S^((1))_"ST" + X_1 = 0$ with the coefficients re-typed from the printed tables. #cmd("python check_st.py", t: [30 s])
/ `SCRATCHPAD/model/` (H): the HepKIT JSON model `sm_unbroken_bmhv.json` and its variant with true ghost numbers `sm_unbroken_bmhv_native.json`; one DOT file per component diagram, `dot/<key>/<assignment>/<tuple>_D<k>.dot`; counts and hand-count comparisons in `ENUMERATION.md`; the encoding in `README.md` (@sec:auto-generation). #cmd("python tests/t_model_consistency.py", t: [1 s]) _Regenerate:_ `python build_model.py`, `python build_model.py --native-ghosts`, `python enumerate_diagrams.py --jobs=12` (about 20 s), `python merge_and_report.py`.
/ `SCRATCHPAD/scripts/automation_draft/` (H): the reproducers `r1` to `r4` of @sec:automation. #cmd("python -I r1_idenso_bmhv.py")
/ `SCRATCHPAD/notes/`, `SCRATCHPAD/targets/`: one digest `lit_*.md` per paper of the literature; literature results in the conventions of this report; notes on the gammaloop $Rop$-operation, on HepKIT, and the design of a BMHV mode for Idenso (M5).
/ `../references/arxiv/`: the arXiv sources of the preprints of the bibliography, one directory per identifier.
]

*Independent verifications of the final result.* @tab:verifications lists what was checked and in how many independent ways, i.e. by a separate code or a separate hand calculation.

#[
#show figure: set block(breakable: true)
#figure(
  placement: none,
  text(size: 8.5pt, table(
    columns: (auto, 1fr, auto),
    align: (left, left, center),
    inset: (x: 4pt, y: 3pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.6pt) } else { (bottom: 0.2pt + luma(200)) },
    table.header([quantity], [how it was checked], [ways]),
    [$N_0$, $N_1$, $X_1$ per class], [two calculations per class, compared term by term, and a referee comparison], [2],
    [$delta S^((1))$], [over-determined assembly: 25 monomials fix the coefficients, the other 72 are reproduced exactly], [1],
    [pole identity (@sec:st-poles)], [per class, and in $D = 6$ with explicit components on all 146 leg assignments], [2],
    [$Delta S^((1))_"ST"$], [hand solution and engine solution, equal on all 97 coefficients], [2],
    [$bop Delta S^((1))_"ST" + X_1 = 0$], [engine on 39 ghost-number-1 leg classes; re-implementation with the coefficients re-typed from the printed tables, at two further coupling points], [2],
    [$bop X_1 = 0$], [engine on 34 leg classes of ghost number 2], [1],
    [literature], [the representatives of #cite(<Belusca-Maito:2020ala>) and #cite(<Lazopoulos:2026two>) differ from $Delta S^((1))_"ST"$ by #bop;-closed terms (@sec:st-lit)], [2],
    [complete finite parts], [gauge two-point function (all gauge factors), the identity with legs $(gh^B, H^dagger, H)$ (@sec:st-checks), and the identity with legs $(gh, overline(Psi)_f, Psi_f)$ for the hypercharge and gluon ghosts], [3],
    [broken phase], [independent expansion (123 checks); massless photon, photon--$Z$ mixing and gluon at zero momentum, cross-checked by a Ward identity], [2],
    [all numbers], [two independent numerical audits, ghost number 0 and ghost number 1; no wrong number found], [2],
  )),
  caption: [Independent verifications of the final result.],
) <tab:verifications>
]
