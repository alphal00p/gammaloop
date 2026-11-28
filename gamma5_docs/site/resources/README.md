# Resources of the gamma5 documentation site

Light inputs and results behind the site, for download. Every file is at most 2 MB; the directory holds 4.94 MB in total. This file and `data/resources.json` are written by `tools/build_site.py`.

Absolute local paths, user and host names were removed from the copies: paths are relative to the directory the file describes (the run directory for the run outputs, the project directory `one_loop_explicit_calculation/` for the report's working material), and `<scratch>` stands for a temporary working directory. Nothing else was changed, except that the report's build script now locates the report next to itself.

| directory | content | origin |
|---|---|---|
| `automation/` | model specification of the automation | `python/hepkit_renorm/hepkit_renorm/models/` of this repository |
| `run_sm_one_family/` | outputs of the end-to-end run of `hepkit_renorm` on that specification: tree-level and renormalised HepKIT models, order store, order data, run report, comparison with the report | the run directory of the automation |
| `final_counterterms.json` | the report's final numbers | `SCRATCHPAD/phase4/` of the project |
| `report_headings.txt` | first page of every labelled heading of the report | the report build |
| `report_classes/<key>/` | final result of each class of Green's functions, `<key>` as in the report's class subsections | `SCRATCHPAD/classes/<key>/` |
| `report_checks/` | assembly of the singular counterterm, pole identity, the two solutions of the restoring condition | `SCRATCHPAD/phase4/` |
| `report_model/` | hand-written HepKIT JSON model of the report and its enumeration summary | `SCRATCHPAD/model/` |
| `toolkit/` | interface and validation notes of the stand-alone toolkit `bmhvkit` | `SCRATCHPAD/toolkit/` |
| `report_src/` | Typst sources of the report, its bibliography and its build script | `report/`, `SCRATCHPAD/scripts/build_report.py` |

The report itself is `../report/one_loop_gamma5.pdf`. Third-party documents (papers, preprints) are not included.

## Files

### Report

- `../report/one_loop_gamma5.pdf` (8.21 MB): The report “One-loop renormalisation of the Standard Model in the unbroken phase with a non-anticommuting gamma5” (PDF, 186 pages).
- `report_headings.txt` (3.4 kB): First page of every labelled heading of the report, used for the page references of this site.

### Automation: model specification and run outputs

- `automation/sm_one_family.toml` (1.6 kB): Input of the automation: gauge group, chiral multiplets, scalar doublet, classical four-dimensional Lagrangian, gauge fixing and regularisation of the one-family Standard Model.
- `run_sm_one_family/comparison.json` (3.2 kB): Comparison of the run with the report's final numbers, quantity by quantity.
- `run_sm_one_family/run_report.json` (31.2 kB): Every step, check and result of the run: graph counts per class, renormalisation constants, beta functions, breaking and restoring coefficients.
- `run_sm_one_family/run_report.txt` (9.5 kB): Readable summary of the run report: steps with timings, the fifteen checks and the results of order one.
- `run_sm_one_family/sm_one_family_order_data.json` (478.7 kB): Matched loop data of order one as exact local functionals: the pole coefficients and the finite breaking.
- `run_sm_one_family/sm_one_family_renormalised.json` (1.24 MB): The renormalised model: classical action, breaking, singular and finite restoring counterterms and their breaking insertions, as a HepKIT model.
- `run_sm_one_family/sm_one_family_store.json` (1.35 MB): Exact order store of the run (all functionals of order one), the input of the next loop order.
- `run_sm_one_family/sm_one_family_tree.json` (236.3 kB): Regularised classical action with ghosts, sources and the breaking insertion, in the Feynman gauge, as a HepKIT model.

### Final numbers of the report

- `final_counterterms.json` (48.1 kB): The report's final numbers: singular counterterm, finite breaking, restoring coefficients and literature comparison.

### Per-class results of the report

- `report_classes/g0_FF/final_result.json` (12.8 kB): Final numbers of the class g0_FF, Fermion two-point functions Ψ̄Ψ (Section 4.1.4): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g0_FFS/final_result.json` (8.6 kB): Final numbers of the class g0_FFS, Yukawa vertices Ψ̄ΨH (Section 4.1.9): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g0_FFV/final_result.json` (9.4 kB): Final numbers of the class g0_FFV, Fermion–gauge vertices Ψ̄ΨA (Section 4.1.5): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g0_SS/final_result.json` (7.2 kB): Final numbers of the class g0_SS, Higgs two-point function H†H (Section 4.1.6): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g0_SSSS/final_result.json` (7.3 kB): Final numbers of the class g0_SSSS, Higgs four-point function (H†H)² (Section 4.1.8): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g0_SSV_SSVV/final_result.json` (12.6 kB): Final numbers of the class g0_SSV_SSVV, Higgs–gauge functions H†HA and H†HAA (Section 4.1.7): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g0_VV/final_result.json` (9.4 kB): Final numbers of the class g0_VV, Gauge two-point functions AA (Section 4.1.1): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g0_VVV/final_result.json` (9.3 kB): Final numbers of the class g0_VVV, Gauge three-point functions AAA (Section 4.1.2): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g0_VVVV/final_result.json` (9.7 kB): Final numbers of the class g0_VVVV, Gauge four-point functions AAAA (Section 4.1.3): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g0_ghost_sources/final_result.json` (11.7 kB): Final numbers of the class g0_ghost_sources, Ghost and BRST-source functions (Section 4.1.10): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g1_cFF/final_result.json` (19.9 kB): Final numbers of the class g1_cFF, Breaking of ωΨ̄Ψ (Section 4.2.4): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g1_cFFS/final_result.json` (10.2 kB): Final numbers of the class g1_cFFS, Breaking of ωΨ̄ΨH (Section 4.2.7): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g1_cFFV/final_result.json` (7.0 kB): Final numbers of the class g1_cFFV, Breaking of ωΨ̄ΨA (Section 4.2.5): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g1_cSS_cSSV/final_result.json` (9.1 kB): Final numbers of the class g1_cSS_cSSV, Breaking of ωH†H, ωH†HA and higher (Section 4.2.6): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g1_cV/final_result.json` (14.2 kB): Final numbers of the class g1_cV, Breaking of ωA (Section 4.2.1): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g1_cVV/final_result.json` (11.9 kB): Final numbers of the class g1_cVV, Breaking of ωAA and the anomaly (Section 4.2.2): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g1_cVVV_cVVVV/final_result.json` (11.8 kB): Final numbers of the class g1_cVVV_cVVVV, Breaking of ωAAA and ωAAAA (Section 4.2.3): graphs, loop functions, local operators, checks and caveats of the two independent calculations.
- `report_classes/g1_ccRF/final_result.json` (9.8 kB): Final numbers of the class g1_ccRF, Breaking with BRST sources: ωωR̄Ψ (Section 4.2.8): graphs, loop functions, local operators, checks and caveats of the two independent calculations.

### Assembly and checks of the report

- `report_checks/ct_assembly.json` (34.4 kB): Over-determined assembly of the singular counterterm from the eighteen classes: monomial coefficients, renormalisation constants, beta functions and evanescent part.
- `report_checks/pole_consistency.json` (52.9 kB): Check of the pole identity for the Green's functions with one breaking insertion, leg assignment by leg assignment.
- `report_checks/x_solve_engine.json` (70.1 kB): Solution of the restoring condition with the plane-wave engine: the finite breaking by class and the restoring coefficients.
- `report_checks/x_solve_hand.json` (8.9 kB): Independent hand solution of the restoring condition under the normalisation conditions.

### Hand-written HepKIT model of the report

- `report_model/ENUMERATION.md` (187.1 kB): Counts of the one-loop one-particle-irreducible graphs per class and their comparison with hand counts.
- `report_model/README.md` (8.8 kB): Encoding of the hand-written model and description of the enumeration scripts (which are not copied here).
- `report_model/sm_unbroken_bmhv.json` (83.7 kB): Hand-written HepKIT JSON model of the regularised unbroken Standard Model with sterile partners, used for the report's graph enumeration.
- `report_model/sm_unbroken_bmhv_native.json` (83.7 kB): Variant of the hand-written model with anticommuting (native) ghosts.

### Stand-alone toolkit

- `toolkit/README.md` (15.9 kB): Interface of the stand-alone BMHV and one-loop toolkit bmhvkit (the code itself is not copied here).
- `toolkit/VALIDATION.md` (9.4 kB): Independent validation of the stand-alone BMHV toolkit bmhvkit: what was checked and how.

### Typst sources of the report

- `report_src/build_report.py` (9.1 kB): Assembles the class section and compiles the report with Typst.
- `report_src/diagrams.typ` (5.5 kB): Typst (CeTZ) macros that draw the report's Feynman graphs.
- `report_src/macros.typ` (5.5 kB): Typst macros of the report's notation.
- `report_src/one_loop_gamma5.typ` (3.2 kB): Title, abstract and the list of included sections.
- `report_src/refs.bib` (6.1 kB): Bibliography of the report (BibTeX).
- `report_src/sections/01_introduction.typ` (6.5 kB): Typst source of the report: Introduction.
- `report_src/sections/02_method.typ` (62.1 kB): Typst source of the report: Systematics of the method.
- `report_src/sections/03_setup.typ` (39.4 kB): Typst source of the report: Lagrangian, notation and Feynman rules.
- `report_src/sections/03b_operator_bases.typ` (46.3 kB): Typst source of the report: Counterterm operators and the action of.
- `report_src/sections/04_classes.typ` (0.9 kB): Typst source of the report: Explicit one-loop calculations.
- `report_src/sections/04a_classes_intro.typ` (6.2 kB): Typst source of the report.
- `report_src/sections/04b_g0_intro.typ` (2.1 kB): Typst source of the report.
- `report_src/sections/04b_g1_intro.typ` (3.3 kB): Typst source of the report.
- `report_src/sections/04z_st_solution.typ` (59.8 kB): Typst source of the report: Solution of the Slavnov--Taylor conditions.
- `report_src/sections/05_conclusion.typ` (22.2 kB): Typst source of the report: Conclusion.
- `report_src/sections/05a_ct_feynman_rules.typ` (20.8 kB): Typst source of the report: Feynman rules of the counterterms.
- `report_src/sections/05b_broken_phase.typ` (22.5 kB): Typst source of the report: Subtleties in going to the broken phase.
- `report_src/sections/05c_broken_lagrangian.typ` (47.0 kB): Typst source of the report: The renormalised action in the broken phase.
- `report_src/sections/05r_recipe.typ` (3.6 kB): Typst source of the report: The recipe in brief.
- `report_src/sections/06_automation.typ` (37.7 kB): Typst source of the report: Towards automation.
- `report_src/sections/07_supporting_material.typ` (8.6 kB): Typst source of the report: Supporting material.
- `report_src/sections/classes/g0_FF.typ` (13.3 kB): Typst source of the report: Fermion two-point functions.
- `report_src/sections/classes/g0_FFS.typ` (18.2 kB): Typst source of the report: Yukawa vertices.
- `report_src/sections/classes/g0_FFV.typ` (17.3 kB): Typst source of the report: Fermion--gauge vertices.
- `report_src/sections/classes/g0_SS.typ` (14.8 kB): Typst source of the report: Higgs two-point function.
- `report_src/sections/classes/g0_SSSS.typ` (20.2 kB): Typst source of the report: Higgs four-point function.
- `report_src/sections/classes/g0_SSV_SSVV.typ` (29.7 kB): Typst source of the report: Higgs--gauge functions and.
- `report_src/sections/classes/g0_VV.typ` (19.0 kB): Typst source of the report: Gauge two-point functions.
- `report_src/sections/classes/g0_VVV.typ` (25.2 kB): Typst source of the report: Gauge three-point functions.
- `report_src/sections/classes/g0_VVVV.typ` (29.0 kB): Typst source of the report: Gauge four-point functions.
- `report_src/sections/classes/g0_ghost_sources.typ` (22.4 kB): Typst source of the report: Ghost and BRST-source functions.
- `report_src/sections/classes/g1_cFF.typ` (20.1 kB): Typst source of the report: Breaking of.
- `report_src/sections/classes/g1_cFFS.typ` (21.3 kB): Typst source of the report: Breaking of.
- `report_src/sections/classes/g1_cFFV.typ` (21.0 kB): Typst source of the report: Breaking of.
- `report_src/sections/classes/g1_cSS_cSSV.typ` (22.4 kB): Typst source of the report: Breaking of , and higher.
- `report_src/sections/classes/g1_cV.typ` (14.0 kB): Typst source of the report: Breaking of.
- `report_src/sections/classes/g1_cVV.typ` (16.8 kB): Typst source of the report: Breaking of and the anomaly.
- `report_src/sections/classes/g1_cVVV_cVVVV.typ` (23.7 kB): Typst source of the report: Breaking of and.
- `report_src/sections/classes/g1_ccRF.typ` (22.5 kB): Typst source of the report: Breaking with BRST sources.
