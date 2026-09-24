= GL638 LU-ray physical replay
<gl638-lu-ray-physical-replay>
The retained LU-ray implementation passes a bounded physical regression at two stored GL638 points: the ordinary stack and forced Arb at hard `hz_03`, plus forced Arb at the `soft22` control. Each call evaluates the actual #strong[936-orientation sum];, all six physical cuts, and the original threshold, local 3D UV and integrated UV terms. All totals and six event weights are finite, all event identities match, and evaluation metadata report no NaN.

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,right,),
    table.header([Case], [Old → new precision], [Relative total change], [Largest relative event change],),
    table.hline(),
    [Hard `hz_03`, stack], [Double → Double], [2.1690e-16], [1.9455e-16],
    [Hard `hz_03`, forced], [Arb → Arb], [1.3995e-294], [3.0385e-294],
    [Double-soft `soft22`, CT on], [Arb → Arb], [1.3671e-281], [2.4737e-296],
  )]
  , kind: table
  )

Relative change is `norm(new−old)/max(norm(new),norm(old))` for the complex value. These are raw momentum-space results: no sampling Jacobian or Monte Carlo weight is applied. The #link("GL638_LU_RAY_REPLAY.json")[machine artifact] preserves both native components, every event identity/weight, precision metadata, original binary64 tokens and exact reproduction files. The hard point has threshold-normal radius `R=0.0002 GeV` as defined in the #link("GL638_X2_PHYSICAL_PILOT.typ")[earlier physical pilot];; both points use the existing precise source\'s exact promotion of the original binary64 values.

The comparison baseline is #link("GL638_X4_ALPHA_REPLAY.typ")[source 9d2bdf483];. Additional commits intervene, so this is a physical regression comparison rather than isolated attribution to the LU-ray change. The hard-point Double-versus-own-Arb distance changes from `5.1540137307e-10` to `5.1540158267e-10`, a slight increase of about `2.10e-16`. No accuracy or variance improvement is claimed. These six simple cut groups all have `max_occurrence=1`; higher raised eta derivatives are covered by the generated raised-cut fixtures, not this GL638 state.

The candidate was linked from base `a8b573f99a4ed2eb0f43982ef8745e2ffdecb262` plus exact tracked source diff SHA256 `6154702655bacc62d4b2faa9ff453dea79686301a8f9b9489f47467c5244ff35`. The executable is SHA256 `db61d1140bc3ade60dfd0c5445db5c64f213607099921e59fd42d4b93bd0dff9`; the newly built API fingerprint `a196f440a26dc6bc` selects core `f3a6c38e9178a92c`. Actual library/source hashes and compiler argv are embedded. After linking, the only source change introduced a private name for the identical ray-energy tuple and moved its comments. Expanding that alias reproduces the compiled source excluding comments/whitespace; the exact diff and independent check are retained.

The one-load batch completed with exit0 in #strong[551.459 seconds];, including #strong[13.163 seconds] for the three evaluation calls. These unoptimized correctness timings do not assess the sampling-cost budget. Exact complete production-key and orientation-vector sets and the original six cut equations matched before evaluation. All #strong[35 saved-state hashes and the complete file set] were unchanged afterward. Events stayed enabled, and native decimal strings preserved values outside the binary64 reporting range.

To reproduce, extract the embedded driver, linker, scripts, card and requests into a scratch directory; relocate their paths consistently and supply the identified checkpoint and frozen baseline results. Link the recorded source snapshot through its verified current API dependency fingerprint, run the single-load batch, then run the decimal comparison and post-state checks. No integration or state regeneration is part of this regression.
