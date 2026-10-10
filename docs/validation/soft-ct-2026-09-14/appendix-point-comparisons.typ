= Appendix B1: matched complete-point comparisons

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<appendix-b1-matched-complete-point-comparisons>
Existing saved full/bare states were loaded read-only and evaluated with `inspect --use_arb_prec` (ArbPrec=1000bits), graph0, summed orientations. The common canonical point is `(0.11,-0.29,0.37,0.59,-0.43,0.71)`. `first`/`second` scale the corresponding three momentum coordinates; `both` scales both loops. These generic canonical rays are a separate diagnostic, not the original seeded candidate-LMB rays.

Values below are serialized f64 imaginary parts after high-precision evaluation. Tiny subnormal real parts are negligible here. CLI raw inspect returns no stability-accuracy metadata. The relative difference is `|projected − explicit| / max(|projected|, |explicit|)` for the full complex values; it is not a reported accuracy bound. Opposite-sign values can therefore give a difference greater than one.

#figure(
  align(center)[#table(
    columns: 5,
    align: (auto,right,right,right,right,),
    table.header([Point], [Direct3D full Im], [Projected full Im], [Relative full difference], [Relative bare difference],),
    table.hline(),
    [both-1e+04], [-4.705346124908e-27], [-4.705346124908e-27], [5.15378e-14], [0],
    [both-1e+08], [-4.704300374451e-59], [-4.704276086773e-59], [5.16287e-06], [0],
    [both-1e+12], [-4.704300374440e-91], [2.424063421950e-88], [1.00194], [0],
    [first-1e+04], [2.552466281142e-18], [2.552466281142e-18], [7.92044e-300], [0],
    [first-1e+08], [2.551141457117e-38], [2.551141457117e-38], [1.24821e-14], [0],
    [first-1e+12], [2.551141472835e-58], [2.551141473155e-58], [1.252e-10], [0],
    [generic], [-5.897721056547e-10], [-5.897721056547e-10], [1.32487e-304], [0],
    [second-1e+04], [1.097630184922e-21], [1.097630184942e-21], [1.81134e-11], [0],
    [second-1e+08], [-5.438850277621e-38], [-5.418968319275e-38], [0.00365554], [0],
    [second-1e+12], [-5.440492108103e-54], [1.988190394117e-48], [1], [0],
  )]
  , kind: table
  )

At the generic point full values agree to displayed precision. All ten bare comparisons agree exactly after f64 serialization. Hard-scaled full comparisons increasingly differ: at scale10^12, the second-loop and both-loop paths disagree in sign and leading size. This corroborates the failed projected UV profiles while locating the discrepancy in the counterterm-dependent result. It is not an exact reconstruction certificate and does not distinguish coefficient rounding, source reconstruction or later evaluation by itself. No numerator expansion was used.

All four corrected commands exited0. Initial comma-delimited inline `--point` commands were rejected by the CLI minimum-argument parser before evaluation; rerunning with whitespace-separated coordinates succeeded. The rejected attempts are retained in command metadata and are not additional test failures. Exact values are in JSON;; raw per-point JSON and command logs remain under `/tmp/soft-ct-validation-2026-09-14/appendix-points-v2-*`.
