# Center-inside check for `axis_x2` im[-] drop

Runtime used for the check:

- `xi = 4 * (p1 + p2)` in the exact-xi helper slots.
- `sliver_width = 30`, threshold `h_function sigma = 1`.
- `SingleParametric` fixed-bin inspect, graph `0`, LMB channel `9`.
- Points checked:
  - before drop: `x2 = 0.9441182287966725`
  - peak: `x2 = 0.9441190487748046`
  - after drop: `x2 = 0.9441185115837143`

The instrumentation evaluates every E-surface in the active overlap group at the
same `rotated_center` used by `solve_rstar()` before solving the local
projection.

Summary:

| point | dumped CT projections | present center evaluations | `inside=false` |
| --- | ---: | ---: | ---: |
| before drop | 86 | 1234 | 0 |
| peak | 86 | 1234 | 0 |
| after drop | 86 | 1234 | 0 |

Across all present surfaces in these three checks, the least-negative center
value is approximately `-2.52621131e+02`; no evaluated center is close to the
outside boundary.

Dominant surfaces around the drop:

| surface dump | edges | `E(center)` | inside |
| --- | --- | ---: | --- |
| `GL05_isr_p2_ct_selected_esurf4_group196_raised4_existing69_center.txt` | `[9, 12, 16]` | `-2.25262113e+03` | true |
| `GL05_isr_p1_ct_selected_esurf34_group192_raised34_existing65_center.txt` | `[9, 10, 15]` | `-2.25262113e+03` | true |
| `GL05_selected_esurf5_group140_raised5_existing49_center.txt` | `[9, 12, 16]` | `-2.25262113e+03` | true |
| `GL05_selected_esurf31_group143_raised31_existing52_center.txt` | `[9, 10, 15]` | `-2.25262113e+03` | true |

Conclusion: the sharp `axis_x2` drop is not explained by the center being
outside any E-surface that uses this center for local threshold projection. The
runtime centers used for these overlap groups satisfy the inside prerequisite
with a large margin.
