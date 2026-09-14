# GL638 screen maxima: replay validation

All three screen-maximum batches reproduce their recorded signed extrema and
pass the unchanged ordinary/Arb numerical criteria. This validates the retained
samples and values. The separate [three-point physical attribution](GL638_MAXIMUM_PHYSICS_ATTRIBUTION.md)
identifies ordinary threshold and moved-star contributions at these samples;
neither result establishes a global maximum or tail bound.

| Batch | Original runs | Distinct Samples | Bit-exact signed extrema | Arb controls | Existing accuracy checks |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baseline seed 11113 | 1 | 3 | 4 | 2 | 51 |
| Six advanced methods, seed 11113 | 6 | 15 | 24 | 12 | 266 |
| All seven methods, seeds 22229/33331 | 14 | 41 | 56 | 14 | 552 |
| Total | 21 | 59 | 84 | 28 | 869 |

The independent audit passes 541 checks. All 336 files in the three frozen
included pilot directories and all 35 generated-state file hashes per batch
remain unchanged. Every ordinary call accepts Double and retains six physical
cuts. The largest ordinary/Arb relative complex-norm errors are 4.81073×10⁻¹²
for a total and 6.28735×10⁻¹³ for a cut contribution. All original norm budgets
pass; 31 near-zero component-relative diagnostics remain visible. This does
not promise separate relative accuracy for negligible components.

The 84 signed real/imaginary extrema identify 59 distinct complete Samples.
Forced-Arb controls cover up to two available unique retained-extremum Samples
per method and batch, ranked by their replayed complex norm. They are not a
global top two over all samples: an unrecorded norm maximum is not recoverable
from component extrema. Original serialized integration checkpoints remain in
the [screen archive](gl638_hosted_joint_gate/mc_screen/artifact_hashes.json).
Replay rows preserve the nested Sample weights, extremum bits, source aliases,
native totals/events and CT decomposition. Returned native events already
contain the outer Sample factor; it is not multiplied again.

The three later scoped Arb traces leave totals, factors and complete six-cut
events exactly unchanged. A separate 356-check audit validates their control
identity and structural component/variant/occurrence joins. Their 1,344/1,400/
1,456 JSONL records retain 34/36/38 native centers. Root-group headers from the
first pass cannot label later centers: the physical reconstruction uses actual
evaluator group, side, local threshold, parent, center and root occurrences.
The physical note states its identity-probe and numerical-reconstruction limits.

The [98-record archive](gl638_hosted_joint_gate/mc_maximum_replay/artifact_hashes.json)
retains all original maximum reports, raw traces, run/resource/log records,
unchanged auditing source, final attribution analysis and build/source links.
Exact gzip copies preserve original hashes; byte-identical prior records share
storage. No generated-process state payload or executable is copied. Calls
used one worker and overlapped confirmation, so elapsed/resources are provenance,
not standalone timing evidence. Source/build5 and physical settings are unchanged.

The [confirmation interruption](GL638_CONFIRMATION_ROOT.md) is separate:
three completed runs are retained without pooling after the ordinary baseline
fails canonical physical root certification at its next seed. Its [endpoint repair](CANONICAL_LU_ENDPOINT_CERTIFICATION.md)
now passes the exact physical replay. A restored confirmation and supported
cross-section estimate remain pending.
