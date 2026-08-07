# GL638 correlated soft/two-threshold cure

This is the physical double-sided GL638 limit used to validate the IR-safe
directives.  It regenerates directive-free and dressed copies of the graph with
one selected problematic orientation and all six process-valid cuts.  Three
approach results compare threshold subtraction off, legacy maximal subspaces,
and the dressed graph's duplicated subspaces and multipliers.

With `K* = sqrt(500^2 - MT^2) (3,-4,12)/13`, the card evaluates

```text
K0(lambda) = K* + lambda (47,29,-31)
K3(lambda) = K* + lambda (258,-128,162)
```

while `K1` and `K2` remain fixed.  Thus `q13 = K0-K3`, `eta(5,10)`, and
`eta(5,12,13)` approach their common soft/threshold point with the same
physical `lambda`. Fifty logarithmic samples on each signed branch span
`10^-6 <= |lambda| <= 10^-2` (100 evaluated points per result); the singular
midpoint is deliberately skipped.

From the repository root, generate all three real CLI outputs with:

```bash
MallocNanoZone=1 CARGO_INCREMENTAL=0 target/debug/gammaloop \
  --clean-state -o -l info -L off -p none \
  examples/cli/epem_a_ttxh/NNLO/analyses/gl638_correlated_soft_threshold_cure/analysis.toml \
  run all
```

Create the curated PDF with:

```bash
python3 examples/cli/epem_a_ttxh/NNLO/analyses/gl638_correlated_soft_threshold_cure/plot.py
```

`correlated_soft_threshold_cure.pdf` shows, on one signed symmetric-log axis:

1. the three kinematic distances;
2. the threshold-off/legacy/IR-safe totals;
3. the three leading Local-Unitarity cut contributions and their sum;
4. the corresponding per-cut threshold-counterterm contributions and their
   sum; and
5. the shared one-loop `(8,12,14)` threshold components across those cuts;
6. the unmultiplied local/integrated pieces of both duplicated `(7,8)`
   variants; and
7. their unmultiplied `LL`, `LI`, `IL`, and `II` Cartesian products with the
   `(5,10)` right threshold.

All pages are rendered by `assets/plot_approach_result.py`; the local script
only computes the graph-specific diagnostic series and merges plotter-produced
PDFs.

The distances scale as `|lambda|^1`. Each leading cut and per-cut threshold CT
scales as `|lambda|^-3`, while its cross-cut sum scales as `|lambda|^-1`. The
threshold-off, legacy, and IR-safe total powers are approximately `-1`, `-3`,
and `-1`, respectively.

This trajectory approaches the right thresholds `(5,10)` and `(5,12,13)`, not
the duplicated left threshold `(7,8)` in effective coordinates. The latter is
reached only in each counterterm's `star` sample. Consequently, the current
own-surface `eta(star, ...)` ratio gives `intrinsic_1l = 0` and
`embedded_2l = 1` there. Those variants remain useful for structural/runtime
coverage, but their weighted intrinsic contribution is exact-zero skipped or
root-residual noise rather than useful soft-cure evidence. The report therefore
keeps the cure evidence (pages 1--5) separate from the unmultiplied variant
anatomy (pages 6--7), where skipped evaluations appear as gaps. Before
weighting, the intrinsic single pieces scale as `|lambda|^-3` and the embedded
single pieces are approximately constant. The cancellation on pages 2--5 is
supplied by the shared `(8,12,14)` one-loop associations.

As in the finite analysis, direct-import integrated-UV generation is disabled
for GL638 because of the independent Vakint reconstruction limitation.  Local
UV and local/integrated threshold CTs remain enabled.  Generated state, JSON,
and PDF artifacts are ignored by Git.
