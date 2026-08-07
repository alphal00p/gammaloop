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
physical `lambda`.  Seven logarithmic samples on each signed branch span
`10^-6 <= |lambda| <= 10^-2`; the singular midpoint is deliberately skipped.

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

`correlated_soft_threshold_cure.pdf` starts with the three kinematic distances
and a threshold-off/legacy/IR-safe comparison.  Further pages show the
weighted and unmultiplied single and iterated components of cut `(2,4,12)`,
followed by occurrence-resolved multiplier contexts.  All pages are rendered
by `assets/plot_approach_result.py`; the local script only computes the
graph-specific diagnostic series and merges plotter-produced PDFs.

The distances scale as `|lambda|^1`.  The original, counterterm sum, and cured
total scale as `|lambda|^-1`; the cure is visible in the leading-coefficient
cancellation and exact event/decomposition closure, not as an artificial
change of this common physical power.

As in the finite analysis, direct-import integrated-UV generation is disabled
for GL638 because of the independent Vakint reconstruction limitation.  Local
UV and local/integrated threshold CTs remain enabled.  Generated state, JSON,
and PDF artifacts are ignored by Git.
