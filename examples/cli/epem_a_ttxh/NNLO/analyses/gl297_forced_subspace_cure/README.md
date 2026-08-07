# GL297 forced-subspace cure

This analysis regenerates one selected GL297 orientation twice: once from the
directive-free `GL297_legacy.dot`, and once from the dressed graph in
`../../graphs/GL297.dot`.  Every one of the nine process-valid cuts is retained,
so each approach evaluates the complete Local-Unitarity cut sum.  The dressed
graph forces the problematic e12 and e13 threshold counterterms into their
physical one-loop subspaces; the legacy graph leaves GammaLoop's maximal-
subspace construction unchanged.

The two approach blocks use the physical parameter `lambda` directly over
`10^-6 <= |lambda| <= 10^-2`, on both signed branches.  In each direction the
soft momentum and the two relevant E-surface equations vanish with the same
parameter.  Full local and integrated UV generation remains enabled.

From the repository root, build fresh state and approach JSON with:

```bash
MallocNanoZone=1 CARGO_INCREMENTAL=0 target/debug/gammaloop \
  --clean-state -o -l info -L off -p none \
  examples/cli/epem_a_ttxh/NNLO/analyses/gl297_forced_subspace_cure/analysis.toml \
  run all
```

Then create both curated PDFs exclusively through the shared approach plotter:

```bash
python3 examples/cli/epem_a_ttxh/NNLO/analyses/gl297_forced_subspace_cure/plot.py
```

The outputs are `e12_correlated_soft_threshold.pdf` and
`e13_correlated_soft_threshold.pdf`.  Each contains:

1. the kinematic soft/E-surface distances;
2. the threshold-off, legacy-maximal, and forced-one-loop totals; and
3. the selected counterterm's summary and expanded-component pages.

The expected fitted powers of the total are approximately `-1`, `-3`, and
`-1`, respectively.  “Threshold subtraction off” is not called “bare”: it is
still the full graph/orientation and nine-cut LU sum.  Within component pages,
“unmultiplied CT” has the narrower meaning of a counterterm before its user
multiplier.

The generated state, JSON, and PDFs are intentionally ignored by Git.
