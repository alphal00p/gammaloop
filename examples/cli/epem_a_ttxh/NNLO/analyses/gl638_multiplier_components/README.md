# GL638 multiplier/component smoke analysis

This analysis direct-imports the final illustrative `GL638.dot`, forces all six
process-valid cuts, and generates exactly the verified orientation.  The
threshold `(7,8)` has independent `[7]` (`intrinsic_1l`) and `[3,7]`
(`embedded_2l`) variants under all four current associations.  Each cut uses a
pair of complementary squared-E-surface ratios built from actual E-surfaces in
that cut's selected-orientation catalog.  Three cuts use `(8,12,14)` as their
secondary surface.  Cut `(2,4,10,13)` has no other active threshold, so its
only legal secondary is the owning cut equation itself.  All unrelated
thresholds remain implicit legacy defaults.

The ratios are deliberately illustrative.  This folder demonstrates finite
multiplier evaluation, expanded L/I component accounting, event cancellation,
and serialization; it does **not** claim a tuned physical cure for GL638.

Run from the repository root:

```bash
MallocNanoZone=1 CARGO_INCREMENTAL=0 cargo run -p gammaloop-api --bin gammaloop -- \
  --clean-state -o -l info -L off -p none \
  examples/cli/epem_a_ttxh/NNLO/analyses/gl638_multiplier_components/analysis.toml \
  run generate inspect_full_cut_samples approach_safe

python3 examples/cli/epem_a_ttxh/NNLO/analyses/gl638_multiplier_components/plot.py
```

The two coordinate-space `inspect` inputs are exactly the finite full-six-cut
samples used by the Rust runtime/save-load regression.  `approach_safe` is a
momentum-space JSON smoke trajectory around a generic point chosen away from
the simultaneous zero of either ratio denominator.  It is component evidence,
not an IR-scaling acceptance criterion.  `plot.py` shows the event-level
original, summed counterterms, and reported total-weight magnitudes without
exposing every expanded component as a separate line.  The schema-v2 JSON also
retains the reconstructed decomposition; its equality to the reported event
weight is checked exactly in the regression tests rather than drawn as a
duplicate curve.

Integrated UV generation is disabled in this direct-import card because the
existing GL638 Vakint reconstruction fails before threshold generation with
`Could not solve Inconsistent`.  Local UV physics and both local/integrated
threshold components remain enabled; this independent CLI gap is recorded
rather than hidden.
