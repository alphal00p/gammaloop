# GL638 cut-(2,4,12) cure and multiplier/component analysis

This analysis direct-imports the final directive-dressed `GL638.dot`, forces all six
process-valid cuts, and generates exactly the verified orientation.  The
threshold `(7,8)` has independent `[7]` (`intrinsic_1l`) and `[3,7]`
(`embedded_2l`) variants under all four current associations.  Each cut uses a
pair of complementary squared-E-surface ratios built from actual E-surfaces in
that cut's selected-orientation catalog.  Three cuts use `(8,12,14)` as their
secondary surface.  Cut `(2,4,10,13)` has no other active threshold, so its
only legal secondary is the owning cut equation itself.  The document also
forces `(8,12,14)` under cut `(2,6,10)` from its legacy `[3,7]` maximal
subspace into the shared `[7]` frame used by the same threshold under the two
other cuts participating in the `q13` cancellation.  All unrelated thresholds
remain implicit legacy defaults.  Multiplier expressions use the mandatory
symmetric edge-set syntax `eta(effective, eset(...))`; Symbolica therefore
canonicalizes the E-surface edge order itself.  In component contexts,
`effective` is the base sample for local pieces and the root sample for
integrated pieces, giving exactly `f(r)` and `f(r_star)`.

This shared-subspace and ratio prescription implements the cure demonstrated
for the selected-orientation correlated `q13` soft/threshold limit of cut
`(2,4,12)`.  This folder also
demonstrates finite multiplier evaluation, expanded L/I component accounting,
event cancellation, and serialization.  It does not claim that the same
prescription has yet been demonstrated to cure every other GL638 cut.

Run from the repository root:

```bash
MallocNanoZone=1 CARGO_INCREMENTAL=0 cargo run -p gammaloop-api --bin gammaloop -- \
  --clean-state -o -l info -L off -p none \
  examples/cli/epem_a_ttxh/NNLO/analyses/gl638_multiplier_components/analysis.toml \
  run generate inspect_full_cut_samples approach_safe approach_q13_correlated

python3 examples/cli/epem_a_ttxh/NNLO/analyses/gl638_multiplier_components/plot.py
```

The two coordinate-space `inspect` inputs are exactly the finite full-six-cut
samples used by the Rust runtime/save-load regression.  `approach_safe` is a
momentum-space JSON smoke trajectory around a generic point chosen away from
the simultaneous zero of either ratio denominator.  It is component evidence,
not an IR-scaling acceptance criterion.

`approach_q13_correlated` targets cut `(2,4,12)`.  With
`K*=sqrt(500^2-MT^2)(3,-4,12)/13`, `a=(47,29,-31)`,
`d=(211,-157,193)`, and `λ=10^-2 t`, its path is

```text
K0 = K* + λ a
K3 = K* + λ (a+d)
q13 = K0-K3 = -λ d.
```

Consequently `|q13|`, `eta(eset(5,10))`, and
`eta(eset(5,12,13))` all vanish linearly at the common soft/threshold point.
The same logarithmic `|λ|` samples are evaluated on both sides of zero, with
the undefined midpoint omitted.  The two right threshold CTs are both solved
in the common one-loop subspace `[10]`.  At the four smallest scales on both
sides, the original contribution, threshold-CT sum, and final total all fit a
power of `-1.000`; the original and CT leading terms cancel by about two
orders of magnitude.

`plot.py` creates `multiplier_decomposition.pdf` for the generic finite path
and `q13_correlated_cure.pdf` for the demonstrated cut-`(2,4,12)` cure.  They
include the event closure, one-sided components, and the LL/LI/IL/II terms for
the two independent `(7,8)` variants.  The correlated report begins with a
double-sided logarithmic fit of the three kinematic distances above, followed
by the mirrored component pages.  The schema-v2 JSON retains every
reconstructed component, and the regression tests independently verify
equality to the reported event weight.

In the expanded-component metadata, “bare component” means the unweighted
threshold CT before its user multiplier.  It does not mean threshold
subtraction is disabled.

Integrated UV generation is disabled in this direct-import card because the
existing GL638 Vakint reconstruction fails before threshold generation with
`Could not solve Inconsistent`.  Local UV physics and both local/integrated
threshold components remain enabled; this independent CLI gap is recorded
rather than hidden.
