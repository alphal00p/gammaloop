# GL638 multiplier and expanded-component decomposition

This finite-path analysis regenerates the dressed GL638 graph from
`../../graphs/GL638.dot` with the verified selected orientation and all six
process-valid cuts.  Its purpose is diagnostic rather than asymptotic: it
exposes the independently named `intrinsic_1l` and `embedded_2l` variants, all
local/integrated Cartesian-product components, their weighted and unmultiplied
values, and every runtime multiplier occurrence. It samples 50 logarithmic
points on each signed deformation branch (100 evaluated points in total).

From the repository root, generate the integrand and real `approach` output:

```bash
MallocNanoZone=1 CARGO_INCREMENTAL=0 target/debug/gammaloop \
  --clean-state -o -l info -L off -p none \
  examples/cli/epem_a_ttxh/NNLO/analyses/gl638_multiplier_decomposition/analysis.toml \
  run all
```

Create the curated multipage PDF with:

```bash
python3 examples/cli/epem_a_ttxh/NNLO/analyses/gl638_multiplier_decomposition/plot.py
```

`multiplier_decomposition.pdf` contains a finite total page and one
occurrence-resolved local-local selector page for both duplicated variants
paired with the `eta(5,10)` right threshold. The correlated analysis carries
the complementary unmultiplied `L/I` and `LL/LI/IL/II` component pages; they are
more informative there because that path reaches the physical soft/threshold
limit. The finite wrapper performs no rendering itself: every page is produced
by `assets/plot_approach_result.py` and only then merged.

This report is structural rather than an IR-scaling demonstration. Because the
current `intrinsic_1l` expression contains its own threshold equation
`eta(star, eset(7,8))`, that equation vanishes by construction at the CT root:
the observed selector is zero (up to root-solver residue), while
`embedded_2l` is unity. The correlated cure report therefore uses the actual
cross-cut soft-cancellation components instead of fitting these flat selector
curves.

Integrated UV generation is disabled because direct-import Vakint
reconstruction for GL638 currently fails before threshold preprocessing.
Local UV subtraction and both local and integrated threshold counterterms stay
enabled.  The generated state, JSON, and PDF are ignored by Git.
