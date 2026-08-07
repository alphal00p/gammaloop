# GL638 multiplier and expanded-component decomposition

This finite-path analysis regenerates the dressed GL638 graph from
`../../graphs/GL638.dot` with the verified selected orientation and all six
process-valid cuts.  Its purpose is diagnostic rather than asymptotic: it
exposes the independently named `intrinsic_1l` and `embedded_2l` variants, all
local/integrated Cartesian-product components, their weighted and unmultiplied
values, and every runtime multiplier occurrence.

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

`multiplier_decomposition.pdf` contains a finite total page, a closure and
component report for cut `(2,4,12)`, and occurrence-resolved multiplier pages
for the two left variants paired with the `eta(5,10)` right threshold.  The
plotting wrapper performs no rendering itself: every page is produced by
`assets/plot_approach_result.py` and only then merged.

Integrated UV generation is disabled because direct-import Vakint
reconstruction for GL638 currently fails before threshold preprocessing.
Local UV subtraction and both local and integrated threshold counterterms stay
enabled.  The generated state, JSON, and PDF are ignored by Git.
