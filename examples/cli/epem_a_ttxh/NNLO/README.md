# ttH NNLO examples

The [threshold metadata guide](../../../../docs/products/gammaloop/content/threshold-subtraction.typ)
and [sampling guide](../../../../docs/products/gammaloop/content/sampling.typ)
explain the settings and their validation limits.

- `epem_a_tth_NNLO_test_GL297.toml` and `epem_a_tth_NNLO_test_GL638.toml`
  provide named generation, inspection, approach and integration blocks.
- `gl297_*.toml` and `gl638_*.toml` compare LMB and surface-channel selections
  in separate states and workspaces. Each card documents its run command and
  parameter overrides; start with a small sample budget.
- `../NNLO_experiment/gl134_*.toml` isolates one tree-level physical cut.
- `monitor_tth_nnlo_integration.py --list-available-graphs` discovers local
  comparison workspaces; `--graph gl638` displays their iteration results.
  The monitor requires Python, `prettytable` and `colorama`.

Run cards from the repository root. Generated states and integration outputs
are local artifacts and are not included in the repository.

The `graphs/` folder contains the complete **166-master-graph unfiltered ttH NNLO
selection**, including raised propagators and cuts, with explicit threshold metadata;
29 graphs have no eligible thresholds. Import with the process specification in the run
cards to select the physical ttH cuts. The graph prescriptions were checked with complete
orientation sums and local 3D plus integrated UV subtraction at 1000 GeV, with
`mu_r = 91.188 GeV` and `m_uv = 50 GeV`. New kinematics need their own checks.

Use `display integrands --show_threshold_subtraction`, `inspect` and `approach` to review
or refine a prescription, as described in the threshold guide. Regenerate states predating
**state format 20**. The scoped soft/threshold checks do not establish bounded weights,
finite variance or global IR safety. Deep imaginary-part stability remains uncertified for
GL474 and GL554; their physical real-part checks pass.
