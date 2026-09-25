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
