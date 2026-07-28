```bash
gammaloop/workspaces/dev is 📦 v0.3.3 via 🐍 v3.13.12 via 🦀 v1.93.1 via ❄️  impure (nix-shell-env)
❯ cargo run -p gammalooprs --example bnl_evaluator_atom_mwe --profile dev-optim -- examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_unfiltered_pre_network.sym --alias-scalars 8
    Finished `dev-optim` profile [optimized + debuginfo] target(s) in 0.31s
     Running `target/dev-optim/examples/bnl_evaluator_atom_mwe examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_unfiltered_pre_network.sym --alias-scalars 8`
input
╭────────┬─────────────────────────────────────────────────────────────────────────────────────╮
│ metric │ value                                                                               │
├────────┼─────────────────────────────────────────────────────────────────────────────────────┤
│ path   │ examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_unfiltered_pre_network.sym │
│ bytes  │ 2038141                                                                             │
╰────────┴─────────────────────────────────────────────────────────────────────────────────────╯
symbolica_parse
╭─────────┬─────────────╮
│ metric  │ value       │
├─────────┼─────────────┤
│ elapsed │ 46.742583ms │
│ terms   │ 1           │
│ bytes   │ 1410951     │
╰─────────┴─────────────╯
network_parse
╭────────┬───────╮
│ metric │ value │
├────────┼───────┤
│ status │ start │
╰────────┴───────╯
network_parse
╭─────────┬─────────────╮
│ metric  │ value       │
├─────────┼─────────────┤
│ status  │ done        │
│ elapsed │ 174.06925ms │
│ nodes   │ 98          │
│ edges   │ 179         │
│ tensors │ 24          │
╰─────────┴─────────────╯
scalar_store
╭────────────────┬─────────╮
│ metric         │ value   │
├────────────────┼─────────┤
│ count          │ 23      │
│ top_level_sums │ 0       │
│ total_terms    │ 23      │
│ max_terms      │ 1       │
│ total_bytes    │ 1408133 │
│ max_bytes      │ 953984  │
╰────────────────┴─────────╯
scalar_aliases
╭─────────────────┬─────────╮
│ metric          │ value   │
├─────────────────┼─────────┤
│ threshold_bytes │ 8       │
│ created         │ 14      │
│ terms           │ 14      │
│ bytes           │ 1408099 │
│ max_bytes       │ 953984  │
╰─────────────────┴─────────╯
network_execute
╭────────┬───────╮
│ metric │ value │
├────────┼───────┤
│ status │ start │
╰────────┴───────╯
network_execute
╭─────────┬──────────────╮
│ metric  │ value        │
├─────────┼──────────────┤
│ status  │ done         │
│ elapsed │ 445.382458ms │
╰─────────┴──────────────╯
result_scalar
╭─────────┬─────────╮
│ metric  │ value   │
├─────────┼─────────┤
│ status  │ ok      │
│ elapsed │ 1.041µs │
╰─────────┴─────────╯
result_aliased_root
╭─────────┬───────────╮
│ metric  │ value     │
├─────────┼───────────┤
│ aliases │ 14        │
│ terms   │ 1         │
│ bytes   │ 102572537 │
╰─────────┴───────────╯
result_alias_resolve
╭─────────┬───────╮
│ metric  │ value │
├─────────┼───────┤
│ skipped │ true  │
╰─────────┴───────╯
summary
╭──────────────────────┬─────────┬──────────────┬───────┬───────────┬──────────────────────────────────────────────────────────────────────────────────────────╮
│ stage                │ status  │ elapsed      │ terms │ bytes     │ details                                                                                  │
├──────────────────────┼─────────┼──────────────┼───────┼───────────┼──────────────────────────────────────────────────────────────────────────────────────────┤
│ input                │ read    │ -            │ -     │ 2038141   │ path=examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_unfiltered_pre_network.sym │
│ symbolica_parse      │ done    │ 46.742583ms  │ 1     │ 1410951   │ -                                                                                        │
│ selection            │ all     │ 28.5µs       │ 1     │ 1410951   │ -                                                                                        │
│ network_parse        │ done    │ 174.06925ms  │ -     │ -         │ nodes=98 edges=179 tensors=24                                                            │
│ scalar_store         │ done    │ -            │ 23    │ 1408133   │ count=23 top_level_sums=0 max_terms=1 max_bytes=953984                                   │
│ scalar_aliases       │ enabled │ -            │ 14    │ 1408099   │ threshold_bytes=8 created=14 max_bytes=953984                                            │
│ network_execute      │ done    │ 445.382458ms │ -     │ -         │ -                                                                                        │
│ result_scalar        │ ok      │ 1.041µs      │ -     │ -         │ -                                                                                        │
│ result_aliased_root  │ done    │ -            │ 1     │ 102572537 │ aliases=14                                                                               │
│ result_alias_resolve │ skipped │ -            │ -     │ -         │ -                                                                                        │
╰──────────────────────┴─────────┴──────────────┴───────┴───────────┴──────────────────────────────────────────────────────────────────────────────────────────╯
```
