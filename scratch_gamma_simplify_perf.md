# Gamma simplify aa_aa standard candle

Command:

```bash
cargo test --package gammalooprs --test aa_aa --profile dev-optim -- aaa --exact --nocapture --include-ignored
```

Date: 2026-05-26

## Baseline before gamma_simplify optimization

Cold run, including compilation:

| metric | value |
| --- | ---: |
| cargo build/test wall (`real`) | 37.72s |
| cargo build/test user | 90.66s |
| cargo build/test sys | 5.15s |
| Rust test harness | 5.53s |

Warm run 1:

| expression | spenso time | symbolica time | average time per evaluation |
| --- | ---: | ---: | ---: |
| Concretized | 2.066083ms | 91.654666ms | 1.749us |
| Simplified | 2.046433958s | 760.570791ms | 16.513us |
| Expanded | 705.361459ms | 677.960459ms | 6.069us |

| metric | value |
| --- | ---: |
| cargo test wall (`real`) | 8.75s |
| cargo test user | 7.74s |
| cargo test sys | 0.51s |
| Rust test harness | 7.96s |
| gamma simplify wall printed by test | 3.333455416s |

Warm run 2:

| expression | spenso time | symbolica time | average time per evaluation |
| --- | ---: | ---: | ---: |
| Concretized | 3.778ms | 87.460458ms | 1.348us |
| Simplified | 2.244166042s | 902.464959ms | 14.109us |
| Expanded | 726.108667ms | 791.817292ms | 4.901us |

| metric | value |
| --- | ---: |
| cargo test wall (`real`) | 9.29s |
| cargo test user | 8.13s |
| cargo test sys | 0.51s |
| Rust test harness | 8.48s |
| gamma simplify wall printed by test | 3.39620575s |

## After gating gamma5/gamma0/projector special passes

Change under test: collect factor kind metadata while parsing chain/trace
factors and only run the gamma5/gamma0/projector-specific rewrite scans when
the corresponding factor kind is present.

Post-change run 1 rebuilt downstream crates, so only the harness and table
timings are useful:

| metric | value |
| --- | ---: |
| incremental rebuild | 1m32s |
| Rust test harness | 7.11s |
| gamma simplify wall printed by test | 2.867065s |

Warm run 1:

| expression | spenso time | symbolica time | average time per evaluation |
| --- | ---: | ---: | ---: |
| Concretized | 5.550583ms | 64.98275ms | 1.161us |
| Simplified | 1.491504208s | 647.043416ms | 8.904us |
| Expanded | 536.831542ms | 680.401917ms | 5.046us |

| metric | value |
| --- | ---: |
| cargo test wall (`real`) | 7.47s |
| cargo test user | 6.11s |
| cargo test sys | 0.40s |
| Rust test harness | 6.20s |
| gamma simplify wall printed by test | 2.481456333s |

Warm run 2:

| expression | spenso time | symbolica time | average time per evaluation |
| --- | ---: | ---: | ---: |
| Concretized | 45.266041ms | 106.924583ms | 1.198us |
| Simplified | 1.996557458s | 785.732625ms | 13.414us |
| Expanded | 726.828958ms | 814.128958ms | 5.481us |

| metric | value |
| --- | ---: |
| cargo test wall (`real`) | 8.90s |
| cargo test user | 7.02s |
| cargo test sys | 0.51s |
| Rust test harness | 8.16s |
| gamma simplify wall printed by test | 3.295693208s |

The avoidable gamma5/gamma0/projector overhead sits in the rewrite phase. The
dominant standard-candle time is still Schoonschip plus evaluator construction,
so this change keeps ordinary-gamma performance within the baseline band and
removes unnecessary special-factor scans from expressions that do not contain
special factors.

## After adding no-epsilon fast path in simplify_epsilon

Change under test: `simplify_epsilon()` first checks whether the expression
contains the Levi-Civita `EPSILON_SYMBOL`. If not, it returns the owned input
without the initial Schoonschip/fixed-point epsilon pass.

Post-change run 1 rebuilt downstream crates, so only the harness and table
timings are useful:

| metric | value |
| --- | ---: |
| incremental rebuild | 18.79s |
| Rust test harness | 8.38s |
| gamma simplify wall printed by test | 3.490408833s |

Warm run 1:

| expression | spenso time | symbolica time | average time per evaluation |
| --- | ---: | ---: | ---: |
| Concretized | 2.05225ms | 61.495542ms | 1.15us |
| Simplified | 1.423377125s | 591.269667ms | 8.37us |
| Expanded | 503.102416ms | 505.496042ms | 4.238us |

| metric | value |
| --- | ---: |
| cargo test wall (`real`) | 6.01s |
| cargo test user | 5.68s |
| cargo test sys | 0.26s |
| Rust test harness | 5.68s |
| gamma simplify wall printed by test | 2.367753833s |

Warm run 2:

| expression | spenso time | symbolica time | average time per evaluation |
| --- | ---: | ---: | ---: |
| Concretized | 4.743541ms | 97.133ms | 1.309us |
| Simplified | 2.03986075s | 793.731958ms | 10.452us |
| Expanded | 769.701458ms | 719.46375ms | 5.943us |

| metric | value |
| --- | ---: |
| cargo test wall (`real`) | 8.65s |
| cargo test user | 7.67s |
| cargo test sys | 0.48s |
| Rust test harness | 7.80s |
| gamma simplify wall printed by test | 3.088369042s |

This fast path matters because the aa_aa expression has polarization symbols
printed as epsilon-like glyphs, but no Levi-Civita `epsilon(...)` tensor for the
epsilon simplifier to rewrite.
