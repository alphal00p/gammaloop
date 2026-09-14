# Verified results — 2026-09-14

| Reproducer | Original fork 4d0a833 | Official main ba373713 |
|---|---|---|
| Complete GL262 evaluator build | Same buffer panic | Same buffer panic, corrected S-format license |
| Nonsymmetric f(x,y) import | Incorrect f(y,x), assertion fails | Correct f(x,y), passes |
| Complete GL00 builder control | Pass | Pass; same four operation counts |

Both packages identify as Symbolica 2.2.0; exact commits distinguish them:

- [Original fork](https://github.com/alphal00p/symbolica/tree/4d0a833eb8e059d1f95bdae5abed2559830b235f)
- [Official main](https://github.com/symbolica-dev/symbolica/tree/ba3737137c2a2ccd7bb39f0441837d38ec867e78)

| Complete standalone GL262 run | Wall, including input load/teardown | Builder reached at | Pre-builder RSS | Sampled peak tree RSS | Exit |
|---|---:|---:|---:|---:|---:|
| Original fork | 179.911 s | 29.703 s | 19.988 GB | 93.856 GB | 101 |
| Official main, corrected license | 186.136 s | 31.310 s | 19.989 GB | 94.133 GB | 101 |

These are diagnostic single runs, not a timing comparison of successful
construction. Both enforce 500 GB and run with H1/CPE5, one core and debug
assertions. No evaluator C++/JIT compilation is attempted. Upstream has a
compatible inlining-policy addition in its function-map binary schema; its
loader verifies the original IDs and raw Atom bytes after adapting this field.
The complete 19,982,016,710-byte expression is unchanged, SHA256
`f3f36f029ae973d75a8238b4654be4a30c8a1e1210d0441f9368bc8b7c31b3a6`.

Exact failure on both revisions:

```text
advance out of bounds: the len is 1721315294 but advancing by 10311249876
```

The upstream backtrace places the failure in:

```text
EvaluatorBuilder::build                evaluate/function_map.rs:461
AtomView::to_evaluator                 evaluate/tree.rs:312
Horner transformation                 evaluate/tree.rs:470
horner_scheme                         collect.rs:1498
horner_scheme_impl                    collect.rs:1526–1530
normalize (outer sum)                 normalize.rs:1542
normalize (product iteration)         normalize.rs:664
ListIterator::next                    atom/representation.rs
bytes::Buf::advance                   panic
```

Optimized frames may name the multiple-Horner search closure even with H1:
the source's H1 route goes directly through `horner_scheme`. This does not mean
that a stochastic search or a different setting ran.

## Strong lead for the Symbolica author

Product/function payload lengths are stored as u32; sum lengths support u64.
[`Mul::extend`](https://github.com/symbolica-dev/symbolica/blob/ba3737137c2a2ccd7bb39f0441837d38ec867e78/src/atom/representation.rs#L1212)
still narrows the actual product payload length with `as u32`. Horner builds
sum-times-variable products before normalizing them. The observed lengths obey:

```text
10,311,249,876 - 1,721,315,294 = 2 * 2^32 - 10
```

This strongly suggests a large Horner-generated product whose stored length
wrapped modulo 2^32, followed by traversal of its intact large sum child through
a truncated product view. The backtrace detects that invalid view; it does not
by itself identify the first overflowing write. No Symbolica patch or forced
small-input simplification was used to obtain the reproduced failure.

## Licenses and earlier attempts

The corrected S-format key is accepted under `SYMBOLICA_LICENSE`; the value is
not embedded here. The complete upstream run above and the final import test
were repeated with it, without the outdated-license warning.
Earlier upstream attempts used the subsequently superseded hash-delimited key.
One was stopped intentionally at 39 s; another reproduced the same panic in
195.305 s under restricted single-thread mode. They remain recorded in the
original workspace, and are not the licensed result reported above.

The final tiny import test source is byte-identical to the original reproducer.
`evidence/import_upstream/` contains its successful separate writer/read logs and
receipt. The old fork's small failure is independently documented in the
original standalone package; the large buffer failure is a different issue.
