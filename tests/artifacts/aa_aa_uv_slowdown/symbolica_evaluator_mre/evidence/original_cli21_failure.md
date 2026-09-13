The original GL262 run failed inside Symbolica's **evaluator construction** call.
Evaluator compilation was disabled. The full portable input is being recovered;
this historical failure has **not yet been reproduced by the standalone package**.
The machine-readable record is [original_cli21_failure.json](original_cli21_failure.json).

The application supplied one valid, contracted physical scalar expression with an
Atom payload of **19,982,016,710 bytes**, assembled from 3,050 independent input
tensor summands. UV generation and tensor preprocessing had completed. The failing
call was equivalent to:

```rust
expression
    .evaluator(&parameters)
    .function_map(function_map)
    .optimization_settings(settings)
    .build()
```

No evaluator source compilation, shared-library loading, or JIT compilation was
requested. The expression came from a three-loop photon-scattering amplitude
with local UV subtraction. Integrated UV and threshold counterterms were disabled.

The exact panic message was:

```text
advance out of bounds: the len is 1721315294 but advancing by 10311249876
```

The reported location was `bytes-1.11.1/src/lib.rs:170`, and the process exited
with code 101. No backtrace was enabled in that original run. The preserved
[panic excerpt](original_cli21_panic.txt) removes ANSI color sequences only.
The precise internal cause remains unconfirmed.

The immutable run used:

| Item | Value |
|---|---|
| GammaLoop revision | `a8b4c9d96c3994f77bba80826e16715d53ddbfc9` |
| CLI SHA-256 | `91c9e3389d9b7c3cd614890d24d46422cb45c4339efa3333db8858f958733fe2` |
| Symbolica version | `2.2.0` |
| Symbolica revision | `4d0a833eb8e059d1f95bdae5abed2559830b235f` |
| Build profile | `dev-optim`, debug assertions enabled |
| Rust compiler | `rustc 1.97.0 (2d8144b78 2026-07-07)` |
| Evaluator compilation | `false` |
| Direct translation | `true` |
| Horner iterations | `1` |
| Common-pair-elimination rounds | `5` |
| Evaluator cores / worker limits | `1` |
| Maximum Horner scheme variables | `500` |
| Maximum common-pair cache entries / distance | `1,000,000` / `1,000` |

`compile: true` in the original JSONL stage events is the **logging tag**
`#compile`; it does not represent the evaluator setting. The frozen card explicitly
sets `compile = false`. Explicit orientation summation also disables the optional
iterative and summed evaluator variants before construction.

The recorded chronology on 2026-09-13 UTC was:

| Event | Time / duration |
|---|---|
| Process started | 14:18:09.984759 |
| Tensor preprocessing completed | 15:12:41.096 |
| Complete preprocessing interval | 2,314.506 s |
| Expression preparation completed | 15:15:37.070 |
| Expression preparation interval | 175.972 s |
| Symbolica `.build()` entered | 15:15:37.070 |
| Process failure recorded | 15:18:08.289353 |
| Whole process wall time | 3,598.305 s |

There was no subsequent successful build event, no numerical-program conversion,
and no saved evaluator. This was a profiled diagnostic with earlier concurrent
build/check activity, so these observations are not final benchmark measurements.

Memory remained below the configured guard:

| Boundary | RSS | Lifetime high-water mark |
|---|---:|---:|
| After tensor preprocessing | 46,050,947,072 B | 68,617,658,368 B |
| At Symbolica build entry | 46,050,967,552 B | 68,617,658,368 B |
| Complete process peak | 120,480,550,912 B | 120,480,550,912 B |

The first two samples were observed respectively 42.844 ms and 8.255 ms after
their logged boundary events. The **500,000,000,000-byte** memory guard and
7,200-second wall guard did not intervene.

Separate controls established that smaller valid inputs build successfully:

| Control | Result | Symbolica build | Whole diagnostic peak |
|---|---|---:|---:|
| One physical GL262 component, 3,389,407,083-byte Atom; direct translation, H1/CPE5 | Passed | 131.555 s | 27,143,188,480 B |
| Small polynomial using the physical parameter context; direct translation, H1/CPE5 | Passed; exact value 36 | 0.001003 s | 558,964,736 B |
| Same small polynomial; tree translation, H1/CPE5 | Passed; exact value 36 | 0.000457 s | 538,660,864 B |

The physical component's portable file is 3,389,426,376 bytes, SHA-256
`66240aab1ee9115af7d8e0c4e5c5291b2cf74ef74a931efe07ce09db84e7fce0`.
Its successful program contains 80,026 additions, 67,213 multiplications,
21 inversions and 11 function calls. The control reconstructs 220 ordered
parameters and 63 function-map definitions, plus the two imaginary-unit aliases.
These controls used the existing GammaLoop-linked constituent-operation probe;
they are separate from validation of this new standalone package.

The large component's alternate tree-translation diagnostic was interrupted by
the supervisor. It is neither a successful control nor a reproduction of the
original panic. The JSON evidence retains this distinction and hashes the source
receipts, frozen card and manifest. No giant log, expression dump or credential
is copied into this evidence directory.
