= Soft counterterm evaluator construction memory

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

== Self-contained reproducer

`mre.rs` generates a small example without GammaLoop or captured input files.
Its embedded Cargo manifest pins the same Symbolica revision. Run each mode
in a fresh process:

```sh
rust-script mre.rs inline
rust-script mre.rs shared
rust-script mre.rs raw
rust-script mre.rs inline 4096 128
rust-script mre.rs inline 16 32 --verbose
```

The last two commands respectively increase the stress and expose individual
CSE passes. The default is 512 function calls and a recurrence depth of 128.
The script builds $H_0(x) = x$, $H_j(x) = 1 / (1 + x H_(j-1)(x))$, and
$F(u, x) = 1 / (u + H_D(x))$. Both function arguments are used. `inline`
constructs $sum_(i=1)^N F(i, x)$; `shared` writes the equivalent terms
$1 / (i + H_D(x))$ directly. Both retain the same alias definitions. `raw`
measures the inlined program before optimization, with matched preprocessing.

The generated source and final program grow linearly with $N + D$, but
ordinary inlining creates $3 N D + 3 N - 1$ arithmetic operations, including
inversions. CSE reduces these to $3 D + 3 N - 1$ through $3 D$ successful
passes when there are duplicated calls. At the default size this is
198,143 operations becoming 1,919, with about 9.5 kB of symbolic input.
Each mode reports construction time and Linux process peak RSS, then checks
three distinct numerical points against an independent recurrence oracle.
This isolates the evaluator mechanism; it is not a replacement physics graph.

A copied single file was built and run with `rust-script` in an isolated
directory, with no GammaLoop packages in its resolved dependencies. The
default inline build took 4.239 seconds and peaked at 67.7 MiB; the shared
control took 0.001488 seconds and peaked at 5.92 MiB. Both produced the same
1,919-operation program and passed all three oracle checks. The raw mode
reported 198,143 operations, and the small verbose case confirmed 96
successful CSE passes. These are single diagnostic measurements. Printed
RSS precedes raw stack compaction and numerical conversion; an external
whole-process memory measurement can include those later stages.
The figures above summarize the original run. Its raw receipts are not distributed.

== Captured GL262 replay

This diagnostic replays the scalar evaluator input captured by the existing
`evaluator_symbolica_input_dump` event. It retains the root, scalar aliases,
tagged function definitions, inlining policies, parameter order, and recorded
optimization settings. It rejects dual-vectorized inputs. Symbolica is pinned
to `578dfcb55fb0662d871456ea8ecc57c1390caa77`.

The comparison changes only `cpe_iterations`: the recorded `None` versus
`Some(0)`. In the pinned direct translator, the latter skips common-pair
elimination while retaining Horner rewriting, function inlining,
common-instruction elimination, and stack optimization. Setting Horner
iterations to zero would bypass several of these passes together and is not
an isolated Horner control.

`replay.rs` is a standalone rust-script. Given a compiled copy, run both modes
in separate processes with the provided supervisor:

```sh
python3 run.py REPLAY_BINARY INPUT.json /tmp/evaluator-replay \
  --seconds 900 --gib 100 --perf /path/to/perf
```

The optional profiler records user-space CPU stacks at 49 Hz. CPU samples
locate work, not live allocations; correlate them with the time-stamped RSS
samples before attributing a memory peak. The supervisor stores binary and
input hashes, commands, logs, resource samples, exit status, and numeric
comparison results. Neither the input JSON nor the experiment outputs are
written into the repository.

Each evaluator is checked at three deterministic complex parameter vectors
at 1000-bit precision. Values are constructed from exact integer ratios.
The acceptance threshold is a relative complex infinity-norm difference
below $10^(-250)$. These vectors test agreement between optimizer variants;
they are not physical momenta and do not replace UV or soft scaling tests.
Report the construction peak separately from the later arbitrary-precision
program conversion and evaluation peaks.

The additional `raw` mode measures the program before CSE, CPE, and stack
compaction. It requires a single function root with numeric arguments, direct
translation, one Horner iteration, and no hot start. For this shape, the
builder's root-based Horner variable selection is empty. The pinned empty
scheme only collects coefficients. Applying that collection explicitly to
the root and every registered body, then using zero Horner iterations,
preserves the actual preprocessing while exposing the raw linearization.
This equivalence is specific to the checked preconditions. Raw mode counts
operations without exporting another copy of the program, and skips numeric
evaluation unless `--check-raw-numeric` is explicitly supplied for a small
input. Small oracle checks cover the three modes and coefficient collection.
The pre-CSE comparison here is validated for the captured all-inline input;
it does not establish that boundary for retained non-inline function bodies.

The original GL262 capture and its provenance are retained under
`/tmp/soft-ct-evaluator-memory-2026-10-01/capture`. Its SHA-256 is
`6a749f060abf887857843e922aefc81ef818d5b86df19824ac71f5cfd66f204a`.
All recorded input-size and function-count measurements match the previous
GL262 run: 494 scalar aliases, 98 recorded functions, and 220 parameters.

== Results

#table(
  columns: 4,
  [Construction endpoint], [Arithmetic operations], [Function calls], [Peak RSS],
  [Raw linearization], [135,694,111], [486,853], [22.99 GiB],
  [CSE, with CPE disabled], [3,822,857], [42], [30.26 GiB],
  [Full optimization], [3,428,224], [42], [30.26 GiB],
)

The input is compact, but inlining produces about 39.6 times the final
arithmetic operation count. Raw construction takes 87.79 seconds plus
0.45 seconds of explicit preprocessing. CPE-disabled and full construction
take 356.43 and 369.71 seconds respectively. These are single diagnostic
runs on a shared host, not repeated timing benchmarks.

Both CPE variants pass the three 1000-bit comparisons; the maximum normalized
residual is approximately $1.53 times 10^(-299)$. The full raw program was
measured without numerical evaluation. Small fixtures separately validate
its preprocessing and arithmetic against the other modes.

CPU profiles place both peaks in `remove_common_instructions`, before CPE.
The observed CSE interval spans roughly 262 seconds; CPE occupies about five
seconds much later. A deeper stack sample resolves hash-table insertion and
instruction-vector allocation/copy under CSE. Clock alignment uncertainty of
less than 0.75 seconds does not affect that phase attribution. These samples
identify active routines, not the allocation sizes contributing to live RSS.

The pinned CSE implementation allocates a hash table and a replacement
instruction vector while retaining the old vector and cloning surviving
operand arrays. Its lookup keys refer to old operand indices; renaming is
applied afterward. Thus equivalence of dependent instructions propagates
through repeated full passes. Avoiding that duplicated intermediate work is
the relevant optimization target; disabling CPE leaves the peak unchanged
and increases the final operation count.

Raw receipts, profiles, and captures are not distributed. Rerun the commands
above with a newly captured input to collect fresh results.
No production defaults or dependency pins changed.
