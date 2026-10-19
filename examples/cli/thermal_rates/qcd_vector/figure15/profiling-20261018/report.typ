#set page(paper: "a4", margin: 20mm)
#set text(size: 10pt)

= Thermal and vacuum QCD NLO CPU profile

These measurements are a historical baseline for `116767bd3`, before the
binary64-first root certification in `75a0ced4e`. Source locations and costs
below refer to that baseline. See #link("../baseline.typ")[the baseline record].

The dominant thermal cost is solving the additional signed H-surface cut
equations. It is already present in native Double evaluations. Fixed
2048-bit directed interval arithmetic certifies these roots; native Quad
escalation then repeats much of the same geometry work.

The experiment profiles the transverse NLO integrand at
$(E/T,k/T)=(1,0.5)$ and $(1.05,1)$, using source commit
`116767bd3ed30491a2b15559d228171d20082abe` and the same immutable binary as
the Figure 15 scan. Thermal and vacuum states have identical models and
effective runtime settings after the explicit overrides. Their medium modes
are `thermodynamic_equilibrium` and `vacuum`, respectively; vacuum generation
was performed separately. The earlier
#link("../vacuum-timing-20261018/report.typ")[matched timing report] preserves
the generation audit and random-sample comparison.

The last shared input boundary is the model, graph topology, projector, and
requested subtraction settings. Medium-specific generation first changes the
retained bath-crossing residues: ordinary cut counts remain three and four
for `GL0`/`GL2`, but thermal orientation counts are 32/32 versus 14/18 in
vacuum, and evaluator counts are 15/21 versus 5/5. Both routes then run the
same sampling, stability dispatch, root preparation, residue/parameter
construction, and expression-evaluation implementations on those different
generated branch sets. The comparison preserves the factorized graph
numerators. The controlled precision ladders and separate graph probes
distinguish branch-content cost from escalation; they do not assert
identical internal residue content between medium modes.

== Measured call-stack attribution

Percentages below refer to CPU samples selected by native evaluation ancestry,
including warmup and one diagnostic inspection. Root, H-root, enclosure, and
MPFR categories are inclusive or nested and must not be added together.
They describe where time is spent within each case, not ratios of absolute
thermal and vacuum costs.

#table(
  columns: 6,
  [Production case], [Samples], [All roots], [H roots], [MPFR/GMP leaves], [Expressions],
  [Ordinary thermal], [526], [94.7%], [93.0%], [57.8%], [3.4%],
  [Ordinary vacuum], [933], [48.7%], [0.0%], [18.5%], [21.2%],
  [Near-cone thermal], [1,060], [89.5%], [87.7%], [57.5%], [7.4%],
  [Near-cone vacuum], [959], [48.4%], [0.1%], [21.9%], [21.4%],
)

The hot path is
`CrossSectionGraphTerm::evaluate` → `Hsurface::positive_radius_roots` →
`EnergySurfaceRay::positive_radius_roots` → interval multiplication,
division, and comparison. Thermal interval helpers occur on 70.7% and 63.4%
of production CPU stacks. The leading self-time leaves are `mpfr_mul`
(17.5% and 17.7%), `mpfr_cmp3` (7.6% and 7.9%), and GMP multiplication,
division, and limb-copy primitives. Allocation/free leaves account for
16.3% and 14.9%, often beneath the same interval work. Square roots occur
in the algorithm but are not among the top thirty sampled self-time leaves.

Parameter building accounts for only 0.4% and 1.6% of production thermal
samples. Evaluated occupation functions are included in the expression
category. They do not account for the large cut/kinematic core observed in
the earlier random-sample timing test.

The source explains the thermal/vacuum difference:

- `cff/hsurface.rs:49` routes signed positive and negative energy terms to
  the radial solver.
- `cff/radial.rs:156` sets enclosure precision to 2048 bits independently
  of the native Double/Quad lane. Starting at line 772, signed rays use
  compactification of the entire positive radius, directed interval bounds,
  and adaptive isolation before native root refinement. The implementation
  permits up to 65,536 visited cells; this experiment does not measure the
  number visited.
- `cff/radial.rs:741` provides a cheaper convex positive-energy shortcut
  when the origin is interior, avoiding subdivision. It retains an interval
  origin check, so some MPFR work also remains in vacuum.
- `integrands/process/cross_section/mod.rs:3144` prepares every active
  E/H cut-group root before the orientation/evaluator loop. It continues
  preparing other cut roots after one fails, retaining observations for
  later precision attempts. Selecting an orientation would not isolate this
  root-solving cost.

These paths are relative to `crates/gammalooprs/src/`. Threshold counterterms
and soft counterterms are disabled in the benchmark. The observed H-root
cost is therefore not evidence of threshold-counterterm evaluation. No
implementation or numerical acceptance policy was changed.

== Separating precision escalation

Each profile uses the same fixed six-dimensional sampling coordinate
`[0.17, 0.31, 0.43, 0.29, 0.53, 0.71]`. The production ladder and the
Double-only and Quad-only diagnostic ladders retain the original stability
rotation and norm checks. Timings are means from eight benchmark batches;
the benchmark disables event retention. “Core” is the exclusive
`physical other` bucket: integrand time minus evaluator, canonical physical
preparation, and event processing. Wrapper/sampling overhead remains in Total.

At the ordinary point, all controlled lanes accept the evaluation:

#table(
  columns: 5,
  [Medium / lane], [Core], [Evaluator], [Total], [Final status],
  [Thermal Double], [51.16 ms], [1.80 ms], [53.41 ms], [Stable],
  [Thermal Quad], [56.93 ms], [8.02 ms], [65.46 ms], [Stable],
  [Vacuum Double], [1.19 ms], [0.36 ms], [1.84 ms], [Stable],
  [Vacuum Quad], [2.42 ms], [3.78 ms], [6.57 ms], [Stable],
)

The accepted Double thermal core is about 43 times the vacuum core. Moving
thermal arithmetic to Quad adds about 11% to core time and 23% to total time
at this coordinate. Thermal Double-only stacks still spend 92.9% in root
solving and 60.6% in MPFR/GMP leaves. Thus MPFR samples cannot be equated
with escalation.

Near the light cone, the accepted production evaluation takes 142.24 ms
in the events-disabled benchmark. A separate metadata inspection measures
64.997 ms in the rejected Double attempt followed by 78.024 ms in the
accepted Quad attempt. The inspection enables events to expose metadata,
so these attempt times are an approximate decomposition rather than exact
subtraction from the benchmark total. Its total is 143.663 ms.

Call stacks independently allocate 44.3% of near-cone production CPU to
Double and 54.9% to Quad, with the small remainder in preparation or
unclassified wrapper work. Root solving accounts for 42.9 and 46.6
percentage points respectively: roughly 97% of Double time and 85% of Quad
time. Escalation adds the second expensive pass, around 55% of the observed
production time at this coordinate.

Both forced single-lane thermal near-cone evaluations are rejected
(`is_nan=true`, `Unstable(0)`). Their displayed zeros are placeholders;
they are not accepted physics results or valid substitutes for the
production Quad cost. The diagnostic rejection names `GL0`, cut group 6,
identity rotation, with a root near 36.6028168038367. Double's residual is
about $4.22 times 10^(-15)$; standalone Quad's is about
$-1.42 times 10^(-30)$. The source's cross-precision roundoff-rescue policy
can use a matching earlier Double observation to certify Quad
(`utils/newton_solver.rs:262–318,403–410,465–507`). This is consistent with
the production/standalone difference; the specific rescue branch was not
directly logged.

== Graph isolation and scope

A second, unprofiled momentum-space probe uses
`[0.17, -0.31, 0.43, -0.29, 0.53, 0.71]`, a four-second target per benchmark, and the
Double-only diagnostic ladder. It isolates graph indices 0 and 1, corresponding
to `GL0` and `GL2`. These coordinates differ from the profiled sampling
coordinate, so the results must not be combined into one decomposition.

#table(
  columns: 5,
  [Kinematics / medium], [All core], [GL0 core], [GL2 core], [GL2 share of isolated sum],
  [Ordinary thermal], [8.94 ms], [3.41 ms], [5.58 ms], [62%],
  [Ordinary vacuum], [0.265 ms], [0.101 ms], [0.159 ms], [61%],
  [Near-cone thermal], [10.05 ms], [4.52 ms], [5.52 ms], [55%],
  [Near-cone vacuum], [0.267 ms], [0.102 ms], [0.161 ms], [61%],
)

Both selected graph paths incur substantial thermal overhead at this
diagnostic probe. The graph benchmarks did not retain full acceptance
metadata. Neither this selector
nor sampled stacks expose CPU time by individual runtime cut-group ID;
group 6 is identified as the rejection trigger, not established as the
dominant time-consuming cut. These are two fixed transverse coordinates,
not a statistical profile of the full scan or a longitudinal profile.

The best-supported optimization target is signed radial isolation and its
interval allocation/arithmetic, followed by repeated certification across
stability and precision passes. This report does not establish that any
reuse or precision reduction would preserve the certification policy.

== Profiling quality and reproducibility

The twelve recordings use Linux `perf`, 99 Hz userspace CPU-clock samples,
and 16 KiB DWARF call stacks. Evaluation loops run synchronously on the main
thread, so per-thread recording covers the benchmark workload. Six workers
were pinned to CPUs 242–247, with Rayon, OpenMP, and BLAS limited to one
thread per process. Offline readers used at most two CPUs.

Across all cases, 12,295 raw CPU samples contain 11,208 selected native
evaluation samples, zero lost-sample records, and zero unknown leaves.
Short/truncated caller chains can still limit inclusive attribution;
the percentages are sampled estimates, not exact instrumented timers.
Per-thread recordings omit timestamps, so startup is filtered by evaluation
ancestry rather than a time window.

The profile benchmark's twelve-second duration is a target used with warmup
to estimate the sample count. It does not guarantee twelve seconds of steady
evaluation. Ordinary thermal production estimated 111 timed samples from
slower warmup, then completed them in about 6.2 seconds; this explains its
smaller profile of 601 raw and 526 selected samples. Its recording has no
throttle/unthrottle or lost-sample records. The recorded sample counts,
rather than the requested duration, establish each profile's resolution.

The native jobs completed, but `perf record` became stuck in optional
build-ID/inline postprocessing of the large debug binary. Recorder/controller
processes were stopped only after native EXIT and all benchmark outputs were
verified. Original recordings remain untouched. Separately recovered copies
complete only the container header's data size and clear unfinished optional
feature flags. Raw event bytes are unchanged, all record boundaries validate,
and there are no trailing bytes. The original controller manifests retain
their failed/interrupted status; recovery is documented separately.

#link("comparison.json")[Machine-readable results] include attempt status,
timing errors, raw artifact paths and hashes.
#link("timing-summary.csv")[Timing CSV],
#link("cpu-attribution.csv")[CPU CSV],
#link("sampled-attribution-summary.json")[full attribution], and
#link("raw-recovery-audit.json")[recovery audit] are committed here.
The per-case benchmark batches, diagnostic event arrays, hot self-time reports,
raw recordings, and compressed stacks are preserved in the verified
#link("../archive.json")[archive], under
`tmp/qcd-cut-profile-20261018T201221Z/`. The original checkout copy also
preserves the compact case files and command drivers under
`profiling-20261018/`.

The archived `profile_matrix.py` preserves exact CLI commands;
`extract_profiles.py` documents recovery, filtering, symbol attribution,
and precision classification. Future captures should add `--no-buildid`
to avoid the recorder's postprocessing hang; offline readers already use
`--no-inline`. Original manifests retain the failed initial mmap attempt,
recorder termination status, and successful graph probes.

All measurement workers stopped. The main scan subsequently finished within
its 60-worker ceiling and original hard cutoff of 19 October 2026 at
07:19:29 Europe/Zurich (05:19:29 UTC).
