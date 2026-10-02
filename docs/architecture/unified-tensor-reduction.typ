= Unified tensor reduction

This record separates the frozen implementation baseline, the initial unified
API measurements, and the subsequent regression correction. The public API and
its notation-conversion contract are documented
in `docs/products/idenso/content/api.typ`. The final source-qualified results,
verification gates and remaining limitations are recorded below. Correctness and
consolidation are the acceptance criteria; measured speedups are not required.
Historical alias-enabled measurements below describe the implementation preserved
on `tensor-aliases-preserved`. The current reduction API returns only `TensorExpression`.

== API and consolidation

The Rust reduction operations are `contract(ContractSettings)` and
`simplify_algebra(AlgebraSettings)`, on the existing `SymbolicTensor` abstraction
returning the same typed carrier with an Atom payload. `contract_ports` remains the explicit
binary or selected-port operation. `ContractSettings` carries representation
filters, metric and rank-one permissions, chain/trace collection preferences,
permission for local arithmetic distribution, the existing optional contraction
order, and an opt-in per-domain change budget.
`AlgebraSettings` independently carries optional `GammaSimplifySettings` and
`ColorSimplifySettings`, the epsilon permission, `AlgebraContraction`, and its
opt-in per-domain change budget. `AlgebraContraction` defaults to `Fully`; `Dots` also
normalizes scalar-product notation, `Minimal` excludes distribution of independent
sum alternatives, `Selected` limits additional connections by representation,
and `None` retains identity prerequisites only. These lower to
work in the same queue; they do not launch another whole-expression pipeline.
Representation filters select the `LibraryRep` base family, including duals;
for example, a Python `Representation.mink(4)` filter does not restrict permission
to four-dimensional slots. Actual connections must still satisfy dimension,
variance and positional port compatibility.
Python exposes both operations through explicit keyword-only parameters on
`TensorExpression`. For example,
`expr.contract(representations=[lorentz], collect_traces=False)` and
`expr.simplify_algebra(gamma=True, gamma_ordering="canonical", color=True)`.
The binding lowers those parameters into the existing Rust settings; it owns no
planning or contraction logic. Family-specific options require the corresponding
family to be enabled, and omitted options retain native defaults. The Python budget
is named `max_steps_per_domain` and defaults to `None` (unlimited); reusable configurations are ordinary dictionaries
passed with `**options`. Python settings classes are removed.

Python's `expr.contract(expand=False)` sets the same distribution permission as
Rust's `ContractSettings.expand=false`. Minimal contraction still substitutes boundary
ports through existing sums and reduces their branches. A connection that needs
a cross-product of independent sum alternatives remains exactly factored and is
outside that request's capabilities, so it need not cause `Deferred`. Enabled
identities retain their necessary prerequisites even under
`simplify_algebra(contract="minimal")`; this option does not forbid identity-generated
sums. The default permits necessary local distribution, and global arithmetic
expansion remains a separate request.

`AlgebraSettings` does not contain contraction settings. Algebra identities authorize their
required structural dependencies, not unrelated identity families.

The explicit notation methods are `to_dots`, `undo_dots`, `undo_chain`, and
`undo_trace`. The audited decision is to retain `to_dots`: supported index-free
metric notation `g(p(R),q(R))` can survive checked construction, and this method
changes its head to canonical `dot` notation. It does not discover repeated
indices. Direct typed multiplication can already normalize a scalar product
intrinsically; contraction filters cannot undo that intrinsic normalization. The schematic
nested spelling `p(q(R))` is rejected by current typed admission and is not an
additional supported notation. Undo methods unfold currently present notation
with fresh scoped indices and do not evaluate identities or traces.

Idenso's superseded high-level `simplify`, family-specific reduction wrappers,
`SimplifySettings`, `ContractionSettings`, public `chainify`, and aggregate
`undo_all` are removed. Family settings retain their existing names and kernels.
The unrelated orthogonal index-cooking settings and APIs remain available.
Production callers request their formerly enabled identities explicitly, with
an explicit ambient structural contraction where the old family wrapper also
performed that work. Former metric-only callers disable chain and trace
collection so migration does not silently enable new transformations.

Both public reductions lower into one domain queue. The expression retains
regional head multiplicities and boundary observations from the existing
candidate scanner. Actual replacements invalidate the affected observations. Cached shallow
structural graphs and unchanged interfaces are reused. The existing contraction
implementation handles structural work requested directly or needed by an
identity. Family kernels retain their specialized algebra. Neither reduction
introduces routine global distribution; explicit expansion remains separate.
The validation below checks both represented tensors and the reuse of these
shared planning and contraction paths.

== Invariant ownership and caller migration

`SymbolicTensor` retains the expression, ordered partial interface, typed-zero
metadata and trusted validation facts. Reduction and collection return that
same typed carrier without creating alias definitions. Spenso's ordered structures, slot matcher and fast
syntactic inference own representation compatibility, variance and logical slot
order. The shared depth-one parser builds planning boundaries and reuses those
interfaces; it is not a second tensor carrier. `SimplificationCandidates` supplies
regional head counts and index observations in its existing traversal. The
single Idenso domain queue owns eligibility, dependency invalidation, budgets
and completion. The existing slot contractor owns metric/vector substitutions;
shared collection primitives own dots, oriented chains and unevaluated traces.
Family kernels own gamma, colour and epsilon identities. Explicit arithmetic
expansion and Python object wrapping remain separate boundaries.

The migration removes repeated whole-family scheduling, global epsilon
permission inferred from registry presence, whole-alias-table progress comparison,
and contraction hidden inside a notation conversion. Trusted continuations pass
observations and interfaces forward; callbacks that can change rank invalidate
the affected proof. Initial domain scans and replacement-region observations are
distinct: the latter are required for new syntax and do not rescan unchanged
regions. Fresh-name discovery also reuses the observed written indices,
including indices inside opaque metadata.

The following caller inventory records permissions and required output, rather
than assuming that every historical method name meant the same work. Here
“ambient contraction” means unrestricted metric/vector contraction with
`collect_chains=false` and `collect_traces=false`.

#table(
  columns: (1.3fr, 2.3fr, 1.7fr),
  [Caller], [Enabled work], [Required output],
  [FeynKit amplitude conjugation],
  [Gamma settings retain `gamma0=true` and trace evaluation disabled; epsilon follows Reduced output; ambient contraction],
  [Resolved Atom with declared physical-leg port order],
  [FeynKit spin and colour completeness checks],
  [Spin identities use the previous gamma conventions plus epsilon and ambient contraction; colour completeness uses metric contraction with rank-one substitution disabled],
  [Exact spin projectors or representation dimension; no colour generator identity from a delta alone],
  [FeynKit TensorReducer],
  [Metric contraction with rank-one, chain and trace collection disabled, then pure dot notation conversion and Lorentz coefficient collection],
  [Selected Lorentz rows; scalar and foreign-representation coefficients stay factored],
  [GammaLoop numerator stages],
  [Colour-only or gamma-plus-epsilon according to the existing stage; ambient contraction after either],
  [Resolved and reversibly uncooked numerator Atom; stage order retained],
  [GammaLoop evaluator preprocessing],
  [When `do_algebra`: gamma, colour with fundamental-dimension invariants, epsilon, ambient contraction; otherwise no reduction],
  [Plain typed tensor, reversibly uncooked at the existing component-network boundary],
  [UV analytic reduction before Vakint],
  [Canonical gamma plus epsilon, ambient contraction and pure dot conversion],
  [Resolved uncooked Atom with compact scalar products; unsupported sigma input remains rejected],
  [UV cleanup after Vakint],
  [Metric-only contraction; then only Bispinor chain/trace collection, with metrics and rank-one substitutions disabled],
  [Collected unevaluated spin notation before the existing dummy-index cleanup],
  [UV comparison and renormalization reporting],
  [Colour-only algebra, ambient contraction and dot notation conversion],
  [Resolved uncooked expression for the existing factor/index comparison or parameter substitution],
  [GammaLoop Python `contract(str)`],
  [Structural metric/vector contraction with chain/trace collection disabled, then pure dot conversion],
  [Printed Atom string; replaces the contraction-performing old `to_dots` endpoint],
  [Native reference/FORM examples],
  [Explicit gamma or epsilon settings and matching ambient contractions; notation conversion is a separately named phase],
  [Each example's original resolved, factored or explicitly expanded result],
  [Notebook typed tensor routes],
  [Explicit gamma/colour/epsilon settings per case; matching ambient contractions; no new global expansion],
  [Same requested scalar or tensor expression; designated cases explicitly expand, with unrelated scalar spectators retained],
)

The caller table above records the migration scope. Historical disabled
`benches_old` and pinned reproducer sources are separate from production callers. GammaLoop's
unrelated numerical network `ContractionSettings` is not the removed Idenso
orchestration type and remains unchanged.

== Current plain-tensor migration

The tensor alias owner, definition registry, materializer, Python wrapper and
alias display controls have been removed. Reduction and collection return the
same typed carrier with an Atom payload, ordered interfaces, regional candidate
observations and cached completed requests. `Complete`, `Deferred` and `Capped`
remain distinct; unfinished results retain their exact expression.

Both Python operations default to `max_steps_per_domain=None`; Rust uses
`max_passes: None`. An explicit limit remains available and returns an exact
`Capped` result when exhausted. Automatic term-count, workspace and recursion-depth
cutoffs have been removed from the shared arithmetic tape and reduction paths.
Checked machine-arithmetic failures and repeated callback states retain exact
`Deferred` results. A separate memoization-size bound may decline a cache entry;
it never stops computation.

The performance audit uses the original generated four-loop numerators with
fully resolved structured indices, their complete scalar prefactors, symbolic
Lorentz dimension and spinor trace unit four. There is no preliminary contraction,
dummy relabelling or global numerator expansion. The historical source fingerprints
begin `63aa5a69` (quark) and `00ade8fe` (gluon). Reproduce ladder inputs and
measurements using `examples/notebooks/fermion_ladder.py` and the protocols in
`examples/notebooks/tensor_benchmarks.typ`. The exact edge-momentum
quark oracle has 525 terms; subsequent notebook routing gives 1,386. The generated
gluon oracle has 6,001 terms. These inputs do not exercise AUTO-index inference.

The final Python release `afa62f3b` completes both inputs in one call, including
with the unlimited default. Every run matches the exact oracle, preserves the
readmitted interface and has a stable completed rerun. Timings below are seconds
from an admitted, projected numerator; explicit expansion is the requested final
output operation.

#table(
  columns: (1fr, 1.6fr, 1fr, 1fr, 1fr),
  [Input], [Build], [Reduction], [Expansion], [Total],
  [Quark], [Archived `344a7942`], [31.201], [0.593], [31.794],
  [Quark], [Baseline `87425e2b`], [10.904], [0.043], [10.947],
  [Quark], [Final `afa62f3b`], [0.1069], [0.0468], [0.1537],
  [Gluon], [Archived `344a7942`], [0.199], [7.959], [8.158],
  [Gluon], [Baseline `87425e2b`], [14.585], [0.869], [15.453],
  [Gluon], [Final `afa62f3b`], [0.2473], [0.0072], [0.2545],
)

Thus the archived gluon's 0.2 s alias result was not its fully materialized result.
The old plain gluon required an exact continuation after an estimated frontier
bound; its 64 MiB estimate was not a measured resident-memory limit. That
15.45 s complete-output regression motivated the correction. Final timings are
means of two warmed samples following a first run; a fourth run checks the
unlimited default. Archived timings use one sample. Final completed reruns took
at most 0.13 ms. Input admission plus dimension assignment took 1.37 ms for
quark and 3.96 ms for gluon. Independent result readmission checks, roughly
53 ms and 75 ms, are outside computation timings. Peak process RSS was
71,244 KiB and 74,456 KiB, including all four samples, setup, oracles and checks.

A controlled quark comparison isolated an unnecessary index-label shape check:
original labels took 11.028 s, or 11.021 s after preliminary contraction; atomic
labels took 0.104 s, or 0.096 s after contraction. All four routes matched the same
525-term result, and their preprocessing costs were recorded separately. The
kernel now treats every admitted index label as an opaque identity. The
intermediate `5b8c459e` build already reduced the original source in 0.116 s and
expanded it in 0.048 s, without a notebook workaround.

The contraction engine now collects generated scalar coefficients with the
existing polynomial implementation while retaining unrelated sums as opaque
leaves. It passes scalar-closure facts forward instead of rediscovering them in
the emitted expression. A native run after removing automatic limits completed
the original gluon in one call: 0.357 s, with a stable 44.9 microsecond rerun.
Its 906,202-byte Atom prints as 2,858,500 characters and matches the exact oracle.
The separate instrumented call took 349.9 ms, including 279.4 ms contraction,
168.6 ms coefficient merging and 16.6 ms emission; these nested measurements
overlap. It recorded 106 candidate-scanner calls taking 0.870 ms. The maximum
estimated live workspace was 8,486,320 bytes, not a resident-memory measurement.
This native build uses optimization level two for Idenso and Spenso; it is
reported separately from Python release-module measurements.

Whole notebook time also includes rendering and additional calculations. With
the final installed module, quark `app.run()` took 7.90 s: diagram rendering
accounted for 6.64 s, while its two `simplify_algebra` calls together took
0.142 s. The gluon app took 2.92 s, including 2.39 s across three algebra calls:
it retains the additional full-numerator calculation as well as the projected
benchmark route. The notebook results contain 1,386 and 6,001 terms, with colour
factors $-1/54$ and $81$. Compilation and process startup are excluded from the
operation timings above; the final release build took 313.76 s. Runs used an
unpinned shared host; sample counts and measurement boundaries are preserved in
the dataset.
The separate full rank-four gluon reduction takes about 2.1 s and completes
without a limit. Its projection matches exactly $81$ times the 6,001-term
spacetime oracle when `color_substitute_cof_dimension_invariants=True` supplies
the oracle's colour convention; logical port order and Lorentz dimensions are
preserved. The timing probe retains symbolic colour invariants, so the convention
matched projection check is recorded separately. These notebook timings and
checks precede the concurrent community change `nqpkrrvpulko`, which adds typed
projector construction and display cells; they do not validate those later edits.
No new FORM timing was performed for these generated source fingerprints; the
historical comparisons below use different authored reference inputs.

Current validation passed 803 ordinary Idenso tests, seven optimized axial and
reference checks, 125 Rust binding tests, six shared-tape tests, 103 installed
Python tests, eight Higgs/generated-ladder tests, 75 selected notebook/API checks,
both ladder apps and four installed production examples. The optimized suite
includes the axial oracle that exceeded the ordinary runner's timeout. Strict
Clippy, formatting, typing, generated-stub checks, this report's Typst compilation
and dataset validation passed. The full documentation builder stops on an
unrelated pre-existing provenance mismatch for `gammaloop-api/src/python.rs`;
the working file and parent blob are identical, and no owner attestation was
changed. The dataset records that limitation and all test-log paths.

The preceding plain-carrier migration's independent tensor-component and
evaluator regression runs remain separate from these current checks. Python
defaults enable gamma and colour; callers that isolate a family explicitly
disable the other.

== Source and environment

The remaining sections record the historical alias-enabled implementation and
its validation runs. Their source paths, phase measurements and completion
claims refer to the named snapshots preserved on `tensor-aliases-preserved`,
rather than the current plain-tensor implementation described above.

The initial Git revision was `0f5a7d89789fa61b439808751bd3a7654c797cb0`, with
preexisting user changes. Their full binary-diff SHA-256 is
`aa647a48d4219a0625b20d92b0a63e09b894b2142c0566df45c04215262f242b`.
The frozen working-tree snapshot is the jj parent `10a6c180`. A separate
source copy was frozen before implementing the new API and built without
implementation edits. The original capture used `/tmp/unified-tensor-baseline`
on the capture host for executable copies, snapshots and logs. These temporary
artifacts are not required repository inputs; use the runner below to create a
new measurement record for a checked-out revision.

The host was `itphlies`, Linux 6.18.45, x86-64, AMD EPYC 9754, with 384
logical CPUs exposed and approximately 1.1 TiB memory. Measurement processes
inherited CPU affinity 0–383. The shared host was not reserved. Compilation
used four jobs pinned to CPUs 20–23; computation was not pinned. Rust and
Cargo were 1.98.1, Symbolica and Numerica 3.0.1, and Graphica 3.0.0. The
profile was `dev-optim`: workspace optimization level 2, dependencies and
Symbolica level 3, no LTO, incremental compilation disabled. A new comparison
should record its revision and dependency lockfile alongside the measurements.

The original example build took 401.684 s; the subsequent FORM example build
took 99.252 s after source-specific Cargo fingerprints were refreshed. These
durations are excluded from all computation timings. Preexisting Sep-27
executables named a different checkout in their dependency files. They were
retained as historical diagnostics and were not used as the source-certified
baseline.

== Reproduction and clock boundaries

The committed runner only orchestrates the existing native examples. It adds
no tensor parser, algebra kernel, or contraction implementation. Build first
in the pinned development shell, then run with a fresh output directory:

```sh
nix develop --command cargo build --locked -p idenso \
  --profile dev-optim --features reference-cases \
  --example metric_contraction_benchmark \
  --example contraction_phase_benchmark --example trace_form

nix develop --command python crates/idenso/benches/unified_reduction.py \
  --source . --binaries target/dev-optim/examples \
  --form /nix/store/c2x9lj2l3x99xrdx1ba1a63lxyxq6ji0-form-5.0.0/bin/form \
  --samples 3 --target-ms 8 --form-max-length 12 \
  --output /tmp/unified-tensor-candidate
```

For the original native build the pinned GCC wrapper was
`/nix/store/z4c6k0mrlkwl3s4w9ysxc8vq1wylm3ms-gcc-wrapper-15.3.0/bin/gcc`,
passed through `RUSTFLAGS=-C linker=...`; the pinned Rust toolchain was
`/nix/store/c3bpv60f5wgp4pmgsifaykxk214h6463-rust-stable-with-components-2026-09-03/bin`.
These paths and binary hashes identify the captured environment. A portable
reproduction may supply its own FORM executable and must record the different
toolchain identity. Compare the same source inputs, permitted transformations,
conventions, and final output requirements on both sides.

The main native suite has three retained samples per case/method and an 8 ms
batch calibration target. Timed typed operations include checked admission,
the requested operation, alias resolution, and result destruction. Input
construction, checking, and snapshot formatting are outside those clocks.
Phase runs are separate executions; overlapping phases are not summed into
an end-to-end estimate. The runner also records wider process clocks including
startup and reporting. Cumulative child maximum RSS is labeled as such and
is not an isolated per-operation allocation measurement.

Inputs and outputs stay factorized unless the original example explicitly
requests expansion. Neither the runner nor the comparison silently converts a
factored result into a fully materialized polynomial. The JSON dataset retains
all native sample times, output term counts, source and input/output hashes,
refused cases, and FORM observations.

== Baseline native measurements

The main suite completed with 432 case/method rows: 91 dot normalization,
71 parsing, 71 parsing plus dot normalization, 91 repeated-index observation,
44 typed contraction, 44 typed epsilon, and 20 complete typed gamma rows.
94 method/case refusals are retained, predominantly untagged generic tensors
in historical fixtures rejected by the current strict admission boundary.

#table(
  columns: (2.5fr, 1fr),
  [Complete typed gamma case], [Median ms],
  [Free ordinary trace, length 8], [6.868],
  [Free ordinary trace, length 12], [295.899],
  [Axial trace, length 6], [1.413],
  [Axial trace, length 8], [6.863],
  [Axial trace, length 10], [118.460],
  [Axial trace, length 12], [831.539],
  [Late-metric ordinary trace, length 10, D dimensions], [131.469],
)

For the closed 32-metric loop, parsing alone took 0.081 ms and complete typed
contraction with admission and alias resolution took 1.408 ms. For free8,
contraction of the already evaluated output took 4.409 ms; free12 took
268.074 ms. These are diagnostic calls on a different input stage and must
not be subtracted from the original trace timings to infer an isolated kernel
cost. All raw rows are retained for later phase comparisons.

The combined phase harness failed its baseline rerun-stability assertion on
`axial_even4_10`: the factorized symbolic expression changes on another full
gamma call. Both forms and the failure are retained; this is not asserted to
be a component-value discrepancy. A separate source-only phase attempt
refused a compact gamma signature. Independent attempts for every original
source retained 19 successes and one refusal (`compact_paired_8`), producing
627 phase rows. Failures were not discarded or replaced by successful retries.

== FORM comparisons

Fresh FORM 5.0.0 comparisons passed for free, paired and alternating words at
each even length 2 through 12 in both `trace4` and `tracen` modes. They use the
same four-dimensional metric convention and Tr(1)=4. Both tools were checked
against a separate exact scalar recurrence after the example's stated vector
projection. This certifies those projected polynomials, not every free-index
tensor component. The native trace example retains all generated FORM programs
and output. Full HEP component tests remain a separate correctness gate.

#table(
  columns: (1.3fr, 1fr, 1fr, 1fr),
  [Free word], [Idenso warm ms], [FORM trace4 process ms], [FORM tracen process ms],
  [Length 4], [0.370], [8.951], [9.339],
  [Length 8], [6.983], [10.097], [9.077],
  [Length 12], [300.534], [15.671], [22.639],
)

Idenso warm clocks are medians of five calls after a separately recorded cold
call. FORM clocks are medians of five fresh processes after one warmup and
include startup, parsing, sorting, and scalar verification. The different
boundaries do not establish kernel-only speed ratios. Free12 has 4383 terms
in the Idenso and FORM trace4 forms and 10395 in FORM tracen before projection.
All 18 rows, exact-check receipts and executable hashes are in the dataset.

== Complete notebook baseline and independent checks

A new combined FeynKit/tensor Python host was built from the same frozen Rust
sources, using Python 3.14.7, Symbolica 3.0.1, Linnet 0.1.0 and Typst 0.15.0.
The nested host's stale lockfile was replaced with the frozen workspace lockfile
and independently resolved offline. No Rust source changed. The resulting host
lockfile hash is `6cb1e9b77fb740e7e036db888a7222f2af4cf500d20eaa04e1797e4fb8090e6e`;
the release core hash is
`90ee7075ddbe537dcd64ab499f342d4ded1a0f7ac11b03921a4e54b5216cf880`.
Release compilation took 514.255 s, excluded from computation. Build processes
used CPUs 24–27. Computation used CPU 28 on the shared host.

The frozen notebook cases needed only an import-namespace correction from the
historical `symbolica.community.hep` to the host's registered
`symbolica.community.feynkit`. The historical cohort recorded source and host
identities, three samples per case, phase diagnostics, output term counts,
process peak RSS and independent checks. The notebook harness and clock
boundaries below describe how to capture a new cohort.
The original 16-case suite completed without failures. A new closed SU(3)
colour case uses the same scheduler and output protocol.

The first numeric column below is a continuous wall clock from fresh case/input
construction through reduction and the requested materialization. Interpreter
startup, imports, checking and output disposal are outside this clock. The
second is process CPU time for prepared-input reduction including requested
materialization. FORM body CPU time uses independent copies and excludes
startup, parsing, initial construction and final validation. Process wall clocks
are retained separately. These boundaries are explicit and are not interchangeable.

#table(
  columns: (2fr, 1fr, 1fr, 1fr),
  [Case], [Fresh complete ms], [Prepared CPU ms], [FORM body CPU ms],
  [free4-trace4], [0.498], [0.249], [0.002],
  [free4-tracen], [0.520], [0.260], [0.001],
  [free12-trace4], [80.112], [97.009], [3.480],
  [free12-tracen], [66.169], [70.771], [5.176],
  [free14-trace4], [324.119], [341.764], [19.900],
  [free14-tracen], [625.842], [620.847], [57.500],
  [axial12-trace4], [790.991], [813.718], [0.735],
  [fermion-3-trace4], [3.685], [1.645], [0.330],
  [fermion-3-tracen], [10.817], [8.346], [1.490],
  [fermion-4-trace4], [17.191], [13.304], [5.783],
  [fermion-4-tracen], [119.110], [115.349], [49.500],
  [gluon-4-trace4], [1629.228], [1573.364], [670.000],
  [gluon-4-tracen], [1693.773], [1626.695], [1028.000],
  [historical-original], [2278.131], [2265.774], [726.000],
  [historical-early], [310.199], [310.496], [165.500],
  [production-aa-aa], [1742.147], [1720.093], [unavailable],
  [colour-fierz], [3.740], [3.597], [0.004],

)

The four-dimensional free12/free14 traces and axial12 have different symbolic
bases in FORM and Idenso. Their checks compare the source Dirac network,
Idenso output and FORM output at three exact HEP vector-component points with
nonzero spanning minors. This is finite exact evidence, not a general symbolic
Schouten certificate. Ordinary tracen polynomials match FORM exactly; where the
D-to-four forms differ, the same exact four-dimensional component check is used.
Fermionic and gluonic ladder scalar polynomials match FORM exactly in both
trace4 and tracen, with D-to-four checks and independent HEP component checks.
The historical ladder also agrees with its separate Symbolica scalar recipe.

The captured mixed colour/Lorentz production numerator has independent exact
HEP component checks at the fixture's documented three momentum assignments,
including its symbolic dimension coefficients. It has no equivalent FORM model
in the existing suite; no FORM time is claimed for that row. This is a fixture
certificate, not a general D-dimensional Clifford or arbitrary-momentum proof.

The new colour case is the fully closed SU(3) word `Tr(Ta Tb Ta Tb)` with summed
adjoint indices and generator normalization `Tr(Ta Tb)=delta(a,b)/2`. Its complete
HEP generator-component contraction, Idenso reduction, and FORM's explicit
Fierz-plus-delta calculation all give `-2/3`. A first fixture-adapter refusal and
a FORM sample below timer resolution were retained. The existing fast-case
200 ms calibration protocol was explicitly extended to this case before its
qualified three-sample cohort. The successful colour rerun of the already
reduced alias owner took a median 0.117 ms in the frozen implementation.

The additional mixed colour/Lorentz case multiplies that closed SU(3) word by
an open four-gamma trace. The legacy implementation applies its colour and gamma
entry points; the candidate requests both families in one planner. Both request
the same fully materialized four-port polynomial. Frozen Idenso and FORM agree
exactly, and all 256 external components agree with the complete original HEP
component network. The frozen median fresh end-to-end time is 4.983 ms,
prepared reduction CPU time is 4.445 ms, and FORM body CPU time is 0.009796 ms.
The fast-case FORM calibration applies to this small case as well. This supplies
a mixed-family FORM comparison; it does not supply a FORM model for the separate
production fixture.

Separate completed-gamma reruns preserve the retained owner and exclude fresh
admission, materialization and checks. The median wall times were free4 trace4
0.025 ms, free12 trace4 11.065 ms, free12 tracen 5.867 ms, axial12 194.557 ms,
fermion3 trace4 0.143 ms and fermion3 tracen 1.286 ms. All six results retained
the exact resolved expression on each of three samples and reported completion.
The separate native axial10 baseline instability remains as recorded above.
The executable Python protocol and sample values are embedded in the dataset.

The axial case retains `(x+y)^8` as one spectator factor and materializes only
the trace. The ordinary traces and ladders request their documented polynomial
outputs explicitly. Raw phase records distinguish gamma/colour kernels,
contraction, alias materialization, scalar momentum routing and normalization;
the independent clocks are never summed to approximate end-to-end time.
Process peak RSS includes initialization, warmup and diagnostics, and is not an
allocation measurement of one kernel.

The existing notebook scheduler can reproduce the cohort with an explicitly
qualified release interpreter and fresh output directory:

```sh
python examples/notebooks/fermion_ladder.py --suite consolidation \
  --interpreter baseline=/path/to/frozen/venv/bin/python \
  --core baseline=90ee7075ddbe537dcd64ab499f342d4ded1a0f7ac11b03921a4e54b5216cf880 \
  --complete-algebra --gate-stage baseline --cpu 28 \
  --form /path/to/form --output /tmp/frozen-notebook-cohort.json
```

Use the preserved frozen frontend for the old API, and the migrated frontend
for the candidate. The scheduler supports a source-hashed `--frontend` manifest
when comparing them in one paired run. The two frontend versions must retain
the same inputs, permitted operations and requested outputs.

== Shared-alias baseline

The separate `--ladder-route baseline=factorized` cohort passed the same exact
FORM and component checks on CPU 29. It retains shared aliases during reduction
and explicitly materializes the same final polynomial. Its timings therefore
include final materialization; an unexpanded alias result is not substituted
for the requested output.

#table(
  columns: (2fr, 1fr, 1fr, 1fr, 1fr),
  [Case], [Fresh complete ms], [Contraction CPU ms], [Materialization CPU ms], [Definitions],
  [historical-original], [92.980], [33.899], [52.068], [163],
  [historical-early], [104.189], [34.723], [53.843], [163],
  [gluon-4-trace4], [271.339], [178.803], [82.429], [488],
  [gluon-4-tracen], [356.439], [186.689], [158.724], [520],

)

Rule application plus admission took 5.376–5.807 ms in these separate diagnostic
calls. The historical words retained 86689 aliased bytes; gluon trace4 retained
277347 bytes and gluon tracen retained 309124 bytes before explicit expansion.
These observations reuse the existing benchmark instrumentation. Materialization
is a substantial baseline cost and must remain in equivalent-output comparisons.

== Initial source-qualified release comparison

This section preserves the initial unified implementation's measurements before
the regression correction below. Its build identity, samples and independent
checks remain unchanged; “final” in these historical tables identifies that
initial implementation capture.

The final release core is
`83079179ede4a59e4231eca6cae9a3b0396fa16999dd875c87bcb77868b0fa6f`.
It includes the final per-term observation proof and the explicit migrated caller
policies. Compilation took 153.977 s, excluded from computation. The manifest
changed only by deletion of an unused `AtomCore` import inside a `cfg(test)`
module. Reversing that one deletion reproduced the exact pre-build SHA-256;
that module is absent from the release host. The dataset records both source
manifests and this qualification instead of claiming identical full source trees.
Subsequent documentation-test fixtures and GammaLoop command Clippy cleanups
are outside the measured runtime dependencies; their exact file hashes are also
recorded. The measured tensor runtime has no subsequent source drift.

The independently checked preceding host completed all 18 typed cases and four
shared-alias cases, with 38 and 10 check records respectively. Those checks
include the exact FORM and HEP qualifications described above. For the final
host, `crates/idenso/benches/unified_replay.py` replays the existing source-bound
child commands, including their warmups, phase diagnostics, completion assertions
and three rotating baseline/candidate rounds. All 114 typed and 30 alias records
match both the independently checked input hash and serialized result exactly.
This transfers the independent certificate by exact equality; it is not a new
component comparison implemented through the modified tensor intake. The replay
preserves all new timings and reruns the frozen baseline. FORM clocks retain the
preceding cohort's timestamp and identical inputs rather than being presented as
new simultaneous measurements.

Fresh end-to-end wall clocks below include construction, admission, reduction
and the requested materialization. Prepared CPU clocks use the already prepared
input and include the same requested final output. FORM body CPU clocks exclude
startup and construction; these columns have different boundaries and do not
justify kernel-only speed ratios. Computation uses CPU 28; alias comparisons
use CPU 29. The shared host was not reserved. The reproduction runners capture
samples, process clocks, phase observations, final term counts and process peak
RSS in a caller-selected output directory.

#table(
  columns: (1.8fr, 1fr, 1fr, 1fr, 1fr),
  [Case], [Baseline fresh ms], [Final fresh ms], [Final prepared CPU ms], [FORM body CPU ms],
  [mixed-colour-trace4], [4.772], [6.225], [5.528], [0.00952],
  [colour-fierz], [3.656], [5.575], [4.324], [0.00354],
  [fermion-3-trace4], [3.246], [3.253], [1.401], [0.32500],
  [fermion-3-tracen], [10.780], [15.475], [13.025], [1.50000],
  [fermion-4-trace4], [17.047], [23.663], [19.220], [5.64706],
  [fermion-4-tracen], [110.759], [622.236], [607.246], [47.20000],
  [gluon-4-trace4], [1386.476], [1437.495], [1471.643], [601.00000],
  [gluon-4-tracen], [1663.714], [1670.387], [1593.990], [870.00000],
  [historical-original], [2310.062], [2301.286], [2285.803], [814.00000],
  [historical-early], [313.307], [313.384], [314.220], [158.00000],
  [free4-trace4], [0.493], [0.513], [0.245], [0.00177],
  [free4-tracen], [0.512], [0.518], [0.266], [0.00135],
  [free12-trace4], [97.113], [382.159], [380.301], [3.13433],
  [free12-tracen], [76.632], [147.661], [120.684], [4.05970],
  [free14-trace4], [355.821], [6339.111], [6815.558], [20.00000],
  [free14-tracen], [612.940], [1079.460], [1078.885], [58.00000],
  [axial12-trace4], [731.440], [117.017], [115.604], [0.65000],
  [production-aa-aa], [2017.386], [3110.005], [3170.475], [unavailable],
)

The axial12 workload improves from 731.440 to 117.017 ms. Regressions remain:
free14 trace4 increases from 355.821 to 6339.111 ms (17.8 times), and the
four-loop fermionic tracen case from 110.759 to 622.236 ms (5.6 times). These
are completed operations with equivalent requested outputs and retained samples.
No speedup threshold is used to accept or exclude a result. The short traces,
historical ladders and shared-alias routes have much smaller changes on this
shared host. Three retained samples do not establish statistical significance
for small percentage differences. The 18-case replay takes 650.9998 s including process startup,
warmups, diagnostic work and serialization; that duration is not an algebra time.

#table(
  columns: (1.8fr, 1fr, 1fr, 1fr, 1fr),
  [Alias route], [Baseline fresh ms], [Final fresh ms], [Final contraction CPU ms], [Final materialization CPU ms],
  [historical-original], [94.398], [93.374], [37.083], [46.388],
  [historical-early], [92.379], [90.060], [36.133], [42.078],
  [gluon-4-trace4], [243.618], [241.503], [160.569], [67.134],
  [gluon-4-tracen], [315.500], [318.326], [167.602], [134.326],
 )

The alias routes retain 163, 488 and 520 definitions for the historical, gluon
trace4 and gluon tracen cases, respectively, before explicitly materializing the
same output. Their intermediate alias weights and final term counts are recorded
by depth. The alias replay takes 101.1935 s including startup and diagnostics.

== Initial completed-owner reruns

These separate CPU-32 measurements retain the completed owner and exclude fresh
admission, result materialization and equality checks. Every sample keeps the
exact resolved expression. The initial implementation reused a same-operation
completion record. The final column instead repeats both operations with the old
gamma wrapper's combined permissions; its latest-policy-only record could not
serve both calls. These historical measurements expose the migrated sequence's
cost before the completion-record correction below.

#table(
  columns: (1.8fr, 1fr, 1fr, 1fr, 1fr),
  [Case], [Old gamma rerun ms], [Final algebra rerun ms], [Final contract rerun ms], [Final combined rerun ms],
  [free4-trace4], [0.023791], [0.002721], [0.001971], [0.043090],
  [free12-trace4], [10.741150], [0.008810], [0.009760], [325.133288],
  [free12-tracen], [6.079513], [0.008200], [0.009100], [68.035945],
  [axial12-trace4], [191.628663], [0.006100], [0.006410], [61.574821],
  [fermion-3-trace4], [0.136921], [0.001149], [0.001080], [0.187400],
  [fermion-3-tracen], [1.179254], [0.001580], [0.001349], [6.956885],
 )

The cached operations remained stable and took approximately 1–9 microseconds in
these samples. Alternating permissions incurred substantial planning work on
large outputs, with exact completed results preserved. The follow-up measurements
below repeat both policies after correcting this limitation.

== Initial native phase and interface measurements

The final `dev-optim` feature build took 1994.294 s including shared Cargo target
lock waits; compilation is excluded from the following clocks. Its only observed
source drift was the same independently verified test-only import deletion as
the release host. The final native suite ran on CPU 30 with the same frozen
baseline binaries remeasured on CPU 30. These native calls preserve the examples'
resolved/factorized output requirements and must not be compared directly with
the release Python polynomial timings above.

There are 432 native case/method records and 627 source-phase records. All 590
public planner outcomes captured during those diagnostics are `Complete` under
their recorded policies. The final per-term observation proof removes the false
`Deferred` result previously produced by free ports repeated across alternative
sum branches. All 18 native projected-polynomial FORM comparisons pass. Nineteen
source phase attempts pass; the historical compact gamma signature remains a
checked-admission refusal. The full phase harness gets past the old axial10
stability failure and then refuses the historical `chain_partner` fixture's
out-of-scope placeholders. Those refusals and 94 main-suite historical method
refusals are preserved; they are not completed operation samples.

#table(
  columns: (1.9fr, 1fr, 1fr, 1fr),
  [Native case], [Baseline median ms], [Final median ms], [Separate diagnostic ms],
  [metric_loop_32], [1.373], [1.314], [1.372],
  [ordinary_even4_10], [9.778], [16.220], [16.278],
  [free_12], [303.158], [626.753], [631.222],
  [axial_12], [837.451], [165.757], [168.582],
  [late_metric_4_12_axial], [226.195], [190.941], [193.792],
 )

The following wall times come from each separate diagnostic call. They are the
outermost time for that category, with different categories still overlapping.
For example, planning contains contraction and kernels, and materialization can
contain interface discovery. They are not additive components of a stacked
end-to-end time.

#table(
  columns: (1.5fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr),
  [Case], [Admission ms], [Scan ms], [Planning ms], [Gamma ms], [Contraction ms], [Materialization ms],
  [metric_loop_32], [0.429], [0.152], [0.759], [0.000], [0.748], [0.001],
  [ordinary_even4_10], [0.323], [2.960], [13.633], [1.262], [8.097], [1.939],
  [free_12], [0.458], [20.113], [382.217], [11.857], [330.074], [243.433],
  [axial_12], [0.249], [8.189], [113.175], [2.943], [25.306], [53.588],
  [late_metric_4_12_axial], [54.903], [26.901], [122.095], [0.000], [122.010], [0.170],
 )

In the free12 diagnostic, planning takes 382.217 ms and structural contraction
330.074 ms across 994 scopes; gamma kernel work is 11.857 ms. Alias materialization
takes 243.433 ms and overlapping interface discovery 233.877 ms across 76775
entries. This locates substantial costs outside the gamma kernel without
pretending the overlapping clocks identify a unique low-level cause. Axial12
uses 140 shallow graph builds, 59 fresh-index reservation scopes (21.453 ms),
and 126 epsilon-kernel scopes (0.081 ms). Initial scan counts include newly
created domains and changed-region observations; they are not evidence of
rescanning one unchanged definition.

A separate CPU-31 process enables the existing Spenso profiler. The table below
reports underlying interface work, not just shallow node counts. Cache event
counters cover classification and interface lookups; they are not the same
population as structure attempts. Miss-side query counts include queries made
without a reusable cache. Recursive visit counts cover ordinary syntax;
chain/trace dispatch has separate rules.

#table(
  columns: (1.6fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr),
  [Case/operation], [Structure queries], [Cache hits], [Cache misses], [Uncached queries], [Recursive visits], [Uncached ms],
  [metric_loop_32], [32], [0], [32], [32], [32], [0.029],
  [late_metric_D_10_ordinary], [3998], [3563], [435], [435], [2791], [3.302],
  [free_12], [0], [0], [0], [0], [1579], [0.000],
  [axial_12], [838], [44], [606], [794], [2592], [4.693],
  [late_metric_4_12_axial], [5772], [3604], [633], [2168], [6925], [9.752],
 )

Free12 performs 1579 ordinary inference visits even though this capture builds
no parser structure query. A shallow or absent graph therefore does not imply
zero underlying inference. The late-metric D10 case demonstrates substantial
boundary reuse while still recording uncached work. All 432 labeled profiler
captures are retained separately; their enabled profiler overhead never enters
the primary comparison.

The existing ten SU(3) reference cases also run through the shared `ReferenceCase`
helper under the feature clock. A small external Rust diagnostic driver links
the certified Idenso rlib; its complete source, compiler command, fingerprint and
hash are retained. It introduces no production helper or alternate algebra path.
Each case has one unmeasured warmup and three separately instrumented captures,
with exact output stability checks. All 60 public planner observations are
`Complete`. These clocks are diagnostic; they do not replace the independent
colour/Fierz FORM and HEP certificate above.

#table(
  columns: (1.9fr, 1fr, 1fr, 1fr, 1fr),
  [Existing colour reference case], [Diagnostic median ms], [Colour kernel ms], [Contraction ms], [Materialization ms],
  [`color_form_one_generator_trace`], [0.257], [0.145], [0.003], [0.002],
  [`color_form_two_generator_trace`], [0.752], [0.589], [0.024], [0.003],
  [`color_form_fierz_generator_contraction`], [3.615], [3.343], [0.731], [0.005],
  [`color_feyncalc_open_chain_separated_casimir_id2`], [2.000], [1.730], [0.310], [0.005],
  [`color_feyncalc_structure_loop_to_adjoint_delta_id3`], [0.762], [0.540], [0.107], [0.003],
  [`color_feyncalc_structure_times_open_chain_id4`], [1.732], [1.423], [0.304], [0.004],
  [`color_feyncalc_delta_closes_doubled_two_generator_chain_id18`], [1.589], [0.895], [0.474], [0.003],
  [`color_feyncalc_three_generator_trace_terminal_id30`], [2.343], [1.927], [0.223], [0.199],
  [`color_feyncalc_repeated_trace_pair_id52`], [1.877], [1.466], [0.327], [0.004],
  [`color_form_four_generator_trace_terminal`], [8.348], [6.973], [1.342], [0.582],
 )

Native cumulative child peak RSS is 569644 KiB. Candidate release child peak RSS
ranges from 122956 to 386096 KiB over the typed samples. These process peaks
include imports, warmup and diagnostics, not isolated allocation peaks of one
kernel. The largest requested final polynomial has 135135 terms. Input/final
term counts and alias-depth intermediate weights are recorded; uniform kernel
peak-term counts and isolated allocator peaks were not instrumented.

== Candidate phase-clock protocol

The existing `reference-cases` feature now offers an opt-in diagnostic capture
around the shared Rust implementation. The native metric and phase examples
emit one separately instrumented invocation after their primary timing samples;
primary samples leave the capture inactive, retaining only the feature build's
inactive capture checks. The release Python host does not enable this feature.
The instrumented categories are
admission, Idenso interface discovery, validation, candidate scanning, planning,
shallow graph construction, gamma/colour/epsilon kernels, structural contraction,
fresh-index reservation, alias materialization, and explicit output normalization.
The same capture records each public planner outcome as `Complete`, `Deferred`,
or `Capped`, including reuse of a previously completed result. These outcomes
qualify each recorded call under its requested policy. Notebook diagnostic calls
also retain the algebra and contraction owners long enough to assert `Complete`
before materialization hides their completion records. Exact output comparison
alone is not used as evidence of completion.

Each category records call count, inclusive wall nanoseconds, and wall time for
its outermost scopes. Outermost scopes avoid counting recursion of the same
category twice. Different categories still overlap: planning encloses kernels;
kernels can enclose contraction; graph construction encloses Spenso leaf
interface walks. These numbers must not be summed or subtracted from separately
measured medians. They measure the capturing thread's elapsed wall time, not
worker-thread CPU totals. Automatic Symbolica function normalization remains
inside its caller's time; the output-normalization category covers the explicit
shared normalization entry points.

The explicit-conversion phase is now labeled `NotationToDots`. The frozen
`ToDots` phase used the old normalizer that could also normalize indexed powers;
the current public conversion only changes notation. No before/after ratio is
computed between those different operations. The separate `normalize_dots` raw
kernel diagnostic continues to exercise the shared default internal normalizer,
and complete-operation comparisons retain their matched structural permissions.

A separate diagnostic process can set `SPENSO_NETWORK_PROFILE=1` to reuse
Spenso's existing parser profiler. The capture resets that profiler immediately
before its diagnostic call and reports afterward with a case/method label.
`parse.structure` times all structure queries, including cache lookup and logical
slot reconstruction. Structure-attempt and opaque-interface cache-hit/cache-miss
counters identify reuse; they must not be relabeled as deep-inference visits.
The separate `parse.structure.uncached` clock encloses only miss-side interface
inference, excluding cache-hit slot reconstruction. `parse.structure.visit`
counts each recursive ordinary-syntax inference entry, including entries made
outside shallow graph construction. Chain and trace dispatch have their own
structure rules, so this counter is not a count of every possible tensor
inspection. These distinctions expose underlying inference work as well as
shallow graph reuse. This process-wide profiler sums nested scope times across
threads, so its totals can overlap and are not CPU time. Its report rounds total
times to microsecond resolution.
The primary performance cohorts leave the environment switch unset. No primary
time from a profiler-enabled process is substituted into the main comparison.

The frozen baseline predates these guards and has no equivalent internal phase
clock data. Its existing parse-only, prepared-kernel, materialization, and
end-to-end measurements retain their recorded boundaries. In particular,
parse-only Atom text parsing is not a measurement of full tensor admission.
The candidate clocks complement that baseline without inventing missing data.
All three feature-specific diagnostic tests pass in the final native suite:
exact contraction is preserved, nested scopes close, inactive capture stays
inactive, panic cleanup permits a subsequent capture, and recorded outcomes
distinguish cached completion from exact budget capping.

== Post-implementation validation

The final Idenso feature suite passes 837 tests with 22 ignored tests. The one
separately filtered exhaustive axial test also passes on the same final runtime,
covering every gamma5 position in its reference set in 136.722 s. The ordinary
suite takes 16.30 s after compilation. Final Spenso shadowing tests pass 335,
Rust bindings pass 130, and the final release host passes 87 installed API tests,
including public signatures, representation filters, notation round trips and
completion status. The binding unit run predates the final engine observation
correction; the final host and Idenso suite cover that correction explicitly.

The final production freshness run passes all 67 FeynKit tensor and amplitude
tests with one skipped test in 0.151 s after a separate 17m53s build and shared
target wait. Generator/amplitude integration passes eight tests in 0.139 s.
The grouped production Cargo check and the GammaLoop tests compilation also
pass. Native doctests pass the one enabled test, with 13 ignored. Documentation
examples pass 88 tests and catalogs pass 23. Isolated generation and check of
Spenso, FeynKit and GammaLoop stubs all pass; the Spenso and GammaLoop exporters
each pass three tests. Rust formatting, strict Clippy across the affected crates
and all targets, the docs builder, changed Python formatting and undefined-name
checks (84 files), and the final documentation compilation pass. The Clippy run
retains Cargo's dependency future-incompatibility notice for
`proc-macro-error2`; no warning-denial check failed. Build and lock-wait durations
are kept separate from test or algebra execution.

Correctness regressions cover unevaluated traces, colour deltas without colour
identities, empty and unrestricted representation filters, mixed-representation
ports, prerequisite contractions, disabled families and newly emitted candidates.
They preserve unrelated sums, scalar-interface internal contractions, independent
powered dummy scopes, callback-induced rank changes, logical port order, typed
zeros, metadata and alias associations. Component evaluators and exact identities
check represented tensors independently when symbolic forms differ. The
notation tests check fresh-index collisions, orientation, surviving explicit
forms and contraction round trips without trace evaluation. Exact capped results
retain the remaining work and aliases; completed reruns are stable.

Architectural counters and regressions check one initial observation per domain,
reuse of unchanged definitions and interfaces, regional head multiplicities,
forwarded candidates and alias dependency invalidation. Trusted continuations
avoid redundant admission; cached depth-one boundaries avoid opening unselected
large leaves. Changed-source continuation tests verify underlying inference
cache reuse, not only graph size. The shared domain queue uses actual changes
and pending regions, with no whole-alias-table equality test for progress.
Fresh-index reservation reuses regional written-index observations, including
opaque metadata; replacement removal updates those facts. Separate replacement
observations remain necessary for newly produced or uncontrolled syntax.

The synchronized tensor-registration regression passes with eight concurrent
callers across 32 new tensor/vector names; conflicting declarations remain
rejected. A separate test launches a fresh process and declares a vector before
any tensor or Symbolica state is warmed. It passes in 0.02 s after correcting
registration order: the existing lazy bootstrap runs before the registration
mutex, because initialization can recursively register other tensor heads.
The production suite independently covers this cold-start path. This corrects
the registration races and cold-start timeouts exposed during migration.

Final measured public calls report `Complete` under their recorded policies.
A deliberately rounded coefficient inside a selected sum,
`(Float(0.1)*p+q)*(r+s)`, remains an exact unchanged `Deferred` frontier; the
collector does not certify the required branch-local operation in that case.
The same coefficient outside the selected region permits completion. A syntactic
head or repeated free port alone is not eligibility: the final observer combines
per-term incidence using maximum across sums and addition across products, then
uses the established external interface to discharge those false frontiers.
No whole-expression rescan is added for that proof. Budget regressions separately
exercise `Capped`, including zero-budget requests, without discarding terms.

== Remaining costs and comparison limits

Planning, contraction and explicit materialization remain substantial costs.
The historical 17.8-times free14 trace4 regression and 5.6-times fermion4 tracen
regression remain in the initial measurement set. Investigation subsequently
identified avoidable per-definition reconstruction of the whole alias registry
and structural cleanup scheduled for definitions already known to have only
free ports. The correction and fresh matched measurements below distinguish
that implementation defect from unavoidable materialization work. Completed
owners now retain exact-policy completion facts across operations that make no
change; an actual root or definition change invalidates those facts. No policy
subsumption or completion of a capped/deferred request is assumed.

Fresh-index reservation no longer routinely walks every observed definition,
and trusted shallow graph rebuilding reuses exact cached leaf interfaces.
The recorded cache misses, recursive inference visits, interface-discovery
entries and materialization clocks show that these changes do not eliminate
all traversal. Initial observations include collision-index discovery; unknown
handles, changed regions and callbacks still require justified local checks.
The parser's standalone admission checks remain available for untrusted input.

The historical refused fixtures are not completed benchmark results. The frozen
baseline's axial rerun failure is retained, while the final suite passes its
independent axial checks. The FORM comparison limits remain explicit: native
projections alone are not full-tensor certificates; finite exact HEP component
assignments are not general symbolic proofs; production-aa-aa has no matched
FORM model here; dimensional axial identities are compared only where defined.
The frozen baseline has no retrospectively invented internal phase clocks, and
uniform intermediate kernel term peaks and isolated allocator peaks are not
available. All requested output forms are labeled as factorized, aliased or
explicitly materialized rather than compared as interchangeable timings.

The replay runner `crates/idenso/benches/unified_replay.py` and the native
examples listed above preserve the reproduction protocol without adding a
tensor implementation. Capture generated measurements and validation receipts
outside the source tree, retaining their revision, settings and clock boundaries
when making a new comparison.

== Regression correction and fresh comparison

The follow-up keeps the same public APIs, settings, kernels and requested
outputs. It removes two sources of avoidable work in the shared planner:

- Structural cleanup authorized by an identity is discharged using the domain's
  existing observations when they already prove no eligible connection remains.
  Free-port trace definitions therefore do not enter the contractor or consume
  a change budget merely because their parent kernel emitted them.
- Domain contraction borrows the planner's alias registry and returns new or
  changed associations. It no longer rebuilds and compares the full registry
  for every definition. The mutable registry representation required by selected
  collection is constructed only when that operation needs it. The planner's
  definition snapshot is lazy and invalidated by actual registry changes.
- Completed owners retain the exact policies proved complete for that immutable
  root and alias DAG. A no-op operation preserves earlier completion facts; any
  root or definition change clears them. A capped or deferred policy still runs
  again and is never turned into a completed fact by another operation.

The preserved initial unified release core is
`83079179ede4a59e4231eca6cae9a3b0396fa16999dd875c87bcb77868b0fa6f`.
The corrected core is
`84687996ae4eb20218bd7af34356c112dcea7b215f00c7de6a2fe9cb6545889c`.
Its separate release build took 162.584 s, excluded from every computation
clock. The measured Rust source hashes matched before and after that build and
again after capture. The original baseline and initial unified wheel,
interpreters, inputs and datasets remain unchanged. The corrected wheel uses a
separate interpreter at `/tmp/unified-regression-fix/fixed-venv/bin/python`.

The existing replay runner executes the same frozen benchmark children and
settings for seven typed cases and three shared-alias cases. Each corrected
input hash and raw serialized output equals the retained independently checked
cohort exactly. All 120 child records across the before/after comparisons pass
that check, including the newly sampled original baseline. There are three
samples of each unified version and six baseline samples, with the baseline
remeasured alongside each version. The table reports medians of complete fresh
operations, including their requested materialization; startup, compilation,
equality checks and serialization are excluded. Typed work remains pinned to
CPU 28 and alias work to CPU 29 on the same unreserved host. No measured samples
are discarded.

#table(
  columns: (1.8fr, 1fr, 1fr, 1fr),
  [Case], [Original baseline ms], [Initial unified ms], [Corrected ms],
  [free12-trace4], [73.862], [382.413], [66.666],
  [free12-tracen], [61.598], [119.149], [60.166],
  [free14-trace4], [331.406], [6827.879], [300.576],
  [fermion-4-tracen], [108.182], [619.525], [96.154],
  [axial12-trace4], [725.554], [117.697], [106.117],
  [mixed-colour-trace4], [4.747], [6.121], [5.842],
  [colour-fierz], [3.606], [4.544], [4.298],
  [historical-original, shared aliases], [90.514], [90.478], [83.454],
  [gluon-4-trace4, shared aliases], [225.804], [227.596], [231.307],
  [gluon-4-tracen, shared aliases], [299.413], [303.388], [298.554],
)

The corrected free14 operation is 22.7 times faster than the freshly measured
initial unified implementation; fermion4 tracen is 6.4 times faster. Both are
back near the original baseline in this cohort. These large changes remove the
identified regression. Small differences on the shared host do not establish
statistical significance, and not every case improves: mixed colour/Lorentz
and colour Fierz still take about 1.1 ms and 0.7 ms more than their baseline.
The gluonic alias routes are effectively unchanged at this sampling resolution.
No result was accepted or excluded using a speedup threshold.

Separate diagnostic CPU medians locate the reduction in the algebra/planning
stage. The “gamma reduction” clock includes the public algebra operation and
its following structural contraction, as it did before; it does not isolate
only the gamma kernel. These are separate calls and are not additive pieces of
the primary wall-clock median.

#table(
  columns: (1.6fr, 1fr, 1fr, 1fr, 1fr),
  [Case], [Before reduction ms], [Corrected reduction ms], [Before materialization ms], [Corrected materialization ms],
  [free12-trace4], [361.072], [42.757], [18.960], [18.541],
  [free14-trace4], [5886.094], [182.757], [99.593], [95.912],
  [fermion-4-tracen], [556.892], [38.558], [20.079], [19.939],
  [axial12-trace4], [105.705], [95.737], [6.051], [6.478],
)

Free14 still emits 26931 terms and the fermion4 scalar still emits 5204 terms.
The unchanged outputs and materialization costs exclude omitting requested
expansion as an explanation for the speedup. Child peak RSS spans
122984–154508 KiB before the correction and 119828–150012 KiB afterward over
this focused cohort; these peaks include imports, warmup and diagnostics, not
isolated kernel allocations.

Completed-owner reruns are separately pinned to CPU 32. Every sample checks all
returned statuses are `Complete` and the resolved expression is exact, outside
the clock. The following table repeats both algebra and ambient structural
contraction on the retained owner, excluding admission and materialization.
All three samples are retained, including the first.

#table(
  columns: (1.8fr, 1fr, 1fr),
  [Case], [Initial combined rerun ms], [Corrected combined rerun ms],
  [free12-trace4], [327.494057], [0.012529],
  [free12-tracen], [69.024713], [0.012100],
  [free14-trace4], [6324.468926], [0.010800],
  [axial12-trace4], [63.132548], [0.009070],
  [fermion-4-tracen], [546.921363], [0.008780],
)

The combined reruns now take roughly 9–13 microseconds. Separate same-policy
reruns remain inexpensive. These results apply to exact policies on an unchanged
owner; a changed domain, new policy or unfinished frontier still requires work.
All 90 retained rerun samples across both builds preserve the exact expression.

The corrected full Idenso suite passes 844 tests with 22 ignored tests and no
filtered tests, including the exhaustive axial test, in 124.41 s after a separate
3m29s build. All six new architectural regressions pass, covering unnecessary
structural scheduling, registry work, alternating completion policies,
root/definition invalidation and exact unfinished requests. The corrected
installed Python host passes all 87 API tests. Rust formatting, the docs builder
and all-target Idenso Clippy with `reference-cases` and warning denial pass.
The preexisting `proc-macro-error2` future-incompatibility notice remains a
dependency notice. Both this report and the updated Idenso architecture document
compile. Reproduce these gates with the tests and replay commands above.

For this historical correction, exact serialized-output equality transferred
the previously established
FORM/HEP certificates; no new FORM timing or broader native-phase capture is
claimed for this focused correction. The historical refusals, exact deferred
case and comparison limits above remain unchanged. Residual costs include
interface discovery, ordinary planning and explicitly requested materialization;
the initial global-registry reconstruction and latest-policy-only rerun costs
are no longer present on these paths.

== Massive routed fermion ladder: 2026-10-01

This follow-up measures the complete four-loop fermion ladder with three gluon
rungs: eight massive propagator numerators, expanded Standard Model couplings,
LMB momentum routing, symbolic Lorentz dimension, and four open external ports.
The bottom-quark mass, scalar prefactor and external Lorentz/colour indices are
retained. It is a different workload from the older massless, projected,
unrouted 525-term scalar benchmark. The request is
`simplify_algebra(gamma=True, color=True, epsilon=False, contract="fully")`,
with no work limit and no arithmetic expansion.

The implementation keeps the shared carrier, planner, contractor and Dirac
recurrence. Selected spinor factors are assembled before reducing a growing
prefix. The recurrence requests occurrence-local pair contractions from the
existing contraction owner and reuses each requested pair. Weighted routed
vectors remain factored; equivalent vector bodies with different external dummy
labels share their identity without identifying those dummy scopes. Selected
collection can consume incident coefficient ports while checking the combined
ordered boundary. Unrelated scalar and tensor spectators remain outside that
work.

Intrinsic replacements forward regional observations and interface proofs.
Successful contraction of a connection does not itself justify rescanning a
large emitted polynomial or revalidating unchanged spectators. Typed zero keeps
its admitted interface, while callbacks and unsupported rewrites retain the
checks needed for their actual result. These changes do not introduce aliases,
a second trace recurrence, or another contraction implementation.

#table(
  columns: (2fr, 1fr, 1.3fr),
  [Measurement], [Reduction seconds], [Completed rerun seconds],
  [Historical Python extension], [23.710816], [Not captured in this sample],
  [Candidate Python, sample 1], [0.848985], [0.025762],
  [Candidate Python, sample 2], [0.740976], [0.003168],
  [Candidate Python, sample 3], [0.681896], [0.009260],
  [Installed Python, sample 1], [0.797951], [0.015423],
  [Installed Python, sample 2], [0.768791], [0.003064],
  [Installed Python, sample 3], [0.886420], [0.017785],
  [Candidate native Rust runner], [0.777609], [0.001558],
)

The Python comparison uses the same admitted full input and requested settings;
all three measured calls and completed reruns are retained. The median falls
from the 23.71 s baseline sample to 0.798 s in the installed module, about 30 times
faster. The earlier candidate's median is 0.741 s. Each call
reports `Complete`, retains the exact logical port order, and reruns stably.
The native row is a separate frontend measurement. The historical extension is
identified by source revision `fe9263cf65e5c67259e010ccf2d01f9378c29453` and
core hash `ccff4edd64869690c92eec7989d15e2bbc8418a671cf12856368383a23a0cb67`.
The native receipt uses the release profile with LTO disabled on the same
AMD EPYC 9754 host. Startup, compilation, text parsing, admission, formatting,
and component materialization are outside the reduction clock.

The older saved projected fermion control remains at 106.1 ms median over three
calls; the saved projected gluon control takes 246.0 ms. Both match their exact
525-term and 6001-term references by factor-preserving cancellation and report
stable completion. Their settings are retained in the receipt; they are separate
workloads from the full massive ladder.

Native admission took 2.922 ms. A separate instrumented reduction took
718.598 ms, including an inclusive gamma clock of 598.080 ms, contraction of
131.272 ms, and candidate observations of 82.350 ms. Nested phase clocks overlap
and must not be added. This diagnostic replay follows the timed call, so caches
and observations can be warm. Regional candidate observations are included in
its scan counter; that counter is not a count of initial domain scans.

The small, source-controlled input is
`crates/idenso/benches/data/massive-fermion-ladder.expr`. The existing runner
replays it without executing a notebook. Record the checked-out revision,
dependency lockfile, build profile and environment with each new comparison:

```sh
cargo build --release --locked -p idenso --features reference-cases \
  --example contraction_phase_benchmark --config profile.release.lto=false

target/release/examples/contraction_phase_benchmark --reduce-input \
  crates/idenso/benches/data/massive-fermion-ladder.expr \
  /tmp/massive-fermion-ladder-replay algebra
```

The release Idenso run passes 844 tests; 22 tests are ignored and the separately
established baseline failure `tensor::tests::parsing::infinite_execution` is
excluded. New regressions check massive routed sums, distinct compound dummy
labels, repeated weighted vectors, coefficient consumption, untouched
spectators, ordered ports and stable completion against independent exact
Clifford components. The full fixture additionally agrees in all 1024 components
at three momentum/mass points with explicit Weyl and Gell-Mann matrices.
Its exact difference from the baseline factors to zero without explicit global
expansion. The final output is byte-identical to that independently checked
intermediate output, with SHA-256
`07a2a503443675b48a5f65d283069f3a3c716c1456c33237c86e39a99fe780ab`;
the receipt explicitly records reuse of those correctness certificates.

The first rebuilt Python result was also checked directly at all three numerical
points, including its actual carrier interface. All 1024 components agree,
with relative maximum-norm error below 1.1e-15. The final installed result has the
same exact plain and canonical output hashes; its carrier ports are checked
again. Its exact difference from the saved baseline factors to zero in 3.65 s
outside the timing window. The different serialization order is recorded
separately from the native output hash.

All 82 Python API tests pass in the actual notebook virtual environment,
including the restored Expression inheritance, tensor-preserving operations,
display, generated signatures and documentation examples. Idenso Clippy with
warning denial and formatting of the modified sources pass. The generated
stubs are synchronized into the bundled and live notebook packages. Both
notebook servers restart with HTTP 200; no notebook source is executed or edited
for this benchmark. The receipt records the loaded extension hash and verifies
that its inode and contents stay unchanged during measurement. The optimized
Idenso sources are unchanged across the display-related extension refreshes.

This report compiles with Typst. The repository-wide documentation gate retains
its preexisting GammaLoop API attestation failure for
`crates/gammaloop-api/src/python.rs`; this check is not reported as passing.

The output tradeoff is substantial: the candidate retains about 24.2 MB of Atom
data in Python (23.4 MB in the native runner), compared with about 0.98 MB in the
historical result, because the
contraction polynomial is left factored differently. Atom byte counts also
vary with process-local symbol IDs. The new representation is exact and complete,
but rendering, evaluator construction and later materialization may have
additional costs; those costs are not included in the reduction measurement.
No fresh FORM comparison or FORM timing is claimed for this fixture.

== Collected coefficients for the massive routed ladder

`simplify_algebra(collect_coefficients=True)` collects generated exact
coefficients by default; `False` retains the nested output. The option currently
controls contracted Dirac traces and bound vector factors. Other identity
families retain their existing coefficient normalization, and free traces retain
their factorization. This output choice belongs to the shared reduction planner
and does not authorize another identity family or a separate contraction engine.

The existing `TraceAlgebra` recurrence accumulates the requested metric/vector
pairings in Symbolica's exact polynomial representation. Its variable ordering
and emitted sectors are reused while equivalent surviving tensor structures
share a coefficient. Selected exterior coefficients participate in this
collection, allowing terms from different generated sectors to cancel. The
collector combines only certified terminal sectors with compatible ordered
interfaces and no remaining identity candidates or uncontracted connections in
their coefficients. Other sectors keep their existing processing path;
unsupported coefficient domains retain the exact Atom arithmetic.

Domain observations identify original scalar sums within vector weights, which
remain opaque polynomial variables. Power bases keep that protection even when
a contraction combines two source powers into a new power. The shared
`PolynomialExpressionExt` emitter restores these original expressions and applies
Symbolica's `expand_num()` to the generated coefficients: numerical factors are
distributed, so `2*(x+y)-2*x-2*y` can cancel, while products and powers of unrelated
scalar sums remain factorized. Unselected scalar spectators remain outside the
collection. This does not request global numerator expansion.

`DomainObservations::compose_terminal_product` now retains the observed nested
scalar scopes of unchanged spectators when attaching a generated terminal
region. Spectators from the same domain share their existing observation cache;
independently admitted spectators contribute their cached exact-syntax facts.
This requires no new candidate walk and lets later selected contractions
recognize scalar powers without rediscovering their interiors. It leaves the
contraction order and distribution policy unchanged. Collection uses the same
contractor and trace recurrence, and the existing tensor carrier preserves
ordered interfaces and typed zeros at emission.

The preceding massive-ladder receipts predate this option and describe the
nested output corresponding to today's `collect_coefficients=False`. Those
historical receipts remain unchanged. New measurements must state the option
explicitly and separate native and Python frontends, reduction and completed
reruns, and computation from construction or rendering. The existing native
runner selects either output form on the same fixture:

```sh
target/release/examples/contraction_phase_benchmark --reduce-input \
  crates/idenso/benches/data/massive-fermion-ladder.expr \
  /tmp/massive-fermion-ladder-collected algebra collected

target/release/examples/contraction_phase_benchmark --reduce-input \
  crates/idenso/benches/data/massive-fermion-ladder.expr \
  /tmp/massive-fermion-ladder-nested algebra nested
```

The final native comparison uses release builds with LTO disabled and the
`reference-cases` feature on the same AMD EPYC 9754 host. Collection work is
recorded in change `qzkuotrv`, based on `5434864a`; the preceding reduction work is
recorded in `ed14addf`. The companion JSON's `coefficient_collection` section
records source and binary hashes, dependencies, every sample, and raw receipt
hashes. Both modes enable gamma and colour identities with `contract="fully"`,
disable epsilon identities, and retain symbolic Lorentz dimension `D` without a
work limit. Gamma traces use a four-dimensional bispinor representation and
repeated-pair ordering; colour uses the standard fundamental SU(3) normalization.
Only `collect_coefficients` differs between the two modes.

#table(
  columns: (1.3fr, 1fr, 1fr, 1fr),
  [Native output], [Reduction median (s)], [Rerun median (µs)], [Atom bytes],
  [Collected], [1.062576], [49.311], [1,049,059],
  [Nested], [0.809103], [1,584.035], [23,377,307],
)

Each mode has three independent process samples. Collected reduction ranges
from 1.059142 to 1.101824 s; nested reduction ranges from 0.808023 to 0.868752 s.
Every call completes in one step and reruns stably. Collection adds about 0.253 s
to the median reduction while reducing Atom storage by a factor of 22.3.
Its plain output occupies 3,137,180 bytes, compared with 63,632,701 bytes for
nested output. Atom size depends on process-local symbol identities. These
clocks exclude startup, compilation, parsing, admission, serialization,
rendering and component evaluation; no isolated native peak-RSS measurement
was captured.

Admission medians are 1.431 ms collected and 1.417 ms nested. A separate
instrumented replay follows each timed call; the instrumented reduction medians
are 1.037457 s and 0.760196 s. Per-phase medians are:

#table(
  columns: (2fr, 1fr, 1fr),
  [Phase], [Collected (ms)], [Nested (ms)],
  [Interface discovery, outermost], [3.000], [2.892],
  [Candidate observations], [83.209], [84.270],
  [Shallow graph construction], [52.276], [51.890],
  [Gamma kernel, inclusive], [1,010.254], [648.924],
  [Colour kernel, inclusive], [13.920], [14.096],
  [Contraction, inclusive], [133.914], [132.299],
  [Output normalization, inclusive], [126.407], [96.066],
  [Coefficient merge, inclusive], [41.929], [0.168],
)

Nested phase clocks overlap and cannot be added; replay observations and caches
may be warm. The candidate counter includes affected-region observations, not
only initial domain scans. Both modes record 645 structural-contraction
frontiers per replay, with maxima of three states, three local terms, and four
generated alternatives. Those frontiers do not measure the entire Dirac
polynomial or global term count. Native reduction does no separate validation
or materialization in these phase receipts.

Both serialized results cancel exactly against the historical baseline using
`collect_factors()` followed by `factor()`, without global arithmetic expansion.
An independent contraction of the original factorized leaf components with
explicit Weyl gamma and Gell-Mann matrices agrees at all 1024 components at each
of the three saved mass/momentum points. These numerical checks specialize
`D=4`, use metric signature `(+, -, -, -)` and coupling one; the maximum relative
error is 1.071e-15 for either mode. Raw serialization does not certify the
original carrier's logical port order: component axes are matched by exact
representation and index identities, and carrier invariants are covered by the
native regressions separately.

The downstream comparison below uses existing Python core `81ef213c607c`,
not the new reduction extension. It parses each native output, specializes its
parameters and admits the result before constructing a direct one-core
Symbolica evaluator, with no JIT or common-pair optimization. Medians cover the
same three distinct physical points, each in a fresh process; startup, the
independent oracle, and reduction are excluded.

#table(
  columns: (2fr, 1fr, 1fr),
  [Downstream phase], [Collected (ms)], [Nested (ms)],
  [Parse, specialize and admit], [173.544], [9,395.301],
  [Construct evaluator], [0.715], [52.255],
  [Construct component points], [13.006], [13.392],
  [Evaluate 1024 components], [7.938], [33.200],
  [Total median], [204.894], [9,503.508],
)

The 46.4-fold downstream difference is dominated by preparing the serialized
expression. It is not an evaluator-only speedup or an end-to-end Python
reduction measurement. Independently computed phase medians need not sum to
the median total.

The native release suite passes 853 tests, with 22 ignored tests and the known
`tensor::tests::parsing::infinite_execution` failure excluded. The nine focused
collection tests and four shared polynomial-emitter tests pass. Clippy passes
with warning denial for `idenso`, `spynso3` and `symbolica-utils`; existing
dependency warnings are recorded separately in the receipt. These native and
independent native-output results remain separate from the installed Python
measurements below. No fresh FORM comparison was made for this fixture.

=== Installed Python result

The Maturin release build with the `native` feature produces core
`44b5f25f94de`, installed into both the source-host package and the notebook
virtual environment. Its full hash, build command, dependency hashes and
installation receipt are recorded in `coefficient_collection.python` in the
companion JSON. Reduction sources stayed unchanged throughout the build.
Concurrent graph and generator edits changed five other source files while
compilation ran; the build receipt records those changes and does not certify
that the extension includes their latest versions. The measured reduction uses
the saved numerator fixture, without regenerating its graph.

#table(
  columns: (1.3fr, 1fr, 1fr, 1fr),
  [Installed Python output], [Reduction median (s)], [Rerun median (ms)], [Atom bytes],
  [Collected], [1.035185], [0.241491], [1,098,752],
  [Nested], [0.779839], [19.237235], [24,213,378],
)

Each mode measures three calls on the same admitted tensor in one process,
using the explicit settings from the native comparison. All six calls report
`Complete`, preserve the exact ordered ports and return stable completed reruns.
Collected times range from 1.034375 to 1.063405 s; nested times range from
0.762212 to 0.827429 s. Collection adds 0.255 s to the median while reducing
the Python Atom size by a factor of 22.0. Admission takes 1.487 ms and 1.519 ms
respectively. Startup, model setup, admission, output conversion, serialization,
rendering and exact comparison are outside the reduction clock.

Each installed result independently cancels against the saved historical
baseline with `collect_factors().factor()`. Those checks take 0.689 s collected
and 3.735 s nested, outside the timing window. The loaded extension's inode and
contents remain unchanged throughout both measurements. The JSON retains every
wall and CPU sample, output hashes, logical slots and cumulative process peak
RSS; that RSS includes repeated reductions, retained results and the exact
comparison, and is not isolated reduction memory.

All 83 API tests pass against both the candidate package and the actual installed
notebook package, including inheritance, algebra, display, generated signatures
and documentation examples. Canonical, bundled, source-host and installed
tensor stubs share the recorded hash. Both notebook servers restart successfully
and return HTTP 200 on ports 2726 and 2727. The ladder notebook's source hash is
unchanged across the reload. These are runtime and server readiness checks;
they do not constitute a fresh execution or visual audit of every notebook.
