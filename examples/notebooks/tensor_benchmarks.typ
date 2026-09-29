#set text(size: 10pt)
= Tensor consolidation measurements

`fermion_ladder.py` remains the entry point for the notebook and command-line
benchmarks. `--suite consolidation` adds the R3/M0 case inventory without
changing the existing ladder algebra. The default is three rounds. Each round
rotates the release interpreters; every observation, subprocess log, expression,
and exact-check result is retained beside the requested JSON report. An existing
output directory is refused, so a later milestone cannot overwrite M0.

== Running a milestone

Use qualified release interpreters and bind their core hashes explicitly. The
M0 release corresponding to the plan's starting Rust sources is the private
sparse-contracted16 release with core `cfe4e655591faf6b871e12d32de18ad00e7ae76a27a793c7b66441149e69de8c`.
The earlier canonical `4eace950` release is a historical comparator, not this M0
baseline. Machine-specific interpreter and FORM paths are command arguments;
no fixture or runner depends on a scratch directory.

```sh
python examples/notebooks/fermion_ladder.py --suite consolidation \
  --interpreter baseline=/path/to/qualified/python \
  --core baseline=CORE_SHA256 \
  --interpreter candidate=/path/to/next/release/python \
  --core candidate=CORE_SHA256 \
  --form /path/to/form --cpu 8 --rounds 3 --output m0.json
```

For the first baseline omit the candidate. The same command and case definitions
apply to later milestones. `--check-only` runs the setup and exact gates without
the three-round cohort. `--cases` selects comma-separated names; the default
inventory is:

- `fermion-3-trace4`, `fermion-3-tracen`, `fermion-4-trace4`,
  `fermion-4-tracen`, `gluon-4-trace4`, `gluon-4-tracen`;
- `historical-original`, `historical-early`: the notebook's eight-vertex,
  eleven-internal-edge, externally polarized gluon fixture;
- `free4-trace4`, `free4-tracen`, `free12-trace4`, `free12-tracen`,
  `free14-trace4`, `free14-tracen`;
- `axial12-trace4`: the existing factored scalar-spectator cleanup control;
- `production-aa-aa`: the existing GL16 integrated-UV-start numerator in
  `crates/idenso/tests/fixtures/aa_aa_2l_gl16_integrated_uv_start_after_simplify_metrics.sym`.

The production capture is used verbatim. Its nested `hedge` and `edge` index
payloads use the existing `CookSettings.indices()` input adapter, outside the
clock. Its measured pipeline is gamma, colour, then metric simplification. It
is a captured post-metric production numerator, not a newly generated process
or a full UV computation. Its bytes and path are recorded. Its check is exact
cross-release output equality and the separately repeated staged pipeline;
it has no independently authored FORM oracle in this suite.

== Boundaries and interpretation

Primary clocks include one complete case reduction and output-list retention,
including final expansion/materialization where the existing case requires it.
Input construction, rule setup, checks, counting, output formatting and disposal
are outside wall, process and thread clocks. Each measured process starts fresh;
one fixed unmeasured warmup is checked and disposed before the clocks. The warmup count is recorded; no measured observation is discarded. `--calls` is fixed before the cohort and is
identical for every interpreter. Phase clocks are separate executions after
the primary interval, with output term counts and Atom byte sizes. They do not
sum to, or replace, the primary observations.

FORM reuses the ladder's independent-expression batch and internal CPU timer. The consolidation scheduler sets its own affinity before launching FORM, so reference processes inherit the requested CPU; each FORM receipt records that affinity and term counts.
One preserved pilot determines a batch targeting 200 ms. Free-four traces calibrate from 8,192 copies, doubling until the body CPU clock reaches 200 ms (at most 1,048,576 copies); all calibration trials are saved and the achieved count is fixed for the three reference rounds. No measured result is selected or retried. `--form-only` can supply a separately labeled reference without resampling completed Python observations. Setup, final checks, output and disposal lie outside its body timer.
Process wall time is also retained and has a wider boundary. A zero body timer
stops the cohort and preserves the observation. FORM medians are separate
reference measurements; they are not paired release comparisons. None of the
reported ratios claims exclusive access to the shared host.

For the historical ladder, the pure Symbolica column reuses the notebook's
outside-in scalar recipe, including its held expansion and local contraction
rules. Conversion of its scalar notation back to the typed metric notation is
an untimed oracle adapter. This compares complete user routes, not an isolated
measurement of interface-proof overhead. Both routes must equal the fresh FORM
polynomial with 9,652 terms. The physical ladder cases also retain the existing
explicit-metric, D-to-4 and independent HEP component checks. Free traces compare
with fresh FORM and D-to-4 where applicable. Axial traces can use a different
Schouten/epsilon basis: literal disagreement must remain recorded, and any
component certificate is described as exact checks at its stated deterministic
points, not a symbolic basis proof. FORM carries the axial scalar spectator as
one symbolic coefficient; Idenso retains its actual `(x+y)^8` factor.

The report keeps side medians, every candidate/baseline paired ratio, median
paired ratios and pair-win counts separately. Three rounds are a baseline
measurement, not a confidence interval. Setup and phase records are diagnostics;
small control differences and host interference must remain qualified.

== Complete production workloads

`--suite production` measures the M5 consumers with the same entry point. It
reuses the existing grouped-diagram smoke cards for `GL000`, `GL053`, and
`GL148`, and the existing `uv::tests::scalars_profile_new` test. These are three
individual three-loop aa→aa graphs and the bubble/sunrise/banana scalar UV
profile. They do not measure the complete 155-graph amplitude. The captured
`production-aa-aa` algebra case above remains a separate measurement.

```sh
python examples/notebooks/fermion_ladder.py --suite production \
  --production-build baseline=/path/to/saved/build \
  --production-build candidate=/path/to/candidate/build \
  --fixture-root /path/to/source/with/canonical/fixtures \
  --cpu 8 --rounds 3 --timeout 900 --output m5-production.json
```

Each build directory contains its `source` tree, `source-manifest.json` with
file hashes, and `build-results.json` with successful `test-lib` and `cli`
artifact paths and hashes. The runner verifies these inputs before and after
each process. Both builds must use the same compiler, profile, dependency
revision, model environment, and common compatibility prerequisites; retain
their build provenance with the report.

There is one fixed warmup per build and case, followed by three alternating
paired rounds. Complete fresh-process wall and child CPU clocks include
startup, generation, and the requested computation. Card construction, input
hashing, output parsing, and comparisons are outside those clocks. Graph cards
use one generation/integration core, twenty samples, and seed 1337. Each run
has a fresh state and integration workspace, a 64 GiB address-space bound,
and a recorded timeout. Failures and timeouts retain their logs and stop the
cohort; no observation is selected or retried.

The graph checks require successful generation, twenty evaluations, finite
integral/error estimates, matching recorded maximum-weight coordinates, and agreement within
relative tolerance `1e-9` and absolute tolerance `1e-12`. These finite numerical
checks supplement the exact tensor and HEP checks; they are not symbolic
identity proofs. The UV case retains its existing seventeen-point profiles,
scale exponents four through eight, and `pass_fail(-0.9)` assertion. FORM
comparisons belong to the independent tensor-algebra cases above, not the
complete production process.

The graph inputs in `fixtures/production/` are complete native runtime DOT
files. Historical compact DOT files omit required vertex numerators. The
saved baseline CLI generated ungrouped diagrams from the restricted
`sm-default` model; exact interaction, particle and directed-topology matches
selected the three inputs. Their local numerator bytes are unchanged,
historical grouped factors are retained, and every transported momentum
signature agrees. `provenance.json` records the mappings, source hashes and
generation commands. Both timing sides consume these identical files. Each
input passed native import, generation and the twenty-sample smoke before the
paired cohort. Python model serialization was not used for production input:
its expression print ordering changed the model fingerprint across frontends.

To prepare the files again, run the recorded generation and export commands
with the saved baseline CLI, also saving its model through `save state -o`.
Then run `tensor_benchmark_production_fixtures.py` with `--generated` pointing
to the exported DOT directory, `--model` to that state's `model.json`,
`--historical` to the original grouped DOT directory, `--generation-card` to
the CLI card, and `--output` to a new directory. The preparer finds the exact
matches and verifies the transport. Its replay reproduced all three saved
DOT files byte for byte.

The saved M1 production build exposed a startup incompatibility: unqualified
`𝑖` is a numerical literal in the pinned Symbolica revision. Vakint's symbolic
placeholder is now explicitly namespaced. The original failed build and log
are retained, and this same prerequisite is applied to both sides of the
production comparison. Production qualification is distinct from completion
of the paired M5 timing cohort.

== Diagnostics carried from the plan

The original scratch scripts are ported into shared diagnostic functions, not
installed as another contraction engine. Run one diagnostic with the same entry
point, release/core arguments, and historical fixture:

```sh
python examples/notebooks/fermion_ladder.py --suite consolidation \
  --cases historical-original --diagnostic partial-parse \
  --interpreter baseline=/path/to/python --core baseline=CORE_SHA256 \
  --rounds 3 --cpu 8 --output m0-network.json
```

- `validation` ports `validation_overhead.py` and `idempotent_cost.py`:
  stages 5, 7 and 8, ordinary/fused Schoonschip, metric simplification,
  expansion, absent-rule replacement and reinference. Exact identity is a
  recorded result, not an assumption.
- `partial-parse` ports `partial_parse_cost.py`. It reports the first 200
  normalized terms through the existing public operations. `list_dangling` is an explicit dangling-index query, not a measurement of the internal partial graph parser; no M2 parser gate is claimed for that row. The R1 rows compare
  a sum of exactly the first 50/100/200 terms with separate evaluation of that
  identical term set, preserving outputs until each clock stops. Construction
  and an exact expanded-result comparison are outside timing. The whole/single
  cost ratio uses total work on that common set; it does not compare the full
  original sum with an unrelated sample.
- `keep-rest` ports the scalar, foreign-tensor, nested-sum and AUTO gamma probes.
  Exceptions are retained as diagnostic outcomes. Required gates fail after the complete record is saved: current successful routes require exact algebra and rerun; `--gate-stage M0` requires the single-owner AUTO fix on candidates, while `--gate-stage M3` additionally requires all keep-rest and scalar-factor checks. The first interpreter is always recorded as baseline, with known limitations retained. Dual colour-port spectators and vector powers are also retained refusal cases, required at M3. The distinct-owner AUTO sum
  remains a known baseline/R2 limitation and a required M3 pass; successful
  single-owner contraction must not be presented as resolving that case.
- `factor-atom` and `factor-polynomial` share the original fixed-fixture
  transition grammar and one DAG builder from `factor_dp*.py`. They report rule
  application, incidence inspection, state/edge counts, contraction, forward
  emission and materialization separately, and verify the final 9,652-term
  polynomial. The Atom variant is an Atom forward pass over the retained DAG;
  it is not relabeled as the old eager-weight DP timing. Edge-coefficient byte
  totals exclude container, state, graph and allocator overhead. Neither route
  is a production API, a generic tensor parser, or an accepted M1 engine.

A diagnostic cohort retains each repeated diagnostic record and publishes
matching-row medians. Error or refusal records are never removed from those
medians by rerunning. The M0 record is descriptive: R1/R2 acceptance, known
limitations, and later milestone requirements must be evaluated explicitly
against the corresponding case and release identity.

== First M0 record

`tensor_benchmark_m0.json` keeps the compact raw clocks, case checks, source/core
identities, all calibration outcomes, and R1/R2 diagnostic rows. It contains 96
paired Python observations (three pairs for each of sixteen cases), six pure
Symbolica historical-ladder observations, and 45 final FORM reference samples.
The baseline is `cfe4e655`; the independently built R1+R2 candidate is `096d809f`.
This is the baseline for the consolidation milestones, not a claim that these
two changes improve every operation.

The following table uses process CPU milliseconds. “Paired” is the median of
three candidate/baseline ratios, not the ratio of the two side medians. “FORM”
uses the candidate side median and the separate pinned FORM body median.

#table(
  columns: (2.4fr, 1fr, 1fr, 1fr, 1fr),
  [Case], [Baseline ms], [R1+R2 ms], [Paired], [FORM],
  [fermion-3-trace4], [1.0053], [0.7972], [0.9665], [2.126],
  [fermion-3-tracen], [2.3636], [2.3645], [0.9562], [1.367],
  [fermion-4-trace4], [6.2585], [6.4362], [1.0611], [1.083],
  [fermion-4-tracen], [46.8716], [50.8886], [1.0857], [0.988],
  [gluon-4-trace4], [773.4732], [1073.9747], [1.3885], [1.220],
  [gluon-4-tracen], [1367.0334], [1477.0379], [1.0690], [1.073],
  [historical-original], [1202.6406], [1428.6503], [0.9953], [1.264],
  [historical-early], [526.7272], [519.2917], [0.9969], [1.714],
  [free4-trace4], [0.1397], [0.1425], [1.1418], [56.081],
  [free4-tracen], [0.1886], [0.2090], [1.0813], [91.309],
  [free12-trace4], [6.5311], [7.6598], [1.1728], [1.550],
  [free12-tracen], [21.3327], [22.7911], [1.0528], [3.275],
  [free14-trace4], [69.1315], [71.6457], [1.0451], [1.902],
  [free14-tracen], [449.3513], [564.7182], [1.2523], [5.591],
  [axial12-trace4], [3.0549], [3.1479], [1.0403], [2.602],
  [production-aa-aa], [1558.5961], [1568.6531], [1.0065], [—],
)

The historical side medians and paired medians differ noticeably, so both are
retained. The larger gluon-four/4D and free-fourteen/D differences also remain
in the record. Three observations per side do not identify their cause or
establish a general speedup or regression. Very small trace controls expose
public Python/typed-call overhead against a native FORM body with a different
boundary.

The R1 identical-termset gate passes all nine rows: candidate median whole-sum
cost divided by individual-term total cost is 0.820–1.022, below two. The
corresponding baseline range is 0.764–1.334. The R2 single-owner AUTO metric/gamma
cases pass. Distinct-owner AUTO sums, vector powers, and dual-colour spectators
remain explicitly recorded limitations for M3.

The initial full cohort stopped when a free-four FORM batch recorded 2 ms and
then 0 ms for 200 copies. These observations and the interrupted reference are
retained. Completed Python rows were not repeated. A subsequent source audit
found that parent-side FORM processes had inherited unrestricted launcher
affinity; all earlier references are therefore superseded for ratio reporting,
without deleting them. Final FORM-only references use CPU 8, calibrate tiny
traces to at least 200 ms, and preserve three fixed reference rounds. Their fresh
polynomial exports match the original exact-check exports byte for byte, with
matching source hashes and final term counts. The primary Python body is
unchanged across these orchestration revisions.

Free and axial trace checks preserve every false literal-equality flag. Where
four-dimensional Schouten/epsilon bases differ, three deterministic, nonzero
exact HEP component points certify those sampled components only. The original
gamma source uses the HEP network evaluator; output metric and epsilon leaves
use cached evaluations from the same library, with exact scalar arithmetic.
This is not a general symbolic proof. Generic-D free traces retain their exact
FORM polynomial comparisons.

The optional `validation` and both-order `factor-polynomial` timing commands
are retained in the JSON but deferred to a later quiet window so implementation
builds can proceed. No timing claim is made for those prototypes in this M0
record.


== Algebra simplification follow-up, 2026-09-29

This follow-up measures the implementation at `73c246d1`. The accompanying
changes update the color notebook, its regressions, and the timing harness;
they do not change the Rust algebra. The primary algebra comparison uses
`78a0a22f` and `73c246d1`, both with Symbolica `6a96c9d7`. It runs all sixteen
cases in three alternating fresh-process rounds, with one fixed warmup per
measured boundary (prepared and fresh). It retains all observations and checks against fresh FORM output where an independent FORM case exists.
The complete ladder route now uses the public contractor's default graph order.
The former fixed fixture order remains available through
`--contraction-order explicit`. To reproduce the current complete route with
qualified releases:

```sh
python examples/notebooks/fermion_ladder.py --suite consolidation \
  --interpreter checkpoint=/path/to/checkpoint/python --core checkpoint=CORE_SHA256 \
  --interpreter current=/path/to/current/python --core current=CORE_SHA256 \
  --ladder-route checkpoint=factorized --ladder-route current=factorized \
  --contraction-order graph --complete-algebra --rounds 3 --calls 1 \
  --form /path/to/form --cpu 8 --output /tmp/tensor-full-stack.json
```

Two clocks answer different questions. Prepared reduction excludes input and
rule construction; complete algebra includes them and the requested output
materialization. Both exclude process startup and correctness checks. FORM's
batched body clock excludes setup and startup, so ratios against FORM use the
prepared native clock. Graph generation, evaluator construction and integration
are outside this algebra cohort. The host is shared; CPU affinity does not establish exclusive
access to memory or other machine resources.

The two historical ladder names select the original and early FORM schedules.
With graph ordering selected, both native cases use the same automatic
contraction strategy. Phase clocks are separate diagnostic executions, not
parts to add to the independently measured complete time. Intermediate output
inspection occurs after the phase clocks, so it cannot resolve aliases before
a later timed phase.

=== Complete algebra results

`tensor_full_stack_timing.json` retains all 102 measured observations, 34 check
processes, 15 FORM references, phase samples, and source/core identities. The
cohort completed with unchanged source/core hashes. Values below are median
process CPU milliseconds. Prepared and fresh clocks are separate executions;
small timing differences between those columns are not incremental costs.

#table(
  columns: (2.4fr, 1fr, 1fr, 1fr, 1fr),
  [Case], [78a0a22f fresh], [73c246d1 fresh], [73c246d1 prepared], [FORM body],
  [3-loop fermion, 4D], [5.698], [5.807], [3.840], [0.325],
  [3-loop fermion, D], [28.633], [27.783], [25.467], [1.580],
  [4-loop fermion, 4D], [40.408], [40.768], [36.564], [5.655],
  [4-loop fermion, D], [299.946], [295.779], [293.545], [47.500],
  [4-loop gluon, 4D], [722.453], [235.601], [236.479], [637.000],
  [4-loop gluon, D], [842.312], [313.120], [337.623], [852.000],
  [Historical, original FORM order], [450.055], [92.416], [90.301], [702.000],
  [Historical, early FORM order], [451.075], [88.526], [94.477], [160.000],
  [Free trace 4, 4D], [1.110], [1.114], [0.936], [0.001740],
  [Free trace 4, D], [1.127], [1.129], [0.935], [0.001354],
  [Free trace 12, 4D], [408.183], [406.151], [403.604], [2.860],
  [Free trace 12, D], [275.974], [276.308], [267.935], [4.134],
  [Free trace 14, 4D], [1800.488], [1756.257], [1799.228], [20.800],
  [Free trace 14, D], [1396.468], [1412.707], [1390.074], [57.750],
  [Axial trace 12, 4D], [1346.019], [1325.472], [1359.298], [0.665],
  [Captured production numerator], [4195.621], [4328.220], [4640.229], [—],
)

The complete default-order gluonic routes improve by about 3.1× (four-loop
4D), 2.7× (four-loop D), and 4.9× (historical original) relative to the paired
checkpoint. Fermionic and free traces remain approximately at checkpoint
performance. The native prepared gluon times are below FORM's body times,
while the fermion ladders remain roughly 6–16× slower. Free and axial traces
have much larger gaps. The original and early historical FORM schedules have
the same exact output but substantially different runtimes.

The older M1 record is materially faster on trace workloads: for example,
its axial control was about 1.9 ms, versus the current prepared 1.36 s. M1
used Symbolica `06906976` and the earlier expanded-result frontend, so these
historical numbers do not isolate the effect of the present amendment.
They do establish that flat timings against `78a0a22f` are not evidence
that the complete consolidation retained the earlier trace performance.

=== Trace and production-capture breakdown

The compact record’s `phase_diagnostics` contains the corrected separate phase
measurements; the original raw report’s earlier diagnostic metadata remains
preserved. These phase measurements locate most current trace time before final
conversion. The four-loop D fermion case spends about 245 ms in gamma
reduction, 33 ms materializing its trace aliases, 8.5 ms converting to a
polynomial, 25 ms routing the polynomial, and 6 ms emitting the scalar result.
The free fourteen-gamma 4D case spends about 1,605 ms in gamma reduction and
130 ms in materialization; generic D spends about 804 ms and 496 ms,
respectively. These diagnostic medians are not a reconstructed total.

The axial control spends about 1,123 ms in gamma reduction. Its explicit
spectator-preserving alias reconstruction takes about 157 ms, followed by
21 ms of expansion. The literal `(x+y)^8` spectator is retained; its
unrelated scalar expansion is not part of the requested result.

The captured production numerator spends about 2,303.5 ms in gamma, 1,482 ms
in color, 0.007 ms in its final no-op index contraction, and 381 ms resolving
the aliases. This case starts from a captured post-metric expression and
measures algebra only. A fast final contractor does not imply a fast complete
simplification pipeline.

=== Where the trace overhead occurs

`tensor_benchmark_trace_diagnostics.json` records an additional read-only
comparison on CPU 20: one warmup and three samples for each existing scheduling
setting. Gamma-only processing takes about 222 ms for free-12/4D, 139 ms for
free-12/D, and 56 ms for axial-12; the default pipeline takes about 389 ms,
243 ms, and 1,193 ms in that same phase. Every variant has the same exact
materialized output as the independently checked primary result. These
comparisons change scheduling through the existing API; no alternative engine
or production optimization was added.

The free-12 results contain 496 typed aliases in 4D and 211 in D. Sampling
locates repeated structure work in `TraceDefinitions::retain`, which infers
the sums before making alias handles, and `with_aliases`, which validates the
definitions. Index observation, interface inference, and alias registration
appear prominently. These are inclusive, overlapping stacks, not percentages
that can be added into a time budget. Production builds typed factored aliases;
the earlier direct expanded trace output is now restricted to tests.

For axial-12, gamma plus metrics takes about 98 ms, while gamma plus epsilon
takes about 925 ms. The alias count stays at 103, but default processing grows
the definition payload from 45,606 to 481,987 bytes, with the largest definition
growing from 2,217 to 225,866 bytes. The global registry-epsilon flag enters
collection for every domain, including metric-only definitions. That collection
path still parses without a depth limit; it is distinct from the amended
shallow contraction planner. The larger payload also increases later interface
reconstruction, as the separate axial phase clocks show.

The perf captures cover short complete diagnostic processes, including setup
and receipt hashing. They identify hot call stacks rather than supplying
another accounting of the primary clock. Together with the within-release
setting comparisons, they identify shared alias admission and domain-pass
scheduling as concrete remaining bottlenecks. Final polynomial conversion
alone does not explain the current trace slowdown.

=== Unfinished full gluon-scattering check

The complete generated `gg → gg` test with symbolic D and SU(N) still does not
finish its two-stage color reduction. A fresh run on CPU 22, retaining every
original operation and assertion, reached the 1,800-second timeout (1,757
seconds of child CPU; about 633 MiB peak RSS). Tensor construction and the
pre-color pipeline reached the first stage after about 0.15 seconds. That
stage, including materialization, finished after about 326 seconds; the
second color stage was still running at timeout. Later polarization,
closed-form, and cross-section checks were not reached, so this is an
incomplete workload, not a passing test or completed timing.

`tensor_full_stack_scattering.json` retains the command, source/core hashes,
stage log, and timeout receipt. The earlier `78a0a22f` run was also stopped in
its first color pass after more than 1,037 seconds; the original checkpoint, migrated
checkpoint, and current inputs at that boundary were byte-identical. The
small color-word notebook fix therefore does not establish that this much
larger color pipeline is fast or fully validated.

=== Color notebook regression

The explicit word and its adjoint are individually valid. Their internal `k`
labels represent independent sums; flattening the two unscoped expressions
creates four compatible occurrences of `k`. The admission check saturates its
count at three and rejects the product before color algebra. This behavior is
identical at both measured releases. The former notebook compacted the word
into a chain before multiplication, which removed the exposed internal label.

The migrated HEP notebook demonstrates both compact chains and explicit words.
It gives the two explicit operands distinct
`wrap_indices(scope, dummies_only=True)` scopes, preserving their shared free
ports. Exact SU(2), SU(3), and SU(5) norm regressions pass; the independent HEP
SU(3) component network gives `16/3`. The weighted complex example retains its
`400/3` component oracle. The complete notebook executes successfully, and raw
unscoped overuse remains rejected. Existing exact trace-equality assertions
compare explicitly materialized alias outputs.

The complete gamma notebook also passes its existing Clifford, trace and HEP
component checks. Its slash example now registers vector momenta and contracts
the typed gamma Lorentz port explicitly. Both HEP notebook regression tests
pass together; benchmark buttons remain opt-in during notebook tests.
