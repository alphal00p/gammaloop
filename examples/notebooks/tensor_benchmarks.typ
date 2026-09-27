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
