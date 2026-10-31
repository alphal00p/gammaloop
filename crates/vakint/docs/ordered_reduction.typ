#set document(title: "Ordered RustRed reduction")
#set page(paper: "a4", margin: 18mm, numbering: "1")
#set text(font: "Libertinus Serif", size: 10.5pt, fill: rgb("243343"))
#set par(leading: 0.65em, spacing: 0.75em)
#show heading: set text(fill: rgb("173f65"))
#show raw: set text(font: "DejaVu Sans Mono", size: 8pt)
#let navy = rgb("173f65")
#let pale = rgb("eaf2fa")
#let panel(body) = block(width: 100%, inset: 10pt, radius: 4pt, fill: pale, body)

= Ordered RustRed reduction

The ordered strategy combines equal, not-yet-reduced integrals across the
request before expanding them. It is a scheduler over the existing checked
reduction operations and Vakint's default strategy. The RustRed library's own
default remains depth-first, which Vakint users select explicitly when they
need completed-result caches.

#panel[
  *Scope.* Existing reduction tables and topology matching; exact rational
  coefficients; a #raw("BTreeMap") frontier; no ordered reduction-result cache.
  No offline rule rewriting, synthesized composite rules, or loop-specific
  scheduling heuristics are required.
]

== Selecting the strategy

Use the same evaluation methods and master modes as for depth-first reduction.
The strategy belongs to the existing RustRed application settings and defaults
to #raw("vakint::DEFAULT_RUSTRED_REDUCTION_STRATEGY"), i.e. ordered. Select
depth-first reduction explicitly, for example to use completed-result caches:

```rust
use vakint::{RustRedReductionStrategy, VakintSettings};

let mut settings = VakintSettings::default();
settings.rustred_application_options.reduction_strategy =
    RustRedReductionStrategy::DepthFirst;
```

#raw("RustRedReductionStrategy") reexports RustRed's strategy enum. Selecting
it changes scheduling, not the definitions, normalization or requested numerical
precision of the masters. RustRed leads the default evaluation order; FORM-based
backends follow it only when FORM is usable. The Python constructor exposes the
same choice:

```python
from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod

engine = Vakint(
    rustred_reduction_strategy="depth_first",
    evaluation_order=[VakintEvaluationMethod.new_rustred_method()],
)
```

The accepted strings are #raw("ordered") (the default), #raw("depth_first"),
and #raw("depth_first_without_memoization"). Unknown values
produce an error. The no-memoization choice, exposed in Rust as
#raw("RustRedReductionStrategy::DepthFirstWithoutMemoization"), recomputes repeated
descendants without completed-result memoization, lookup/insertion bookkeeping,
or persistent caches, including dependent factor reductions. A zero cache budget
or subtraction of bookkeeping time is not this execution path. Immutable
program preparation remains shared. Python instances configured with ordered
or no-memoization reduction reject explicit cache operations, with an error
explaining the #raw("depth_first") opt-in, and skip prepared caches during data
loading. In Rust, #raw("OfficialExternalDataOptions::default()") and
#raw("load_external_data") likewise skip prepared caches;
#raw("OfficialExternalDataOptions::for_strategy(RustRedReductionStrategy::DepthFirst)")
and #raw("load_external_data_with_prepared_caches") load them. Rust evaluations
in these modes never consume attached DFS completed-result backings. Ordinary
numerical, symbolic and retained-master evaluation retain their existing public
interfaces.

== What is ordered

For routed applications, define the support of a canonical integral by
$s(I) = abs({i : a_i > 0})$: the number of distinct active denominator slots.
Extra powers of an existing propagator do not add slots. The drawn incidence
edges of the original Feynman diagram need not equal this count: repeated
quadratic denominator forms can coalesce.

The saved-application schedule processes support sizes downward. At a given
support size it finishes routing into canonical owners before processing those
owners. Same-owner IBP children strictly descend in the saved local ordering;
an exit into another reducible topology has smaller support. Declared masters
are terminal sinks, including masters on the same support.

The phase schedule is fixed for one request and advances once through these
stages. It does *not* require each individual parametric rule to be used once.
An owner stage can apply a template to many different integrals. Inside that
stage the greatest accumulated operational key is expanded first, ensuring
that contributions from larger keys can meet before that key is expanded.

Sealed lower-loop artifacts use their checked rule ordering and symmetry
canonicalization through the existing artifact adapter. Preferred factorization
remains a checked complete action: the ordered scheduler does not rewrite its
internal product reduction or its isolated dependency memoization. A direct
closed/candidate artifact whose saved order is not support-primary uses one
certified exact-order phase, rather than inventing downward topology stages.
The scheduler is generic; which families are
available depends on the loaded or vendored registry, not on new loop-count
heuristics.

== Where cancellation happens

Suppose the checked identities give $X -> 2 U + V$ and $Y -> -2 U + W$,
with both $X$ and $Y$ greater than $U$ in the reduction ordering. A request
containing $X+Y$ accumulates $2U-2U$ at one frontier key before expanding $U$.
The zero entry disappears, leaving only $V+W$ to reduce. Depth-first reduction
can reuse a completed reduction of $U$, but generally computes it before the
two contributions cancel.

Two integrals combine only when their *canonical operational identities* agree:
family, transport/owner state, signed powers and compatible coefficient channel.
Equal-looking raw power vectors from different families are insufficient.
Vakint keeps mass scaling, external factors and normalization in the existing
owners and batches compatible unit-mass channels together. It does not group
terms by printed master names or floating-point approximations.

Exact complex rational coefficients use the canonical carriers $1$ and $i$;
their rational magnitudes enter the native weights. Supported rational
combinations of native parameters are split in the same way when their
denominators are real and the native polynomial representation can hold them.
Other scalar factors remain opaque, including floating coefficients with their
original precision. This avoids fragmenting compatible frontiers by a numeric
phase without adding complex poles or changing the input.

The engine reuses checked rule selection, specialization, transports, nonzero
conditions and descent checks. Conditions from applied identities remain
required even when their contributions later cancel. Canceling a key before
applying its rule can avoid that rule and its conditions. Coverage warnings,
preserved unresolved integrals and numerical capability errors retain their
existing meaning.

== Weight arithmetic and worker cores

The scheduler's control flow depends on the weights only through exact-zero
erasure of merged keys, so how the weights are computed is a separate choice.
#raw("settings.ordered_weight_arithmetic") (#raw("vakint::OrderedWeightArithmetic"),
a selection of RustRed's #raw("OrderedWeightArithmetic"), reexported as
#raw("vakint::RustRedOrderedWeightArithmetic")) selects it for the closed one-
to three-loop artifacts, the four-loop candidate programs and the saved
five-loop application:

- #raw("exact") keeps exact rational functions of the dimension on the
  frontier.
- #raw("modular_sampled") records the weight-independent trace once, with
  small modular probes deciding the zero erasures, replays it over 64-bit prime
  fields at sample points and interpolates every terminal coefficient per prime
  (Thiele continued fractions). Every replay re-checks each recorded erasure and
  fails the request on a disagreement.
- #raw("modular_polynomial") runs scheduler passes whose weights are rational
  functions of the dimension over prime fields, without interpolation. Until
  RustRed implements it, a request fails with a
  #raw("modular_polynomial weight arithmetic is not implemented yet") error.

Both modular modes lift their prime images over the rationals (Chinese
remaindering and rational reconstruction) and verify every coefficient on fresh
primes. They need a coefficient context with at most one variable, which holds
for Vakint's native requests. Every mode expands the same keys and returns the
same terms and nonzero conditions; only the monitor's frontier monomial counters
differ (modular weights count one monomial per key). The default is
#raw("vakint::DEFAULT_ORDERED_WEIGHT_ARITHMETIC") (#raw("exact")), a single
constant that also sets the Python default.

The text form is RustRed's canonical syntax #raw("exact"),
#raw("modular_sampled[:key=value,...]") or
#raw("modular_polynomial[:key=value,...]"); #raw("modular") is an alias of
#raw("modular_sampled"). The keys are #raw("threads"), #raw("seed"),
#raw("verification_primes"), #raw("max_primes") and #raw("cross_check_exact")
for both modular modes, #raw("batch_points"), #raw("verification_points") and
#raw("max_points_per_prime") for the sampled mode, and #raw("lanes"),
#raw("denominators") (#raw("factorized") or #raw("dense")) and
#raw("max_degree") for the polynomial mode. The historical spellings
#raw("batch"), #raw("checks"), #raw("max_points") and #raw("cross_check") are
accepted. #raw("settings.rustred_application_options.ordered_weight_arithmetic")
is derived from this selection and must keep its default.

#raw("settings.n_cores") is the single source of every Vakint and RustRed worker
count. It defaults to #raw("vakint::default_n_cores()"): the cores available to
the process (affinity masks and cgroup quotas included), at most ten, and one
without a Symbolica license. Zero, more cores than the process may use and
more than one core without a license are rejected before any work. One core
runs every computation on the calling thread, without a thread pool or a
spawned thread; only an enabled progress monitor adds its presentation thread.
A modular arithmetic runs on #raw("n_cores") workers unless its selection sets
#raw("threads"). The same count verifies the five-loop routes while the saved
application loads (#raw("preparation_workers"), derived and otherwise left at
its default; data loaders take it from
#raw("OfficialExternalDataOptions::n_cores")). Results are identical for every
core count.

```rust
let mut settings = vakint::VakintSettings::default();
settings.n_cores = 8;
settings.ordered_weight_arithmetic = "modular_sampled".parse()?;
```

```python
vakint = Vakint(n_cores=8, ordered_weight_arithmetic="modular_sampled")
```

The rung fixtures and the
#raw("ordered_reduction_rungs") example read #raw("VAKINT_N_CORES") and
#raw("VAKINT_ORDERED_WEIGHT_ARITHMETIC") (without a #raw("threads") option the
modular workers follow the cores) and record both in their report.

== Reading the monitor

Enable the existing terminal monitor or progress snapshot as usual. The
ordered fields describe *one native request*. A public evaluation can contain
several incompatible channels or family requests, so the local phase counter
can restart for a new request; the request identifier makes that scope explicit.
The separate ordered-batch and compatible-cohort counters survive nested
activities. Read their progress lexicographically: the outer batch advances,
then its cohort counter advances, while the native phase counter is local to
the current request. Different batches need not have the same cohort total.

#table(
  columns: (0.43fr, 1fr), inset: 7pt, stroke: 0.4pt + rgb("d5dee7"),
  fill: (x, y) => if y == 0 { pale } else { none },
  table.header([*Field*], [*Meaning*]),
  [Ordered batches completed / total], [Loop-count/master-mode groups collected
    by this evaluation. The total is fixed before dispatch. Other backends and
    final master evaluation are outside this counter.],
  [Compatible cohorts completed / total], [The current collected batch's fixed
    family/program, mass and spectator carrier groups. It resets only when the
    outer batch advances. This is not a total over a heterogeneous whole job.],
  [Completed / total phases], [A monotone count within the fixed request
    schedule. Empty stages may complete immediately. It is not a fraction of
    elapsed or remaining time.],
  [Topology stage], [Current action kind, positive support size and canonical
    owner, when defined. Local rank, total powers and $D=A-R$ are not individually
    required to decrease.],
  [Input / current / peak frontier keys], [Distinct pending operational keys.
    This is an integral-count census, not the full byte size of coefficient
    polynomials or of the caller's Symbolica expression.],
  [Frontier / terminal coefficient monomials], [Separate current/peak sums of
    numerator and denominator monomial counts in the already-materialized
    rational coefficients. This size proxy reveals coefficient swell or
    shrinkage even when the integral count is unchanged.],
  [Processed keys], [Operational keys whose checked action was processed.
    This does not supply an advance bound on the number still undiscovered.],
  [Emitted RHS contributions], [Nonzero frontier and terminal contributions
    produced after the checked step's local coalescing and before global
    grouping. It is an event count, not a current expression size.],
  [Nonzero merges / exact cancellations], [Insertions that combined with an
    existing key and remained nonzero / removed the resulting zero key.
    Input and intra-rule coalescing are excluded.],
  [Terminal keys and merges], [The separate Native core output map, including
    preserved unresolved integrals and its own merge/cancellation counts.
    This precedes terminal substitutions and final catalogue canonicalization;
    it is not the final declared-master count.],
)

Watch frontier size together with processed keys, merging, cancellations and
coefficient work. A shrinking frontier indicates less pending integral structure;
a growing frontier can still be normal progress while a large rule fans out.
One cancellation can avoid a huge subtree, while many cheap merges may save
little time. These counters do not measure avoided CPU seconds. Existing RSS,
machine-available RAM, phase age and heartbeat displays remain separate.

The coefficient-monomial census updates incrementally without CAS or frontier
scans. It excludes the popped active coefficient, guards, temporary arithmetic,
keys, Atom encodings and bytes/RSS. Terminal counts precede saved substitutions.
Integer coefficient bit sizes are not counted either.

Owner indices in telemetry are zero-based. The fixed schedule includes empty
stages; completion advances only once a stage has been drained or skipped.
The topology-action hierarchy follows the useful outer structure of
#link("https://liacs.leidenuniv.nl/~plaata1/theses/ThesisBenRuijl.pdf")[Ben
Ruijl's thesis, section 3.7 and figure 21]; it does not assert that each raw
parametric template is a once-only action.

== Replaying native requests

Engine development does not need Vakint in the loop. Setting
#raw("VAKINT_DUMP_ORDERED_REQUEST=DIR") makes every native reduction request
(one per compatible mass/spectator channel, at one to five loops) write
#raw("DIR/<stem>.request.json") before it runs: the loop count, the program
identity (family fingerprint and external-data files, or the saved 5L
application path), the strategy, Ordered weight arithmetic, master mode and
channel label, and the exact input keys with coefficients as Symbolica
canonical strings of rational functions in the family dimension. Caller-side
nonzero conditions are recorded for reference; RustRed receives only the
inputs.

When no monitor is attached, #raw("<stem>.phases.jsonl") receives the native
ordered counters at every phase end, flushed per line, so truncated runs keep
them. #raw("<stem>.result.json") records BLAKE3 digests of the exact native
output: its terms in native order and the sorted set of its guards (layout
#raw("vakint.ordered-native-output.v1")). #raw("DIR/index.jsonl") lists all
files. With #raw("VAKINT_DUMP_ORDERED_REQUEST_ONLY=1") the requests are written
and their native reductions skipped; the evaluation then fails explicitly.
Without these variables nothing is written or allocated.

RustRed's #raw("ordered_replay") example (crate #raw("rustred-app")) loads the
same program from an external-data checkout, replays a request and prints one
JSONL record per phase with its wall time and counters. It compares the counters
and the output digests with the recorded ones, can stop after a phase and can
capture the complete weighted frontier at the start of a phase. It replays with
exact weights unless #raw("--arithmetic SPEC") selects another arithmetic: the
output digests of a request recorded under modular weights still compare, but
its monomial counters do not, and saved five-loop requests recorded with
non-default options are refused. The recorded preparation workers (Vakint's
#raw("n_cores")) do not count as such an option.

== Cache isolation and limitations

Ordered reduction neither consumes nor produces the depth-first engine's
completed-result cache. An attached DFS cache is not queried by ordered
evaluation. Python ordered loading skips prepared caches; explicitly requesting
cache operations in an ordered Python instance produces an error. Immutable
program preparation and the active external-data registry remain reusable:
they are not reduction-result caches.

Global cancellation deliberately avoids computing some individual unit-integral
reductions, so a completed result for each such integral does not exist to save.
Future ordered caching needs its own semantic type and binding. No partial
frontier or weighted amplitude is written as a completed DFS identity.

Keeping more unfinished contributions visible can increase frontier memory.
Exact coefficient combination can also cost more than the reductions it avoids.
The ordered strategy therefore has no universal speedup promise. New comparisons
measure three arms directly: ordinary DFS with intrinsic per-run memoization,
DFS without that memoization/bookkeeping, and Ordered. Prepared/saved caches and
previous-run reduction state are unloaded in all three arms. Cached historical
measurements remain separate cohorts.

== Reproducible physical comparisons

The collector #raw("examples/ordered_reduction_rungs.rs") uses the same native
QCD graph, color/Dirac contraction, routing and scalar-lowering fixture for both
engines. All new rung fixtures, numerical references and measurements use
#raw("LoopNormalizationFactor::MSbar"). One loop has no vertical gluon rungs;
each successive loop adds one.
Both zero-momentum metric and complete quadratic transverse projections are
tested. Prepared scalar inputs are authenticated and reused independently of
reduction, preserving full mass and numerator conventions.

From the repository workspace root, for a functional comparison:

```sh
cargo run --locked --release -p vakint \
  --features feynkit-ingress,experimental-rustred --example ordered_reduction_rungs -- \
  --work-dir /absolute/workspace/rung-L1-q0 \
  --loops 1 --projection zero-momentum --engine all --build-profile release
```

Through four loops, use separate #raw("--engine dfs"),
#raw("--engine dfs-no-memo") and #raw("--engine ordered") invocations
with fixed affinity, #raw("--no-monitor")
and a release build for qualified fresh-process timing. Declare the actual build with #raw("--build-profile");
this is a label, not compiler introspection. The
#raw("--engine all") mode creates three fresh engines in one process and is useful
for exact parity; it is not three fresh-process measurements. Existing report
files are never overwritten. Pass #raw("--artifacts PATH") to select an
authenticated checkout or installation; explicit four- and five-loop data is
needed. #raw("--prepare-only") saves input without reduction.
Five-loop physical runs use the ordered engine only, zero momentum followed by
quadratic. The accepted saved zero-momentum result supplies the exact historical
reference after conversion from its original pySecDec normalization through the
existing normalization API, before any Laurent truncation. Its raw payload and
original timing remain unchanged. No new five-loop depth-first run is required.

After separate runs,
#raw("--work-dir OUT --compare-reports DFS_DIRECTORY DFS_NO_MEMO_DIRECTORY ORDERED_DIRECTORY")
authenticates their saved payloads and compares exact retained expressions and
serialized numerical centers/error estimates without reducing again. These are
mathematical bounds only for components marked #raw("Bounded"). It writes a bound
#raw("parity.json") suitable for the benchmark adapter.
For current performance qualification, run each DFS variant once and Ordered
three times. Reuse those same saved baseline reports for the three cheap exact
comparisons. Their repeated references are not repeated timing samples: the
report labels each DFS baseline as one observation and gives only Ordered a
three-run median/range. Prior balanced measurements retain their original counts.

The comparison headline is #raw("reduction_math_seconds"): measured
non-overlapping math spans covering IBP/routing and exact weighted terminal
coefficient assembly/grouping. Necessary terminal substitutions, projections
and lower-factor reductions are included, as is request-local reducer/dispatch
state construction for both engines. Immutable artifact loading, deserialization
and program import, graph/data preparation, matching,
initial scalar lowering, MATAD normalization, numerical evaluation and
serialization are excluded. This is broader than a Native-kernel IBP timer.
It sums disjoint math blocks, including later exact coefficient combination;
normalization intervals themselves remain outside. Natural request-state cleanup
inside the current call is included, while engine-lifetime teardown after the
call is excluded. No forced cache clearing changes the execution path.

The collector retains gross math time and subtracts the *measured* lazy owner
preparation wall time, with one owner-preparation worker. The checked net,
interval count, fixed scope and completeness flag remain in the report. Invalid
telemetry makes the net unavailable without invalidating the scientific result.
The optional capture requires #raw("experimental-rustred"); a build without it
reports no math duration. The public #raw("retained_reduction_seconds") timer
remains separate and includes frontend lowering/finalization. Missing math
times never fall back to that public timer. In particular, the historical
40-hour public observation supplies no isolated reduction-math duration.

Prepared caches and previous-run result caches remain unloaded for these
comparisons. Each new report states *primary-engine completed-result*
memoization as a boolean, true only for ordinary DFS. Ordered's opaque checked
factorization may use a fresh private memo internally; it does not feed a
reusable ordered-result cache. The no-memo DFS arm disables these dependent
result memos too. Reusing the common authenticated scalar input is frontend
preparation, not a reduction-result cache. Data loading and
the subsequent numerical reference calculation remain separate.
Qualified comparisons disable monitoring in every arm. Dashboard/snapshot
checks, including the large five-loop runs, retain their observer overhead and
are diagnostic observations; independent RAM guarding remains active.
For one through four loops the collector chooses the smallest validated finite
source cap among surviving dependencies, or an explicit 300-digit target if
all are arbitrary precision. It records achieved component accuracy and errors
in MSbar at unit UV mass, unit renormalization scale and unit coupling. Five-loop masters
are uncomputed, so only exact retained comparison is available there.

Qualification proceeds through both projections at increasing loop count.
Only complete comparisons measured below fifteen minutes qualify for ordinary
CI; longer controls remain explicit. The benchmark source and measured records
live in #link("benchmark/benchmark.typ")[the benchmark report], with append
instructions in #link("benchmark/maintenance.typ")[its maintenance contract].
