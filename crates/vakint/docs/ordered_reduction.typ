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
request before expanding them. It is an alternative scheduler over the existing
checked reduction operations. Depth-first reduction remains the default.

#panel[
  *Scope.* Existing reduction tables and topology matching; exact rational
  coefficients; a #raw("BTreeMap") frontier; no ordered reduction-result cache.
  No offline rule rewriting, synthesized composite rules, or loop-specific
  scheduling heuristics are required.
]

== Selecting the strategy

Use the same evaluation methods and master modes as for depth-first reduction.
The strategy belongs to the existing RustRed application settings:

```rust
use vakint::{RustRedReductionStrategy, VakintSettings};

let mut settings = VakintSettings::default();
settings.rustred_application_options.reduction_strategy =
    RustRedReductionStrategy::Ordered;
```

#raw("RustRedReductionStrategy") reexports RustRed's strategy enum. Selecting
it changes scheduling, not the definitions, normalization or requested numerical
precision of the masters. RustRed must still be selected in the evaluation
order. The Python constructor exposes the same choice:

```python
from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod

engine = Vakint(
    rustred_reduction_strategy="ordered",
    evaluation_order=[VakintEvaluationMethod.new_rustred_method()],
)
```

The accepted strings are #raw("depth_first") and #raw("ordered"). Unknown
values produce an error. Ordinary numerical, symbolic and retained-master
evaluation retain their existing public interfaces.

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

The engine reuses checked rule selection, specialization, transports, nonzero
conditions and descent checks. Conditions from applied identities remain
required even when their contributions later cancel. Canceling a key before
applying its rule can avoid that rule and its conditions. Coverage warnings,
preserved unresolved integrals and numerical capability errors retain their
existing meaning.

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
The ordered strategy therefore has no universal speedup promise. Comparisons
must retain normal DFS request memoization and disable persistent caches for
both engines; cache-enabled historical results are a separate cohort.

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
  --features feynkit-ingress --example ordered_reduction_rungs -- \
  --work-dir /absolute/workspace/rung-L1-q0 \
  --loops 1 --projection zero-momentum --engine both --build-profile release
```

Through four loops, use separate #raw("--engine dfs") and
#raw("--engine ordered") invocations with fixed affinity and a release build
for fresh-process timing. Declare the actual build with #raw("--build-profile");
this is a label, not compiler introspection. The
#raw("--engine both") mode creates fresh engines in one process and is useful
for exact parity; it is not two fresh-process measurements. Existing report
files are never overwritten. Pass #raw("--artifacts PATH") to select an
authenticated checkout or installation; explicit four- and five-loop data is
needed. #raw("--prepare-only") saves input without reduction.
Five-loop physical runs use the ordered engine only, zero momentum followed by
quadratic. The accepted saved zero-momentum result supplies the exact historical
reference after conversion from its original pySecDec normalization through the
existing normalization API, before any Laurent truncation. Its raw payload and
original timing remain unchanged. No new five-loop depth-first run is required.

After separate runs, #raw("--work-dir OUT --compare-reports DFS_DIRECTORY ORDERED_DIRECTORY")
authenticates their saved payloads and compares exact retained expressions and
serialized numerical centers/error estimates without reducing again. These are
mathematical bounds only for components marked #raw("Bounded"). It writes a bound
#raw("parity.json") suitable for the benchmark adapter.

The retained timer covers the complete public retained evaluation, including
routing, scalar lowering and finalization. It is *not* a pure-IBP timer.
Data loading and the subsequent numerical reference calculation are separate.
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
