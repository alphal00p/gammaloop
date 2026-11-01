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
reduction operations and Vakint's only reduction strategy: depth-first
reduction, its completed-result caches and their Python operations were removed
from Vakint. The RustRed library keeps its depth-first reducer for its own users.

#panel[
  *Scope.* Existing reduction tables and topology matching; exact rational
  coefficients; a #raw("BTreeMap") frontier; no reduction-result cache, only a
  cache of one-step records (see "Ordered step cache").
  No offline rule rewriting, synthesized composite rules, or loop-specific
  scheduling heuristics are required.
]

== Configuration

There is no strategy setting: every RustRed evaluation, at one to five loops,
uses the ordered scheduler. #raw("VakintSettings::rustred_application_options")
is RustRed's #raw("SavedApplicationOptions"); Vakint sets its
#raw("reduction_strategy") to #raw("Ordered") by default and ignores any other
value, as it ignores the depth-first #raw("cache_policy") and
#raw("cache_representation"). Its reduction limits, preparation workers and
#raw("ordered_weight_arithmetic") (see below) apply. The scheduler changes the
order of the checked identities, not the definitions, normalization or
requested numerical precision of the masters. RustRed leads the default
evaluation order; FORM-based backends follow it only when FORM is usable.

The Python constructor has no strategy argument. The former
#raw("rustred_reduction_strategy") and #raw("rustred_cache_policy") keywords and
the depth-first cache operations (#raw("prepare_cache"), #raw("load_cache"),
running caches, #raw("save_cache"), #raw("import_cache"),
#raw("download_cache"), #raw("cache_statistics") and the related clear/detach
methods) were removed. #raw("load_official_external_data") no longer takes
#raw("load_prepared_caches"). #raw("VakintSettings::ordered_cache") (Python:
the #raw("ordered_cache*") keywords) configures the Ordered step cache below.

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
The zero entry disappears, leaving only $V+W$ to reduce. The removed
depth-first reduction expanded each parent to its terminals separately, so it
generally reduced $U$ before the two contributions could cancel.

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
- #raw("modular_polynomial") records the same trace and replays it once per
  64-bit prime with weights that are rational functions of the dimension over
  that prime field, without sample points or interpolation. It is the fastest
  arithmetic: at five loops (quadratic projection) it needs about 28 times less
  replay CPU than #raw("modular_sampled").

Both modular modes lift their prime images over the rationals (Chinese
remaindering and rational reconstruction) and verify every coefficient on fresh
primes. They need a coefficient context with at most one variable, which holds
for Vakint's native requests. Every mode expands the same keys and returns the
same terms and nonzero conditions; only the monitor's frontier monomial counters
differ (modular weights count one monomial per key). The default is
#raw("vakint::DEFAULT_ORDERED_WEIGHT_ARITHMETIC")
(#raw("modular_polynomial")), a single constant that also sets the Python
default; #raw("exact") and #raw("modular_sampled") stay selectable.

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
#raw("threads"). The same count sets the worker threads of RustRed's Ordered
scheduler for four- and five-loop requests (closed one- to three-loop programs
always use the sequential scheduler), the five-loop frontend workers (witness
screening and per-contribution scalar lowering), the retained post-processing
workers (guards, coefficients and the master catalogue combination) and the
workers verifying the five-loop routes while the saved application loads
(#raw("preparation_workers"), derived and otherwise left at its default; data
loaders take it from #raw("OfficialExternalDataOptions::n_cores")). Results are
identical for every core count.

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
  [Step cache hits / misses], [Cache panel: expanded keys of this request whose
    checked step was replayed from the Ordered step cache / ran (and was
    recorded when the cache populates), for the request, the current phase and
    the last completed phase. The JSON snapshot lists the most recent completed
    phases under #raw("ordered_reduction.step_cache"). Both stay zero without
    a cache.],
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
application path), the strategy (always #raw("Ordered")), weight arithmetic, master mode and
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

== Ordered step cache

Every expanded key runs one checked step: the selected route or rule, the
ordered children with their exact step coefficients, emitted terminal
coefficients and the step's nonzero conditions. The step is a deterministic
function of the key and the immutable program, so RustRed's
#raw("OrderedStepCache") records it once and replays it on every later
expansion of that key. A replayed step is identical to the checked one: the
scheduler performs the same merges, cancellations and guard interning as
without a cache, in every weight arithmetic and every modular pass. Records
never contain frontier or terminal weights, zero decisions or final answers,
so cold, warm and uncached runs have identical frontier digests at every
phase end and identical results. Closed-artifact product factorizations are
recorded as opaque and always recomputed, because their memoized lower-factor
reductions retain guards only at first use.

Vakint keeps one cache per reduction program: the closed one- to three-loop
artifacts, the four four-loop candidate programs and the five-loop saved
application. Records are sharded by stage, support and owner sector (Route by
the owner it transports to, Apply by owner) and are shared by all requests,
channels and arithmetics of an engine and its clones.
#raw("VakintSettings::ordered_cache") selects how evaluations use them:

#table(
  columns: (0.3fr, 1fr), inset: 7pt, stroke: 0.4pt + rgb("d5dee7"),
  fill: (x, y) => if y == 0 { pale } else { none },
  table.header([*Field (Python keyword)*], [*Meaning*]),
  [#raw("mode") (#raw("ordered_cache"))], [#raw("read-write") (default):
    replay recorded steps and record new ones (live population).
    #raw("read-only"): replay the steps saved in the directory (without a
    directory: the shipped caches only), record nothing. #raw("off"): no
    cache.],
  [#raw("directory") (#raw("ordered_cache_directory"))], [Root directory; none
    (default) keeps records in memory for the lifetime of the engine. Shard files
    are read at their first lookup.],
  [#raw("max_resident_bytes") (#raw("ordered_cache_ram_bytes"))], [Ceiling on
    the estimated resident bytes of each program's cache (default 8 GiB). New
    records and shards beyond it are not kept and their steps are recomputed;
    the Symbolica RSS budget of #raw("rustred_resources") applies as well.],
  [#raw("autosave") (#raw("ordered_cache_autosave"))], [When records of a
    directory-backed read-write cache are written: #raw("after-evaluation")
    (default; at the end of every public evaluation call, also a failed one),
    #raw("on-drop") (when the last clone of the engine is dropped) or
    #raw("never") (only explicit saves). An automatic save never fails an
    evaluation; its last error is reported in the statistics.],
  [#raw("shipped") (Rust only; default on)], [Use the caches the loaded
    external data ship (see "Shipped caches" below); #raw("false") ignores them,
    e.g. for a cold preparation.],
)

Operations (Rust #raw("Vakint") methods, Python methods of the same name):
#raw("load_ordered_cache(settings)") attaches the caches of every loaded
program now (a five-loop application once it is loaded) and checks their
manifests; #raw("save_ordered_cache()") writes the unsaved records of every
directory-backed read-write cache; #raw("ordered_cache_statistics()") reports
hits, misses, inserts, records, resident and file bytes in total and per
program and shard; #raw("clear_ordered_cache()") drops the resident records
(files are kept); #raw("warm_ordered_cache(settings, integrals)") evaluates each
expression with RustRed in retained mode (results discarded, the cache used
read-write) and saves. The rung tests and examples read
#raw("VAKINT_ORDERED_CACHE") (#raw("off"), #raw("read-only"),
#raw("read-write")), #raw("VAKINT_ORDERED_CACHE_DIR"),
#raw("VAKINT_ORDERED_CACHE_SHIPPED") (#raw("on")/#raw("off")) and
#raw("VAKINT_ORDERED_CACHE_MAX_RESIDENT_BYTES"), and report the cache counters
(with each program's shipped-cache status) in #raw("ordered_step_cache").

A cache directory holds one subdirectory per program,
#raw("<loops>L-<program>.<tables>.<stamp>") (e.g.
#raw("4L-four-loop-h.5c1d....1f2e...")): a manifest and content-addressed shard
files written atomically, so an interrupted save keeps the previous cache
valid. Validity is decided by versions and cheap identities only, never by
digests of engine sources. The stamp is a BLAKE3 digest of RustRed's binding:
its cache schema version, its hand-bumped step-semantics version and cheap
identity fields (adapter, family fingerprint, coefficient variables,
unresolved policy). The tables identify the program's rule tables with digests
the loader has already verified, so no data file is hashed for the cache: the
sha256 of a 1L–3L #raw("ibps/rules.rrbin"), the declared sha256 values of a 4L
program's four #raw("ibps") files, and the identity of the 5L saved application
(its #raw("application.json"), which pins every payload, with the pack's
families and masters). The external-data release version is reported
(#raw("data_version") in the statistics) but is not part of the key, because
checkouts change their tables without changing it. New rule tables (a data
update or a local override) or a RustRed semantics bump therefore write a new
subdirectory and never replay foreign records, while a data update that keeps a
program's tables keeps its records; stale subdirectories can be deleted. When
loaded data or a five-loop application replace a program, its cache is detached
(saved first unless its autosave policy is #raw("never")). RustRed rejects a
manifest or shard whose stamp or payload digest differs; Vakint reports this as
an error naming the directory. A directory
has one writer at a time: two engines or processes writing the same
directory lose records of one of them (never mix them); readers are safe.

=== Shipped caches

External data may ship the cache of each program of a loop pack in the same
layout: #raw("<L>L/ordered-cache/<label>.<tables>.<stamp>/") (the RustRed
manifest, the shards it lists and optionally the tail table
#raw("ordered-tail.rrtail")), described by #raw("<L>L/ordered-cache.json")
(schema #raw("vakint.external-data.ordered-cache.v2"): producer, workload, and
per program its label, directory, rule tables, stamp, records and bytes, and
whether it is complete or trimmed) and listed file by file with their sha256 in
the loop manifest (section #raw("ordered_cache")). Loading does not read them.
The same versioning rule applies as to every cache: a shipped cache is used iff
it was recorded for the program's rule tables under this build's stamp (step
cache schema and semantics version and the binding's identity fields; a tail
table must carry this build's tail schema and semantics) and its files match
the loop manifest's sha256. Otherwise Vakint prints one line on stderr
(#raw("vakint: shipped Ordered step cache ... not used (<reason>); <program>
runs without it")) and the program runs cold. There is no fingerprint,
attestation or restamping: a RustRed release that bumps the step semantics
needs regenerated caches, any other release keeps them.

Without a cache directory (the default) an accepted cache is opened in place:
its shards are read at their first lookup, new records stay in memory, and
nothing is ever written into the data directory; #raw("read-only") replays the
shipped records only. With a read-write cache directory the files are copied
into the program's subdirectory there if it holds no cache yet (an existing
cache is kept) and the cache then grows there as usual; #raw("read-only")
never writes, so it opens the shipped cache in place when the directory holds
no cache of the program. Each program's statistics report
the outcome (#raw("shipped"): #raw("in-place"), #raw("installed"),
#raw("already-present") or #raw("rejected") with the reason);
#raw("shipped_ordered_caches()") lists every loaded program with its
directory name, tables, stamp and the verified state of its shipped cache.
The vendored 1L–3L data ship no caches (their cold reduction takes well under
a second).

The example #raw("ordered_cache_maintenance") maintains them. #raw("prepare")
runs a deterministic workload cold (shipped caches ignored) into a fresh cache
directory: the zero-momentum and quadratic QCD-ladder rungs at one to four
loops, the four-loop benchmark integrals and the five-loop benchmarks of the
data, then exports the caches of the loaded programs in the pack layout
(#raw("export_ordered_caches"); the bytes depend only on the records, so the
export is reproducible) to be staged by vakint-external-data's tooling.
#raw("check") prints, for every program loaded from a data directory, its
directory name, tables, stamp and whether its shipped cache is accepted, and
fails if one is rejected. The five-loop zero-momentum and quadratic workloads
reduce the same saved application, so they share one program cache; a pack
may ship it trimmed (#raw("complete: false")) when the full quadratic cache
exceeds the size budget, and the full cache can be placed unchanged into a
cache directory under its directory name.

The depth-first caches stored a complete terminal decomposition per integral.
Consuming such a record replaces a frontier key by its terminal image and skips
every grouping and cancellation below it, so ordered reduction never used them.
External-data releases up to v0.1.2 still ship them; Vakint neither loads nor
checksums them (an explicit read of such a file still verifies it), which
removes a large part of the five-loop loading time.

== Saved five-loop route preparation

Loading the five-loop saved application verifies its 8,179 momentum routes
exactly: about 136 s on one core and about 64 s at any larger
#raw("n_cores") (RustRed's verification scales only to about two workers).
That work depends only on the saved application in #raw("5L/ibps"), so RustRed
records it once into a route preparation add-on
(#raw("rustred_app::write_saved_application_route_addon")) and restores the
identical routes from it in under a second; the whole native load then takes
about 22 s at every core count. The add-on is an extra file next to the
unchanged artifact: a data pack declares it in #raw("5L/manifest.json") as

```json
"route_preparation": {"path": "route-preparation.rrbin",
                      "format": "rustred.saved-application.route-preparation.v1"}
```

and lists the file in #raw("files"), so its sha256 is checked like every other
file. Vakint uses it only for the official pack's own #raw("ibps/")
application; RustRed accepts it only when its format, its verified-map
semantics version and the application identity (family fingerprint and the
digest of #raw("selection.json")) match. A missing, corrupt or stale add-on is
rejected with a warning and every route is verified, as without it. A user
artifact passed through #raw("rustred_artifact_path") uses RustRed's default: a
#raw("route-preparation.rrbin") next to its #raw("application.json"), when
present. Produce or refresh the add-on with

```sh
cargo run --release -p vakint --example five_loop_route_preparation -- \
    DATA_ROOT route-preparation.rrbin 6 --check
```

which verifies every route once, writes the file, checks that loading from it
gives identical routes and prints the manifest entries. Regenerate it whenever
the 5L #raw("ibps/") files or RustRed's #raw("SAVED_VERIFIED_MAP_SEMANTICS")
change.

== Limitations

Keeping more unfinished contributions visible can increase frontier memory.
Exact coefficient combination can also cost more than the reductions it avoids;
the step cache saves the checked steps, not the frontier arithmetic (a warm
exact run is dominated by the latter). Ordered reduction therefore has no
universal speedup promise over the historical depth-first measurements, which
remain separate cohorts in the benchmark report.

== Reproducible physical comparisons

The collector #raw("examples/ordered_reduction_rungs.rs") uses the same native
QCD graph, color/Dirac contraction, routing and scalar-lowering fixture for both
engines. All new rung fixtures, numerical references and measurements use
#raw("LoopNormalizationFactor::MSbar"). One loop has no vertical gluon rungs;
each successive loop adds one.
Both zero-momentum metric and complete quadratic transverse projections are
tested. Prepared scalar inputs are authenticated and reused independently of
reduction, preserving full mass and numerator conventions.

From the repository workspace root, for an exact/modular functional comparison:

```sh
cargo run --locked --release -p vakint \
  --features feynkit-ingress,experimental-rustred --example ordered_reduction_rungs -- \
  --work-dir /absolute/workspace/rung-L1-q0 \
  --loops 1 --projection zero-momentum --arithmetic both --build-profile release
```

#raw("--arithmetic") selects #raw("exact"), #raw("modular[:key=value,...]") or
#raw("both"); without it the collector uses
#raw("VAKINT_ORDERED_WEIGHT_ARITHMETIC") (Vakint's default when unset). #raw("both")
creates two fresh engines in one process and checks that they retain the same
expression; it is not two fresh-process measurements. Use separate invocations
with fixed affinity, #raw("--no-monitor") and a release build for qualified
fresh-process timing. Declare the actual build with #raw("--build-profile");
this is a label, not compiler introspection. Existing report files are never
overwritten. Pass #raw("--artifacts PATH") to select an authenticated checkout
or installation; explicit four- and five-loop data is needed.
#raw("--prepare-only") saves input without reduction. Five-loop physical runs
reduce zero momentum followed by quadratic.

The historical depth-first results remain the exact references. The rung tests
compare the exact Ordered result of every one- to four-loop case with the
depth-first retained expression frozen in
#raw("tests/resources/ordered_rungs/<case>.depth-first-retained.symbolica")
(produced at gammaloop 31a8b166e) and with the recorded high-precision
numerical reference, and require modular weights to retain the identical
expression. The accepted saved five-loop zero-momentum depth-first result
supplies the five-loop reference after conversion from its original pySecDec
normalization through the existing normalization API, before any Laurent
truncation.

#raw("--work-dir OUT --compare-two-reports REFERENCE_DIRECTORY ORDERED_DIRECTORY")
authenticates two saved reports, an Ordered run and either another Ordered run
or a preserved historical depth-first report, and compares their exact retained
expressions and, unless #raw("--no-numerical-reference") is given, serialized
numerical centers/error estimates without reducing again. These are
mathematical bounds only for components marked #raw("Bounded").

The comparison headline is #raw("reduction_math_seconds"): measured
non-overlapping math spans covering IBP/routing and exact weighted terminal
coefficient assembly/grouping. Necessary terminal substitutions, projections
and lower-factor reductions are included, as is request-local reducer/dispatch
state construction. Immutable artifact loading, deserialization
and program import, graph/data preparation, matching,
initial scalar lowering, MATAD normalization, numerical evaluation and
serialization are excluded. This is broader than a Native-kernel IBP timer.
It sums disjoint math blocks, including later exact coefficient combination;
normalization intervals themselves remain outside. Natural request-state cleanup
inside the current call is included, while engine-lifetime teardown after the
call is excluded.

The collector retains gross math time and subtracts the *measured* lazy owner
preparation wall time, with one owner-preparation worker. The checked net,
interval count, fixed scope and completeness flag remain in the report. Invalid
telemetry makes the net unavailable without invalidating the scientific result.
The optional capture requires #raw("experimental-rustred"); a build without it
reports no math duration. The public #raw("retained_reduction_seconds") timer
remains separate and includes frontend lowering/finalization. Missing math
times never fall back to that public timer. In particular, the historical
40-hour public observation supplies no isolated reduction-math duration.

No prepared or previous-run result state is loaded for these comparisons; each
report records #raw("intrinsic_memoization") as false. Ordered's opaque checked
factorization may use a fresh private memo internally; it does not feed a
reusable result cache. Reusing the common authenticated scalar input is frontend
preparation, not a reduction-result cache. Data loading and the subsequent
numerical reference calculation remain separate. Qualified comparisons disable
monitoring. Dashboard/snapshot checks, including the large five-loop runs,
retain their observer overhead and are diagnostic observations; independent
RAM guarding remains active.
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
