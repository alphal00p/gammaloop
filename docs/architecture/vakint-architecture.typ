= Vakint implementation architecture

#quote(block: true)[
#strong[Status:] Current implementation architecture, audited against the Vakint source on
2026-10-27.

This note describes the Rust engine and its optional Symbolica community-module wrapper. Backend
availability and numerical coverage depend on the selected topology, epsilon depth, and installed
external tools.
]

== Boundaries

Vakint is an in-process Symbolica transformation engine for single-scale vacuum integrals. The
Rust crate owns these boundaries:

- `Vakint` and `Topologies`: construction of the known topology library and entry points for
  canonicalization, tensor reduction, parametric evaluation, and numerical evaluation;
- `VakintExpression` and `VakintTerm`: the structured numerator/topology representation used
  between stages;
- `VakintSettings`, `EvaluationOrder`, and `EvaluationMethod`: precision, normalization,
  dependency, and backend-selection policy;
- FORM and pySecDec adapters: template rendering, subprocess execution, and result decoding;
- RustRed's immutable data registry, reduction programs, master capabilities, and strict caches;
- `precision` and `inspection`: tracked Laurent coefficients and snapshot-owned constant/master
  records, using Symbolica/Numerica's shared numerical error propagation;
- the optional Symbolica community Vakint module: Python classes that own a Rust engine and
  settings while exchanging Symbolica expressions.

The modules `alphaloop_numerics`, `matad`, `matad_numerics`, `fmft`, and `fmft_numerics` contain
backend-specific reductions and master-integral substitutions. `symbols` owns the namespaced
Symbolica symbols. `graph` is a small vacuum-topology graph used for validation, contraction, and
loop-momentum-basis work; it is not Linnet's generic half-edge graph.

== Core representation

`Vakint::new` initializes the symbol registry and generates an in-memory `Topologies` collection.
The collection contains canonical one- through four-loop families and selected contractions.
Explicitly loaded five-loop applications supply their own catalogue. An `Unknown` matcher permits
pySecDec when the caller allows unknown integrals. Each `Integral`
stores its canonical and optional short expressions, match patterns and conditions, graph,
loop/propagator counts, and applicable backend set.

Input is represented as a Symbolica `Atom`. `VakintExpression::split_integrals` separates a sum
into `VakintTerm` values, each containing:

- one `topo(...)` expression;
- its coefficient as the numerator;
- the loop and external vectors discovered in that numerator.

`ReplacementRules` is the result of matching one input topology. It carries the selected
canonical topology, the canonical-to-input edge map, substitutions for masses, powers, momenta,
and nodes, and the substitutions required to rewrite the numerator consistently. It is an
ephemeral transformation plan, not a persisted cache.

`VakintSettings` controls epsilon and scale symbols, executable paths, input rationalization and
runtime precision, numerator verification, normalization, unknown-integral policy, temporary
files, epsilon depth, dot-product notation, and `EvaluationOrder`. A `Vakint` engine is reusable;
settings are passed explicitly to Rust operations. The Python `Vakint` wrapper stores the engine
and settings together for convenience.

== RustRed data and numerical capabilities

Each engine owns an immutable `RustRedData` snapshot shared through `Arc`. Its programs, master
coefficients, topology identities, coverage predicates and graph references feed matching,
reduction, product factors, master evaluation and cache attachment. Ordinary construction remains
offline with cache-free one- through three-loop material. Four- and five-loop RustRed material
requires explicit loading.

`load_official_external_data(N, options)` validates and activates one cumulative one- through
N-loop snapshot. It recognizes a local checkout before considering network access; managed
installations use checksummed per-loop release archives, locking and atomic installation. Existing
validated packs are reused offline. Activation is transactional and replaces incompatible caches
and numerical context. Explicit IBP, saved-application and cache paths remain available. See
the crate's `external_data.typ` for loading and release conventions.

Prepared caches attach by default when the selected data supplies them. Official loading can
set `load_prepared_caches` to false to activate the validated registry without importing those
caches. A saved-application path inside a recognized release uses that release's cumulative
snapshot; a bare application supplies terminal declarations without inventing coefficient or
graph records.

RustRed master modes are `Retain`, `Symbolic` and `Numerical`, with numerical substitution the
default. Symbolic mode expands exact coefficients and retains labeled numerical-only coefficients;
the originating engine can resolve those symbols later. After exact combination, evaluation checks
only the required Laurent powers and numerical ingredients. Exact formulas have no finite source
cap unless they contain finite constants; numerical-only data carries validated significant-digit
limits. Missing coefficients are explicit errors. Finite data remains capped at its validated
precision; numerical evaluation reports precision shortfalls without inventing source digits.
Authenticated pack-level constant bindings carry finite values, significant digits and certified
absolute errors. Names such as `K5Moment(1)` and `K5LogMoment(3)` are finite sources, not exact
arbitrary-precision functions. Their payload and metadata contribute to the snapshot identity;
retained formulas use `RustRedConstant(snapshot, identity)` so another registry cannot silently
rebind the same name to different data.
Symbolic reduction binds these finite leaves after exact combination, using the same period
registry without choosing a numerical fallback. Loading rejects declarations that shadow native
evaluators, reserved Vakint symbols, existing period aliases or master basis identities.

Campaign certification is advisory outside its recorded predicates. Safe partial reductions retain
family-qualified unknown integrals. After exact combination, Retain mode warns with sorted,
deduplicated undeclared identities and authenticated candidate parent-topology graph references;
a common family does not imply a unique parent. Symbolic and Numerical substitution both fail
gracefully if surviving unknown integrals or declared terminals lack required coefficients.
Uncomputed declared pack masters remain valid retained results. Prepared caches accept strict
completed reductions, including valid captured descendants, and never store partial results.

Both FORM-backed and RustRed evaluation use the same tagged sixth-root HPL evaluator. FORM
reduction scripts and hard-coded master formulas remain intact. Opaque finite master data does
not become an arbitrary-precision constant through a symbolic name.
At endpoint one, the shared GPL adapter accepts words through weight six over zero and the sixth
roots of unity: `0`, `1`, `-1`, and `(±1 ± i√3)/2`. For a nonempty word, the first letter must
differ from one and the last must differ from zero. These use the ordinary `G` kernels
`dt/(t-a)` along the real path from zero to one. Exact algebraic expansion recognizes equivalent
letter spellings before constructing Symbolica's checked tags. The real and imaginary projection
wrappers preserve this restricted endpoint interface; general arguments, analytic continuation
and trailing-zero regularization are outside its contract. Constant discovery exposes concrete
real projections and loaded finite identities, while ordinary expression operators retain their
normal evaluation rules.
Before tracked evaluation, exact complex constant denominators use an unfactorized rational
polynomial over `Q(i)` so conjugate products cancel algebraically. This normalization preserves
finite source identities and skips raw floating coefficients; it does not discard uncertain
imaginary components or relax checked conversion to a real value.

== Main control flow

The complete Rust entry point is `Vakint::evaluate`:

```text
AtomView
  -> validate selected settings and external dependencies
  -> split into VakintExpression terms
  -> match and canonicalize topology plus numerator
  -> tensor-reduce each numerator through FeynKit by default, or selected FORM reduction
  -> choose the first supported EvaluationMethod for each term
  -> run the backend and decode its Symbolica expression
  -> combine RustRed contributions before master capability checks and substitution
  -> combine the transformed terms into an Atom
```

`to_canonical`, `tensor_reduce`, and `evaluate_integral` expose narrower stages. Canonicalization
tries registered topologies in order, derives replacement rules, normalizes momentum routing,
and can emit the short topology form. The FeynKit tensor reducer works in process with explicit
index scopes and exact coefficients. Before tensor collection, factorized Lorentz closure
validates original index multiplicities and contracts only connected raw-index components.
Closed loop invariants remain scalar coefficients through projection, preserving powers and
preventing spurious tensor-rank growth. Native tensor components retain their slots and pass
through the checked reducer; its projection, ambiguity and rank guards remain authoritative.
The optional FORM tensor path maps vectors and indices into
FORM syntax, runs the embedded reduction program, and maps the result back to Symbolica.
Saved index atoms retain the caller's epsilon symbol so vector mappings and returned
metadata agree. Literal native indices reuse exact registered Symbolica declarations;
conflicting attributes or tags still fail parsing. This transport does not alter FORM
programs, master formulas or scalar epsilon coefficients.

Integral evaluation matches again so that it has the canonical topology and numerator mapping.
It then walks `settings.evaluation_order` and selects the first method whose topology, loop count,
and requested epsilon depth are supported. AlphaLoop is limited to its registered topologies,
at most three loops, and at most four requested terms; MATAD covers its registered topologies up
to three loops and five terms; FMFT covers registered four-loop topologies up to five terms;
pySecDec covers topologies registered for that numerical path. RustRed uses loaded programs and
checks actual coefficient availability after reduction. Failure to find a method is an
explicit `NoEvaluationMethodFound` error, not an implicit fallback.

Optional request resources impose caller-selected process-RSS and elapsed-time ceilings;
omitted ceilings and count diagnostics remain disabled. Cooperative checkpoints cannot provide
a hard deadline inside an opaque call. `evaluate_with_monitor` uses the ordinary reduction path
and native observer, reporting phases, counters, bounded frontier samples, and independent
process and system memory observations. Any ETA remains a diagnostic estimate. Detaching the
terminal display does not cancel the calculation.

Five-loop retained evaluation also reports guard and coefficient materialization,
normalization and notation conversion, both outer and catalogue exact-combination
passes, and residual catalogue classification. These caller phases continue after
native reduction finishes; they inspect existing Atom headers rather than creating
another expression for monitoring. Each recorded stage duration is one completed
span, with repeated calls retained separately. `detailed_evaluation_total` includes
its reduction and materialization children plus callee cleanup;
`retained_finalization_total` includes the catalogue combination and classification
children. Add only the explicitly listed disjoint spans when comparing their sum
with public evaluation time. Never add a parent to its children or label an
unassigned residual as pure IBP time.

Symbolic sizes follow the objects actually stored by each owner. The caller reports the
encoded byte length and stored summands of its existing Atom at phase boundaries, without
formatting or constructing another expression. Native reports unstarted coalesced roots,
deferred rule maps and their child-integral occurrences, and collected output keys. Deferred
maps include the current map while its coefficients are being combined. These partitions
can overlap; they are not summed into a distinct-integral count or the size of an expanded
global expression, which the depth-first reducer never materializes. Native coefficient
byte censuses remain clone-owned bounds, distinct from Atom bytes and process RSS. Each
sample retains its scope and age; caller forwarding and fresh OS memory samples do not
make old native observations current. Metadata overflow gives an unknown count and never
introduces a reduction limit.

The terminal and supervisor snapshot viewer clamp scrolling to each panel's actual
wrapped content. Short panels remain visible while longer panels scroll; resizing or
shrinking observations recompute the navigation extent. Panel titles show their
current offset and final offset. This is display state only and cannot limit reduction.
OS memory is displayed in explicitly labeled IEC units; exact byte counts remain
available in the snapshot's `resources` object. Display rounding never changes a
reported guard threshold or its separately authenticated coverage. Restrained metric
accents distinguish current integrals, work counts, time and memory; unknown values
remain labeled and terminal color suppression leaves their content unchanged.

Numerical evaluation is a later boundary. Parametric backend output remains a Symbolica
expression. `precision` prepares native numerical sources and contracts external momenta before
delegating to Symbolica/Numerica's tracked evaluator. Each `NumericalCoefficient` owns its complex
central value and independent real/imaginary accuracy: significant digits, absolute error,
exact-zero status and bounded/estimated quality. `NumericalEvaluationResult` owns the ordered
coefficient map, requested decimal precision and structured source warnings. The optional second
result remains the backend's statistical uncertainty series, separate from arithmetic and
source-data uncertainty.

Raw coefficient records retain the evaluated midpoint and its enclosure. Checked symbolic
export enforces requested and achieved precision without lowering the raw midpoint first.
Inspection lookups explicitly cap their returned numerical values to the request and source
availability. Represented external dot products use the same component source identities and
reject incompatible dimensions or missing momenta.

Public MATAD/FMFT numerical master-substitution helpers prepare the same finite-source identities
and exact Gamma expressions used by the complete engine. Laurent expansion and normalization
precede tracked materialization; this preserves exact cancellations involving Gamma factors.
MATAD's existing exact master definitions can establish rational negative-order coefficients.
A numerical table entry is replaced by that exact value only when it matches nearest rounding
of the proved rational at its own stored precision. Nonrational poles, mismatches and nonpole
entries retain their finite-source identities. Rounded numerical agreement establishes no new
exact definition or source precision.

Evaluation uses one working-precision pass with the fixed guard allowance; it does not retry
after cancellation. Native `Float` and `ComplexFloat` inputs retain their stored precision, and
ordinary Python floats retain binary64 source limits. Exact integers and expression parameters
are substituted symbolically. Inexact expression leaves first receive immutable source identities,
so independent inputs with equal central values cannot cancel into certified exact zeros.
Authenticated runtime input tags retain their exact binary centers and original bit precision,
including one-to-four-bit cancellation residues with no complete decimal digits. They enter the
tracked binary-input path; decimal source caps for archived tables and external data are unchanged.
External components use private scalar identities after the common contraction boundary.
RustRed preserves finite input sources in all master modes. Input lifting happens before
topology collection and distribution: a shared factor keeps its identity, while independent
occurrences and evaluation requests have distinct identities. The same request counter serves
the lower public entry points, and already retained sources are left intact. Partial numerical
materialization also lifts raw finite literals before its first algebraic combination.
When a complete coefficient resolves to an exact real or Gaussian rational, partial materialization
retains that rational across later algebraic operations. Finite Float midpoints cannot enter this
shortcut, and final numerical-result export still uses checked precision-bearing Floats.
Laurent collection validates intermediate integer degrees before expansion and uses the configured
epsilon symbol throughout the Rust and Python entry points. Oversized degrees and unsupported
fractional or symbolic regulator powers produce descriptive errors instead of collector panics.
Parameters cannot shadow registered constants. Checked `to_atom` / `to_expression` conversions
encode each component's achieved precision and reject uncertain zeros or components without
reliably encodable digits; raw values and metadata remain inspectable.

`get_constant` and `get_master` inspect only the loaded or vendored RustRed registry. Typed
identities and discovery share its qualified names; unknown or ambiguous names return explicit
errors. Records own the original `Arc` snapshot and remain coherent after a release switch.
Master records expose every available epsilon coefficient, per-coefficient symbolic availability,
validated source precision and a native FeynKit diagram with authenticated powers, numerator and
loop-momentum routing. Declared uncomputed terminals remain inspectable with empty coefficient
maps. Inspection does not download data or use FORM's private tables.
Bare saved applications may announce master identities without authenticated graph payloads;
`get_master` then reports that the matching external-data pack is required.

== External tools and feature boundaries

AlphaLoop, MATAD, and FMFT invoke FORM. The pySecDec method invokes both FORM and pySecDec through
the configured Python executable. RustRed and the default FeynKit tensor path require no FORM or
Mathematica process. `validate_settings` checks only dependencies required by the
selected evaluation order, validates the loop-normalization expression, and enforces FORM >=
4.2.1 and pySecDec >= 1.6.4 when those tools are selected.

FORM programs, headers, and run templates are embedded in the Rust binary with `include_str!`.
Input serialization follows the exact expression tree, retaining symbol identities and argument
order instead of invoking display callbacks. Complex coefficients have explicit real and
imaginary components and parentheses, including a numeric coefficient for the imaginary unit.
For each run, Vakint renders those resources into a uniquely named temporary directory, launches
the configured process, reads `out.txt`, and removes the directory when `clean_tmp_dir` is true.
pySecDec can reuse a caller-selected output directory; the implementation labels this unsafe for
normal operation because stale generated sources or results can be mistaken for the current
calculation.

The Cargo build script has one platform-specific boundary: on macOS, setting
`EXTRA_MACOS_LIBS_FOR_GNU_GCC` opts into linking `gcc_s.1.1` from MacPorts' libgcc directory when
that library exists. This affects native linking only; it does not select an evaluation backend.

The core Rust API has no feature requirement. `symbolica_community_module` enables the PyO3
wrapper and registers the engine, expression, evaluation-method and numerical-result classes,
plus frozen inspection and precision records on the community module. `python_stubgen` adds stub metadata and
also enables Symbolica's Python export. The wrapper delegates transformations to the Rust engine
and converts `VakintError` into Python `ValueError`; it does not implement a second evaluation
pipeline. Its numerical input adapter preserves native precision and delegates contractions and
tracked evaluation to the same Rust core. One shared presentation model supplies structured
terminal and notebook views. Terminal output respects color settings; HTML reuses Symbolica's
formula rendering and FeynKit's existing diagram rendering. Display reads stored values only;
missing optional graph-rendering dependencies do not prevent metadata inspection.
Presentation hides transport tags when a formula belongs to one data snapshot, showing the
mathematical constant names while leaving the stored expression and provenance untouched.
Formulas containing several snapshots retain explicit tags to avoid apparent false cancellation.
The pySecDec adapter rejects retained uncertainty sources before generating input, because that
backend cannot propagate their source errors; ordinary numeric pySecDec inputs keep their existing
route.

== Persistence and reproducibility

Vakint does not serialize its topology library or settings. Known topologies are regenerated by
`Vakint::new`, which is why one engine should be reused for multiple operations. Durable values
cross the API as checked Symbolica expressions or precision-bearing Laurent records. Temporary FORM/pySecDec
workspaces are implementation artifacts unless cleanup is disabled or a reuse directory is
explicitly supplied. RustRed's external packs and native prepared caches are explicit durable
data, validated against their schemas and authenticated program bindings before attachment.

A running five-loop cache may append strict completed native Route/Apply nodes during ordinary
reductions. Explicit flush reports acknowledge durable records; freezing records a generation
descriptor for later read-only attachment. Loading an immutable snapshot has no writer and
captures no new records. Partial reductions never enter distributed strict caches.

Single-file cache admission uses `RustRedCacheFileOptions` (`CacheFileOptions` in Python) for
compressed, decoded and owned-index byte grants and an optional decoded-store directory. Omitted
grants retain the compact distribution defaults; larger independently hosted files require
explicit grants. These do not enlarge runtime RAM, work, queue or eviction policies. File-backed
loading authenticates transport and framing before attachment, shares a lazy coefficient table,
and checks complete record semantics at each solver hit.

Five-loop `FullSnapshot` export combines completed attached-base records and resident completed
nodes into one compressed native file with a shared coefficient table. An attached running
writer must be closed first, and the destination must be new. One- through four-loop exports
retain their native base-copy and checkpoint workflow. Python forwards the same file options
and snapshot operation to Rust; it does not implement another cache format.

A reproducible result therefore depends on more than the input expression: record the topology
form, numerator, normalization convention, epsilon depth, decimal precision, backend order and
options, numerical masses and external momenta, executable versions, and whether a pySecDec
workspace was reused. For RustRed, also record the resolved data version and load report.

== Maintained invariants

- Each input term contains exactly one first-power `topo(...)` factor; scalar-only terms and
  powers of a topology are rejected by the structured-expression boundary.
- Propagator ids are contiguous from one, and loop-momentum ids are contiguous from one through
  the detected loop count. Vacuum-graph nodes must have degree greater than one.
- Canonical topology substitutions and numerator substitutions are applied as one mapping;
  sequential replacement is avoided because one replacement target can be another source.
- When `verify_numerator_identification` is enabled, no unrecognized loop-momentum dependence may
  remain after canonical mapping.
- `EvaluationOrder` is policy: the first supported method wins. Unsupported methods are skipped,
  but execution errors from the selected method are returned.
- Loop normalization may reference only the documented symbols and is expanded to the requested
  epsilon order before backend results are normalized.
- Requested working precision never increases the declared accuracy of an input or finite source.
  Exact cancellation and a numerically unresolved zero have distinct metadata and conversion rules.

== Verification and related documentation

The integration-test modules separate topology/input matching, tensor reduction, analytic
evaluation, free-form evaluation, backend comparisons, and direct pySecDec evaluation. Analytic
cross-checks compare AlphaLoop with MATAD, while numerical suites compare supported analytic
results with pySecDec. Tests that need pySecDec detect or explicitly skip that dependency; FORM
and backend-enabled coverage still requires the corresponding external tools. The crate also has
focused unit tests for expression encoding, source limits, component accuracy, cancellation and
checked conversion. The Python precision/inspection suite verifies native high-precision inputs,
exact and inexact external kinematics, frozen records, typed lookup failures and native graph/HTML
rendering against one coherent Symbolica extension.

For supported workflows, see the
#link("../../../products/vakint/latest/tutorial/")[matching tutorial],
#link("../../../products/vakint/latest/guides/evaluation/")[evaluation and reproducibility guide],
#link("../../../products/vakint/latest/reference/interfaces/")[Rust/Python interface guide], and
#link("../../../products/vakint/latest/reference/topologies/")[generated topology table]. Exact
Rust signatures are in the
#link("../../../products/vakint/latest/reference/rust/vakint/")[native Vakint Rustdoc].
