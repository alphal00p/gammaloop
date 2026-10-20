= Vakint implementation architecture

#quote(block: true)[
#strong[Status:] Current implementation architecture, audited against the Vakint source on
2026-10-19.

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

RustRed master modes are `Retain`, `Symbolic` and `Numerical`, with numerical substitution the
default. Symbolic mode expands exact coefficients and retains labeled numerical-only coefficients;
the originating engine can resolve those symbols later. After exact combination, evaluation checks
only the required Laurent powers and numerical ingredients. Exact formulas have no finite source
cap unless they contain finite constants; numerical-only data carries validated significant-digit
limits. Missing coefficients and insufficient precision are explicit errors.

Campaign certification is advisory outside its recorded predicates. Safe partial reductions retain
family-qualified unknown integrals and warn with parent-topology graph references. Numerical
evaluation lists unresolved masters and uncomputed declared terminals. Prepared caches accept
strict completed reductions, including valid captured descendants, and never store partial results.

Both FORM-backed and RustRed evaluation use the same tagged sixth-root HPL evaluator. FORM
reduction scripts and hard-coded master formulas remain intact. Opaque finite master data does
not become an arbitrary-precision constant through a symbolic name.

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

Integral evaluation matches again so that it has the canonical topology and numerator mapping.
It then walks `settings.evaluation_order` and selects the first method whose topology, loop count,
and requested epsilon depth are supported. AlphaLoop is limited to its registered topologies,
at most three loops, and at most four requested terms; MATAD covers its registered topologies up
to three loops and five terms; FMFT covers registered four-loop topologies up to five terms;
pySecDec covers topologies registered for that numerical path. RustRed uses loaded programs and
checks actual coefficient availability after reduction. Failure to find a method is an
explicit `NoEvaluationMethodFound` error, not an implicit fallback.

Numerical evaluation is a later boundary. Parametric backend output remains a Symbolica
expression. Numerical parameters and optional external momenta are converted to the configured
binary precision, and `NumericalEvaluationResult` represents the result and optional error as
Laurent coefficients `(epsilon power, complex value)`.

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
wrapper and registers `Vakint`, `VakintExpression`, `VakintEvaluationMethod`, and
`VakintNumericalResult` on the community module. `python_stubgen` adds stub metadata and
also enables Symbolica's Python export. The wrapper delegates transformations to the Rust engine
and converts `VakintError` into Python `ValueError`; it does not implement a second evaluation
pipeline.

== Persistence and reproducibility

Vakint does not serialize its topology library or settings. Known topologies are regenerated by
`Vakint::new`, which is why one engine should be reused for multiple operations. Durable values
cross the API as Symbolica expressions or Laurent-coefficient lists. Temporary FORM/pySecDec
workspaces are implementation artifacts unless cleanup is disabled or a reuse directory is
explicitly supplied. RustRed's external packs and native prepared caches are explicit durable
data, validated against their schemas and authenticated program bindings before attachment.

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
- Runtime numerical precision is converted to a binary precision of at least the equivalent of
  17 decimal digits; Python convenience inputs are converted at that boundary.

== Verification and related documentation

The integration-test modules separate topology/input matching, tensor reduction, analytic
evaluation, free-form evaluation, backend comparisons, and direct pySecDec evaluation. Analytic
cross-checks compare AlphaLoop with MATAD, while numerical suites compare supported analytic
results with pySecDec. Tests that need pySecDec detect or explicitly skip that dependency; FORM
and backend-enabled coverage still requires the corresponding external tools. The crate also has
focused unit tests for expression encoding and numerical-result behavior.

For supported workflows, see the
#link("../../../products/vakint/latest/tutorial/")[matching tutorial],
#link("../../../products/vakint/latest/guides/evaluation/")[evaluation and reproducibility guide],
#link("../../../products/vakint/latest/reference/interfaces/")[Rust/Python interface guide], and
#link("../../../products/vakint/latest/reference/topologies/")[generated topology table]. Exact
Rust signatures are in the
#link("../../../products/vakint/latest/reference/rust/vakint/")[native Vakint Rustdoc].
