= Precision and inspection validation

This report records completed checks for the precision-aware evaluation and
constant/master inspection changes. The implementation is still being
validated; the pending gates below have not been claimed as complete.

== Numerical contract

Each complete Laurent coefficient is evaluated once with Numerica complex
balls through Symbolica's direct interpreter. Exact rational algebra precedes
approximation. Finite master values, floating inputs, masses and normalization
factors retain their source identities and uncertainty through that boundary.
The result records real and imaginary accuracy separately, including absolute
error, exact-zero status and whether the error is bounded or estimated.

Working precision does not enlarge a finite source's validated capability.
Warnings identify limiting sources, master coefficients, epsilon orders,
requested/achieved accuracy and data versions. Checked symbolic conversion
encodes each component's achieved precision; an uncertain zero remains
inspectable but cannot be exported as an exact zero.

== Completed implementation checks

#table(
  columns: (2fr, 1fr, 3fr),
  [Check], [Result], [Scope],
  [Optimized Vakint units], [313 passed],
  [All enabled library tests with Python, inspection, experimental RustRed and
    weighted diagnostics; 20 explicitly skipped, no retries, 144.535 seconds.
    Includes finite primitive admission, symbolic snapshot isolation, low-bit
    input accuracy, exact partial coefficients, checked Laurent orders,
    stale-release fixture reuse and malformed-source error handling.],
  [Binary inputs and nested integration], [44 + 4 passed],
  [44 precision tests and four core integration tests over eight graphs. Original
    bit precision survives low-precision Float residues; exact rational partial
    coefficients survive nested forest algebra. Literal physics references,
    forest counts and existing numerical tolerances are unchanged.],
  [Precision and inspection], [60 passed],
  [Process-isolated numerical, lookup, input-identity and public API tests.],
  [Numerical result records], [6 passed],
  [Finite imports, explicit precision suffixes, repeated Laurent orders and
    independent statistical-error quadrature.],
  [Four-loop scientific inventory], [15 cases passed],
  [H/FG/BMW/X independent FORM reductions versus current-data RustRed at the
    original `1e-30` tolerance. Historical FORM table-center references retain
    their strict check, with a separate tracked-error containment check.],
  [FORM-free rank-four input], [Passed],
  [Original numerical reference and `1e-30` tolerance, with the authored binary
    reference-point coordinates preserved exactly.],
  [Legacy reference migration], [54 passed],
  [Original reference values, relative tolerances, epsilon targets and reduction
    workloads retained; no retries, with the original 240-second limits. Authored
    binary reference coordinates are reconstructed exactly in test fixtures.],
  [Python bindings], [17 passed],
  [Native high-precision inputs, exact and uncertain cancellation, immutable
    records, graceful lookup failures, finite moment bindings, snapshot isolation
    and reserved-constant validation.],
  [Complete inspection inventory], [1,132 values passed],
  [All 115 computed masters and their 739 coefficients, plus all 393 discovered
    constants at 16 digits; zero warnings. All 378 named real constants used by
    677 exact source formulas are represented. Equivalent-master provenance is checked.],
  [FeynKit graphs], [76 passed],
  [DOT/model semantics, routing, loop bases, half-edge order, dotted
    propagators and numerator annotations.],
  [CFF], [34 passed],
  [Includes explicit rejection of unsupported propagator powers.],
  [Master graph inventory], [1,020 passed],
  [Full native semantic DOT round trips for 115 computed master entries and
    905 declared, uncomputed five-loop terminals.],
)

The focused suites overlap and must not be added into a whole-repository test
count. Vakint passes Clippy for all targets with the inspection, Python,
experimental RustRed and weighted-diagnostic features enabled. FeynKit graph
and CFF Clippy checks also pass with warnings denied.

The low-level MATAD and FMFT master-preparation helpers preserve exact algebra
until the common tracked evaluation boundary. MATAD promotes an archived
negative-order coefficient to a rational only when its existing exact formula
simplifies to that rational and its stored Float equals the correctly rounded
value at the original precision. Finite parts, nonrational poles and mismatches
retain their source uncertainty. The original exact Gamma-series assertions and
the convergent two-loop MATAD/RustRed comparison pass without relaxing their
checks.

Runtime input tags preserve exact binary centers and original bit precision,
including one-to-four-bit cancellation residues with no complete decimal digits.
They use the tracked binary-input path without padding a decimal source cap.
Partial materialization preserves resolved exact real and Gaussian rationals;
finite Float values remain identifiable sources, and final numerical records
still export checked Floats. The focused core gate preserves all 92 original
literal references and the six two-loop forest counts.

The dependency integration changed printer parentheses and canonical cyclic
gamma representatives in 16 tensor snapshots. Seventeen exact old/new formula
comparisons verify mathematical equality. Updated expectations and a direct
disabled-trace assertion pass the complete Nix Spenso group: 578 tests, with
no production algebra change.

GammaLoop's generated evaluators now use native pi constants at the active
precision. Three internal f64 overrides previously disagreed with Vakint's
high-precision pi coefficients. User-supplied f64 points and settings keep their
existing interpretation. State layouts are unchanged; evaluator trees already
saved in `gammaloop_state/` retain their embedded constants and must be regenerated
to obtain the corrected Quad/Arb evaluation.
The corrected scalar-bubble values are independently derived from its closed
CFF expression using bounded rational Machin pi and 380-digit decimal sqrt/log
arithmetic. Quad and Arb comparisons retain their original `1e-27` and `1e-60`
limits. The final strict bubble test passes in 1.82 seconds. The GL20 snapshot
retains its exact assertion, metadata and kinematics; 61 measured numerical
leaves are refreshed after a separate 1,000-bit replay and independent audit.
Its final strict test passes in 5.79 seconds.

The package inventory contains cache-free 1L–3L runtime data. It contains no
4L/5L pack, prepared cache or former artifacts directory. Ordinary engine
construction remains offline; explicit external-data loading supersedes the
vendored loops.

Cargo's package inventory check passes. Creating a crates.io archive remains
unsupported by the existing Git-only RustRed dependency, which has no registry
version requirement; the repository build uses its authenticated commit pin.
No registry version was invented to bypass that packaging check.

The optimized prepared-cache parity, four-loop factor-provider and portable
five-loop application tests pass with their original 240-second limits and no
retries (27.83 s, 33.34 s and 90.24 s respectively). The underlying 9,969
IBP/cache/control files are unchanged; the complete cache workloads contain
5, 77, 2,320 and 164,352 roots through four loops.

The refreshed reduction benchmark uses optimization level three and fat LTO,
48 fresh-process samples, exact rational parity and verified cache hits. Cached
RustRed is 5.81–16.60 times faster than FORM/FMFT for the four prepared base
cases and 1.42–3.92 times faster for the held-out cases. Engine/data loading
and cache attachment are outside the timer. The canonical one-page benchmark,
Typst source and raw JSON record the measured source, data and environment.

== Interfaces and presentation

Rust and Python expose typed constant/master identities, available source
precision, achieved numerical accuracy and provenance. Master inspection
returns every recorded Laurent order and the shared FeynKit diagram; declared
but uncomputed terminals have empty coefficient maps. Unknown or ambiguous
identities return descriptive errors.

All 19 Python exports and 180 members match the generated canonical stub. Documentation
reference generation and architecture-owner checks pass. Actual Marimo output
was rendered and visually inspected with native Symbolica formulas and
FeynKit's SVG graph hook. Rendering uses stored values without evaluation or
downloads.

The final Python host also passes 13 Laurent-boundary controls using a custom
epsilon symbol: five valid cases include orders 128 and -129 and exact rational
cancellation; eight oversized, fractional or symbolic-power cases return
descriptive errors before expansion. No panic is used as an error path.

== Dependency integration

Symbolica, Numerica and Graphica share the published integration revision
`9a3abbeaac5c6499048bc4600a348762f356cec9`; no development path override is
required. The coordinated tree passes 892 Symbolica tests, with four ignored.
The precision-preservation fix and tracked evaluator are submitted as
#link("https://github.com/symbolica-dev/symbolica/pull/112")[Symbolica PR 112]
and #link("https://github.com/symbolica-dev/symbolica/pull/113")[PR 113].

RustRed remains at the compatible owned-program/partial-reduction revision
`a666e705a92b06ae640a0d9e69b2dfe5bbe748ef`. Its guarded cache, persistence and
partial-reduction checks pass 129 distinct tests with the new Symbolica cohort,
including ten additional durable-codec checks. Native codec/export/import
versions and the compatible distributed cache payloads are unchanged.
Focused example-build fixes against current RustRed main are submitted in
#link("https://github.com/alphal00p/rustred/pull/5")[RustRed PR 5].

The restricted endpoint-one GPL extension passes 22 independent component
comparisons at 385 digits and a separate 20,000-digit gate. All 13 newly
discovered GPL projections expose arbitrary precision. The complete PR14
epsilon-one coefficients for both P31 and P58 achieve 240 bounded digits and
agree with their archived 200-digit numerical references. Exact denominator
algebra establishes their realness before tracked evaluation; no uncertain
imaginary component is discarded.

== Finite periods and scientific evidence

PR4fin and PR4dfin, together with their first epsilon coefficients, have
independently validated 20,000-digit payloads. Separate 20,200-digit production
inputs are preserved for the remaining period calculation. Their proof uses
independent Bessel-moment partitions, rigorous truncation bounds, independent
exact rational reconstruction and Mathematica development references. A
transparent additive receipt corrects the original reference-anchor inventory;
the original sealed payloads and certificates remain unchanged.

New upstream exact PR4-family formulas are preserved. Their named
`K5Moment(1|3)` and `K5LogMoment(1|3)` dependencies bind directly to independently
certified finite moment payloads in the active snapshot. The logarithmic
normalization and all 24 new formulas were checked independently; no circular
master inversion or process-global callback is used to supply these values.

Public constant and master lookups at 20,000 and 20,001 requested digits pass
exact rational certificate-containment and identity checks for the finite
moments and all three primitive periods, including the three direct PR11d
aliases. Derived PR11d formulas inherit the 20,000-digit capability while
epsilon-one coefficients retain their 200-digit limits. This final gate passes
in 19.82 seconds. The current
conservative source/error conversion reports 19,998–19,999 achieved digits for
these returned midpoints, with explicit rounding warnings. This is distinct from the
20,000-digit source availability.

PR11dep0 now has an independently accepted 20,000-digit payload, with 20,008
stored digits and total absolute error bounded by
`3.3489120584697589184106699319561242E-20008`, including publication rounding.
Two complete bounded production profiles give 20,072 and 20,368 relative-radius
digits and the same unique 20,000-digit rounded value. Two isolated Python
standard-library verifier runs replayed both final interval convolutions and
authenticated the complete portable proof. Its large source objects remain
outside the runtime loop packs. Only the three direct epsilon-zero aliases gain
the new source capability; higher orders retain their individual limits.

The independent exact replay proves all three target identities and 21 of the
42 differential rows. The remaining 21 retain an explicit authenticated Kira
trust boundary. Full independent differential-equation validation is not claimed.
The practical Mathematica/AMFlow solver refinements agree through 125 observed
complex digits; their real central value agrees with the final period through
154 observed digits. This is a development crosscheck rather than a 20,000-digit
AMFlow computation. Mathematica remains unnecessary for the build, required
user tests, runtime and portable verifier.

== Pending delivery gates

- Finish the repository-wide selected CI checks; publish only after native
  and numerical compatibility pass.
- Publish the complete immutable external-data release before updating
  `data-latest`, then exercise actual downloads and offline reuse.
- Push the coordinated Vakint changes directly to `vakint_rustred` after the
  required final checks and upload.
