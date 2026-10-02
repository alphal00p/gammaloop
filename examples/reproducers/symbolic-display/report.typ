= Symbolic notebook display

This records the original bounded-preview experiment and its measurements.
The preview policy below was superseded on 2026-09-30: symbolic notebook
outputs now render the complete expression. Lazy backend caching and the
chain-collection optimization remain; the historical timings below are
retained as measured, not presented as full-render timings.

The chain-collection optimization records compatible matrix channels during the
existing regional observation pass. Metrics, vectors, and completed scalar dots
do not by themselves authorize matrix-chain matching. This changes scheduling,
not the permitted contractions or the resulting tensor.

Automatic tensor notebook displays inspect the encoded expression size before
presentation-port reconstruction, attachment discovery, render-tree encoding, or
Typst compilation. Ordinary formulas retain their full display. Large scalar
sums show a bounded prefix; other large expressions show an explicit summary.
The exact tensor, logical ports, metadata, and completion status are unchanged.
Explicit `to_html`, `to_svg`, `to_latex`, `to_typst`, and `format_tensor` still
export the complete expression.

`TensorExpression.formatted` retains the existing immutable tensor owner and
constructs no display backend. Symbolica's existing `FormattedOutput` caches
each requested backend independently. Concrete component tensors retain their
existing snapshot behavior; this change does not introduce a second display
wrapper or tensor carrier.

== Shared dependency

The accompanying `symbolica-lazy-formatting.patch` applies to released Symbolica
3.0.1, revision `e9a0d35490a834afb2fa40d97c8c0b1e58ef5a5e`. It extends the existing
`PythonFormattedOutput` owner, retaining its Python constructor and methods.
The tested local dependency commit is
`7d157e04f4ce4dce226ef87ca781aa4917d3735a` on `lazy-formatted-output`.
Errors remain retryable; unavailable optional representations are cached, and
recursive requests fail instead of deadlocking.

This owner is independent of tensor notation. Any Rust producer can supply a
`PythonFormattedOutput::lazy` callback for text, HTML, and LaTeX. Its existing
Python constructor still accepts ready strings. Standard `_repr_html_`,
`_repr_latex_`, and `_repr_pretty_` methods support Jupyter/IPython and Marimo;
each representation requested by a frontend is computed and cached separately.
Preview limits are the tensor producer's policy, not part of this shared owner.

For this development change the Cargo roots select a sibling checkout at
`/common/dev/symbolica-lazy-formatting`, using relative paths in the manifests.
To reproduce in another workspace, check out the same base revision at the
corresponding relative path and apply the accompanying patch with
`git apply --unidiff-zero`. All consuming
Cargo roots must use the same Symbolica source. The community host and bundled
notebook host are configured accordingly. No Cargo registry cache is modified.

This remains a local development dependency. `just ci-update` cannot resolve
the sibling checkout inside Nix's isolated source tree; the graph-generation
check stops at Cargo metadata, before running a CI suite. Native Cargo and
Hakari checks pass. Before CI publication, publish and pin the dependency
commit so the existing Git vendoring handles it. No custom Nix source-rewriting
path is introduced. CI remains disabled, and these changes have not been pushed.

== Checks

Native regressions count chainifier calls and candidate scans, verify surviving
regional channels, and exercise generic and mixed-representation chains,
compact gamma/vector products, trace preferences, filters, and scalar scopes.
The installed Python regression
`crates/spynso3/tests/installed_symbolic_preview.py` checks lazy backend calls,
caching, fallback, bounded large sums/products/powers, full explicit export,
logical ports, typed zeros, metadata, and exact expression preservation.

The formatter dependency includes tests for independent caches, clone sharing,
error retry, missing text, recursive rendering, and concurrent requests.
The existing tensor print and component-explorer tests remain applicable.

Final validation passes 815 Idenso tests, 125 Spynso tests, six formatter-owner
tests, and 37 installed Python regressions. Strict Clippy, formatting, generated
stubs, the API catalog, and documentation examples pass. IPython's actual MIME
formatter independently verifies that text and LaTeX requests do not compile
Typst, and repeated HTML requests compile once. Light/dark browser checks at
320 and 1280 pixels verify readable previews and unchanged ordinary MathML.
The Idenso run retains 25 existing skips; it is not an all-features exhaustive
test run.

The rebuilt extension is installed and the two owned notebook servers have
restarted. A separate fresh, read-only Marimo run completes the saved gluon
notebook in 11.70 s including startup, all cells, and display. It shows two
bounded rank-four previews and the expected 6,001-term result without JavaScript
or MathML errors. The active editor's unsaved changes and cached output snapshots
were preserved; rerunning its cells is required to update those snapshots.

== Measurements

The historical measurement used the input and settings below. To repeat it,
record the source revision, environment, samples and completion checks in a
separate output directory. The input is the
unprojected four-loop, three-gluon-rung ladder from the community notebook:
`Model.standard_model`, the V_36 vertex restriction, four loops, eight vertices,
`ladder_filter(outer_quark=False)`, numerator routing with `in_lmb=True`, expanded
couplings, and Lorentz dimension `three_gluon_rung_ladder::D`. It retains four
external ports. This is distinct from the faster projected scalar benchmark.

Release timings exclude imports, graph generation, admission, and compilation.
The timed call is `numerator.simplify_algebra(color=True, contract="fully")`.
The structural continuation times start from the same color-only result.
The full output is exactly equal to the saved baseline expression; its encoded
size remains 3,315,315 bytes. Disabling chain and trace collection still produces
that same result, but now takes approximately the same time as allowing them.

Preview times measure Python HTML generation, not browser painting. A large
tensor receives a summary; it is not fully rendered more quickly. Lazy
`formatted()` construction likewise postpones rendering: the first requested
HTML backend still pays the ordinary Typst cost for a small formula. Later
requests reuse that output. Explicit full exports and raw Symbolica expression
printers retain their existing behavior. Network formatting snapshots its
mutable source, and concrete component formatting retains its existing eager
snapshot path.

The primary three-sample comparison gives 2.081 s before and 0.843 s after for
the complete reduction. Structural contraction after color reduction falls
from 1.816 s to 0.595 s. A later warm run gives 0.781 s for the complete call;
the separate samples are retained rather than combined. Color-only reduction
remains about 8 ms.
The final rebuilt binary independently gives 0.818 s over three samples,
preserves the exact saved full tensor, and takes 0.156 ms for a completed rerun.

A CPU profile of eight complete reductions no longer samples the ineffective
chain collector. Instrumented regressions independently verify that it is not
called for finished metric/vector/dot coefficients. Remaining inclusive sample
costs include replacement observation maintenance (38.1%), structural contraction
(23.9%), and color-candidate discovery (21.5%). These overlap and must not be
added. This change leaves those paths and the algebra kernels intact.

Automatic HTML for the 3.3 MB result is a 1.3 kB summary generated in about
4 microseconds. For an ordinary small formula, lazy `formatted()` construction
falls from about 9 ms to 7 microseconds; the first HTML request still costs
about 10 ms and subsequent requests take about 5 microseconds.

A separate explicit-export comparison measures LaTeX source generation at
18 microseconds and Typst HTML generation at 9.08 ms for the same small formula.
Of the latter, 8.34 ms is inside `typst.compile`; preparation, temporary files,
and output processing account for the remainder. These compare serialization
with compilation, and exclude the browser's subsequent LaTeX typesetting or
MathML layout. A fresh Typst project is still used for each uncached formula.

== Complete-output follow-up: 2026-09-30

Automatic symbolic displays now render every term. The compact/explicit
invariant display, hollow unresolved-index placeholders, and color-symbol
migration reuse `tztusnmp`, `xmzwzlsl`, and `ywkyrvvr`. Lazy backend caching and
the chain-collection optimization remain in place.

The live gluon notebook constructs its projector from two metrics bound to the
generated ports and reduces the complete projected tensor in one call. It
displays the explicitly expanded tensor already needed by its scalar-coordinate
calculation. This exposes 6,001 top-level terms to the existing line wrapping;
the preceding factored result otherwise forms one enormous, unbreakable MathML
row. Neither representation is truncated.

Marimo independently rejected outputs above 8,000,000 bytes. The notebook host's
`tool.marimo.runtime.output_max_bytes` is now the platform-sized integer maximum.
The release extension was rebuilt and installed, the two notebook servers were
restarted, and the saved session was regenerated without either kind of preview.
The installed core has SHA256
`70434f75ef341916c07c42eb67787c147154feaee2134fbb8d3a5215b0178192`.

A fresh Chromium run matched all three complete formula fingerprints against
the executed session, including all 6,001 result terms, with no truncation,
MathML errors, JavaScript errors, or document overflow. It completed in 38.72 s,
including startup, calculation, and display. The large display cell took
11.81 s in the kernel. Projected algebra reduction alone took about 0.26 s,
explicit expansion about 12 ms, and a completed reduction rerun about 0.15 ms.
These are full-output measurements; the earlier 11.70 s notebook measurement
used bounded previews and is not an equivalent rendering workload.

The installed display/API suite passes all 45 tests; the generated-ladder suite
passes all five. Formatting, strict Spynso/Idenso Clippy, generated stubs, the
API catalog, and affected production-library checks pass. The broader native
selection passes 301 of 302 tests: its remaining component validator rejects
zero sampling calls after the tested identity has already produced 64 exact
zero components. An isolated diagnostic preserved that assertion and proved
the cancellation; the test was not relaxed. The broader all-targets caller
check also encounters seven type-inference ambiguities in unchanged GammaLoop
graph-parser tests. These failures prevent claiming a fully green workspace.

The full-formula audit also exposes two remaining notation issues: the small
raw-expression excerpt uses a different display style, and powers of dot
products retain their MathML grouping but lack visible parentheses. These are
separate from the removed output limits; this follow-up does not claim a
complete notation-polish audit.
