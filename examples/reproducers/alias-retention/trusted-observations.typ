#set page(paper: "a4", margin: 20mm)
#set text(size: 10pt)
#set heading(numbering: "1.")

= Passing terminal-trace facts to the tensor planner

The ordinary closed-trace kernel now publishes what it knows about its output:
only intrinsic scalar arithmetic and compact metrics remain in the visible
expression. The collector carries these facts through scalar sums and products;
the shared contraction implementation and explicit dot conversion preserve them.
This avoids reconstructing candidate inventories, ancestor leaf sets and regional
expression keys for the generated polynomial. The trace recurrence, retention
policy, public API and contraction implementation are unchanged.

This is a certificate for local syntax, not a declaration that every referenced
alias is finished. Dependencies retain their own domains and pending work. Exact
literal references are discovered lazily after normalization, so cancellations do
not retain nonexistent dependencies. Detailed leaf inventories remain absent,
rather than falsely empty: a later arbitrary collector still inspects the region
it selects. Unrelated identity heads prevent scalar certification and retain the
ordinary regional observations.

Callbacks and uncontrolled rewrites retain checked observation updates. Free-port
and axial outputs are not certified by the terminal evaluator. Wildcard slot
labels or dimensions are also excluded: construction can deliberately retain
inert diagonal metrics for those templates, even when a label occurs twice.
These exclusions restrict the bookkeeping proof. The existing trace evaluator
continues to accept the same index spellings.

== Matched before-and-after measurement

The input is the same production-generated quark ladder used by the preceding
#link("report.typ")[retention investigation]: symbolic Lorentz dimension $D$,
spinor trace unit four, the full scalar prefactor, and four bound graph labels
renamed to symbols. This last transformation keeps both modes on the same
terminal-trace route; the independent graph-label dispatch issue remains.
Gamma is enabled, colour and epsilon identities are disabled, contraction output
is `dots`, and the per-domain budget is 128.

Four release binaries were compared on the same host: the saved before versions,
and new alias-enabled and diagnostic never-alias versions. There were two rounds
with the second process order reversed, forty fresh reductions per alias-enabled
binary and twenty-four per never-alias binary. Compilation, model generation,
admission, warmup, profiler overhead and result destruction are outside the
reduction clocks. Local builds and test runs were paused for these measurements;
the host was shared and CPU affinity unrestricted.

#table(
  columns: (2fr, 1fr, 1fr, 1fr, 1fr),
  align: (left, right, right, right, right),
  inset: 5pt,
  table.header([Mode], [Before reduction], [After reduction], [Before + expansion], [After + expansion]),
  [Aliases], [91.24 ms], [31.31 ms], [192.72 ms], [128.42 ms],
  [Never aliases], [779.34 ms], [116.03 ms], [826.96 ms], [163.37 ms],
)

The last two columns add separately timed reduction and materialization means;
they request the same final expanded polynomial. Expansion itself averages
101.48 / 97.12 ms before / after with aliases and 47.63 / 47.34 ms without them.
Completed reruns after the change average 16.7 and 23.1 microseconds respectively.

The never-alias reduction is now 6.7 times faster than its previous implementation.
The between-mode reduction ratio falls from 8.5 to 3.7. Including expansion, the
new never-alias result costs 1.27 times the alias-enabled result. Thus the earlier
large discrepancy predominantly reflected processing of the generated expression,
not additional Clifford-algebra work.

The output representations are unchanged: 576 definitions and 168,894 stored
bytes with aliases; zero definitions and 1,387,316 bytes without aliases. Every
measured result reports `Complete` and materializes to the saved 525-term
polynomial exactly. The never-alias extension has a release-active trap at the
alias factory, in addition to the zero-alias assertion.

== Remaining work

Separate native profiles use `cpu-clock:u`, 65,528-byte DWARF captures and a
1,024-frame decode limit. Captures at requested 499 Hz and 99 Hz give category
shares within 0.52 percentage points; the table uses the 99 Hz repeat. FIFO
controls restrict capture to fresh reductions. Categories are exclusive,
assigned to the nearest matching frame from leaf to caller.

Both rates have zero lost samples, with exact agreement between record logs,
raw sample counts and the classifier. Kernel throttle notifications occur at
both rates, including the preceding baseline captures: requested cadence is not
fully achieved. The 99 Hz repeat contains 1,405 / 1,632 samples, approximately
87 / 90 samples per wall second, and reaches the Python binding in 99.7% / 95.8%
of weighted stacks. Rate sensitivity and these limitations are retained in the
JSON; phase estimates must not be treated as exact clocks.

#table(
  columns: (2fr, 1fr, 1fr),
  align: (left, right, right),
  inset: 5pt,
  table.header([Work after the change], [Aliases, estimated ms], [Never aliases, estimated ms]),
  [Observation and dependency bookkeeping], [2.76], [15.93],
  [Dot notation, excluding observations], [11.41], [87.10],
  [Trace kernel], [7.73], [6.47],
  [Other], [9.40], [6.54],
)

These are sampled CPU shares multiplied by independently measured mean wall
times, not instrumented phase clocks. In the never-alias profile, observations
fall from the preceding profile's 81.4% to 13.7%; dot notation now accounts for
75.1%. The kernel timings remain small and do not establish a meaningful
between-mode advantage. Exact lazy dependency inspection is still present;
expensive generated-result regional reconstruction has been removed from this
certified route. Dot normalization still traverses repeated unaliased payloads.
Other kernels and outputs without this certificate retain their existing
observation updates. In particular, the gluon ladder still has its large
never-alias cost: an additional correctness run retained the exact 6,001-term
polynomial and unchanged 118,522,487-byte unaliased representation. Its single
smoke reduction took 71.94 s, versus 68.79 s in the saved baseline. This run
overlapped test builds and is not a controlled timing comparison. The alias
representation also remained unchanged (616 definitions, 258,689 bytes).
This trace-focused change does not resolve the general contraction bottleneck.

== Correctness and reproduction

The updated extension also matches an independent integer Clifford-basis
component calculation for three momentum assignments, with exact traces
$-209664$, $125952$ and $46848$. This oracle performs explicit four-dimensional
index sums and blade multiplication, independently of Symbolica's trace
recurrence and tensor intake. The all-$D$ polynomial comparison is separate.

The regressions verify absent generated-result scans, preserved scalar
factorization, stable completed reruns, surviving epsilon candidates, later
arbitrary leaf collection, public metric coefficient collection, exact lazy
literal dependencies and wildcard exclusions. Existing callback/interface tests
remain part of validation.

All 853 enabled library tests passed: 852 under the regular nextest runner and
the unchanged exhaustive axial fixture separately in release mode (103.07 s).
The latter exceeds the regular runner's timeout on the baseline as well. The 22
normally ignored tests were not run. Strict Clippy (`-D warnings`), Rustfmt,
Ruff and both investigation-document builds passed. Cargo retains an unrelated
future-compatibility advisory for the dependency `proc-macro-error2` 2.0.1.

#link("trusted-observations-results.json")[`trusted-observations-results.json`]
records all individual timings, exact checks, commands, source and dependency
hashes, compiler settings, extension hashes, warmup phase timings and process
peak memory. The build provenance identifies jj change `srlrlonw`; the Rust source
manifest distinguishes its implementation from later reporting-only changes.

Use #link("profile.py")[`profile.py`] with `--materialize` to record expansion and
completed-rerun timings outside the profiling window. Supply `--python-root` to
select an isolated extension, `--compare` for the saved reference, and the same
settings/input options recorded in the JSON. The updated diagnostic patch is
#link("trusted-observations-no-aliases.patch")[`trusted-observations-no-aliases.patch`].
It is used only for the experiment; the working source was restored byte-for-byte.
The measured extensions were isolated builds; the notebook installation was
unchanged. This focused experiment includes no fresh FORM timing. It preserves
the existing trace algorithm and uses independent components and exact polynomial
comparisons to check the changed bookkeeping.
