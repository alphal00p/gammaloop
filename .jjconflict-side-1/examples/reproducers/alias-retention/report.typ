#set page(paper: "a4", margin: 20mm)
#set text(size: 10pt)
#set heading(numbering: "1.")

= Never creating aliases: ladder experiment

Diagnostic measurement, 30 September 2026. Automatic retention was disabled
at its creation sites, keeping sums factorized. A release-mode assertion in
`SymbolicTensor::alias_handle` rejects any missed path. It was independently
triggered to verify that the trap is active. Neither successful ladder run
called it; both returned zero definitions and reported `Complete`.

== Results

Wall seconds, measured in fresh processes. “Total” means reduction plus explicit
expansion into the final polynomial; it excludes setup, generation, admission,
and output serialization. All phase measurements, including setup and CPU time,
are in #link("results.json")[`results.json`]. These are diagnostic observations,
not a statistical benchmark or a claim about all tensor expressions.

#table(
  columns: (auto, auto, auto, auto, auto, auto),
  align: (left, left, right, right, right, right),
  inset: 5pt,
  table.header([Ladder], [Retention], [Reduction], [Expansion], [Total], [Aliases]),
  [Gluon], [On], [0.216], [8.269], [8.485], [616],
  [Gluon], [Off], [68.792], [2.632], [71.424], [0],
  [Quark], [On], [32.046], [0.639], [32.685], [8,831],
  [Quark], [Off], [14.382], [0.042], [14.424], [0],
)

The gluon result grows from 258,689 stored bytes to 118,522,487 bytes when
identical subexpressions must be copied into their parents. Process peak RSS
rises from about 600 MiB to 2,664 MiB. The eight-level dependency graph therefore
provides substantial sharing despite its modest depth.

The quark result shrinks from 2,535,604 stored bytes to 1,759,507 bytes, but process
peak RSS rises from about 283 MiB to 599 MiB. Removing retention accelerates this
particular reduction while increasing transient memory. RSS includes generation,
validation and the exact comparison, and is not a kernel-only memory measurement.

Completed reruns take approximately 11–34 microseconds in these runs. The
whole measured pipeline, including generation and input preparation, takes
8.581 versus 71.522 seconds for gluons and 32.836 versus 14.609 seconds for quarks.
Compilation and Python startup are excluded from these phase totals.

== Inputs and correctness

Both cases come from the production generated-ladder helpers: four loops,
eight vertices, eleven internal edges and three gluon rungs. The quark ring is
massless. Lorentz dimension is symbolic $D$; the spinor trace unit is four.
The external metric/colour projection is applied before splitting factors.
Only the spacetime factor is reduced here; colour algebra and the quark result's
subsequent loop-momentum routing are outside the measurement. The gluon input
already uses loop-momentum routing, matching its production notebook.

Settings are `gamma=False` for gluons and `gamma=True` for quarks, with
`color=False`, `epsilon=False` and `contract="dots"`. The gluon comparison uses
16 steps per domain. With no aliases the quark run exhausts that budget, because
its work now stays in one domain instead of thousands of definition domains.
Both quark modes were therefore measured with the same budget of 128. The
incomplete 7.9-second run is not counted as a speedup.

The final 6,001-term gluon polynomial and 525-term quark polynomial compare
exactly equal to their saved baseline expressions. Sums are not globally
expanded before reduction. Expansion is an explicit, separately timed request.
The no-alias budget-16 quark result can also be resumed to the same exact result.

== The fast quark path and the historical benchmark

The quark slowdown is triggered by a scalar-context eligibility check in
`crates/idenso/src/dirac/simplify.rs`: it requires explicit index labels to be
variables. Generated graph labels such as `hedge(23,1)` are functions. A
standalone trace accepts these labels, but adding its scalar prefactor selects
the stricter check and loses the fast terminal trace evaluation.

A controlled diagnostic retains the full scalar prefactor, collects the trace,
and changes only its four bound dummy labels to fresh symbol names. Both builds
then take the fast path and produce exactly the same 525-term polynomial:

#table(
  columns: (auto, auto, auto, auto, auto, auto),
  align: (left, right, right, right, right, right),
  inset: 5pt,
  table.header([Retention], [Reduction], [Expansion], [Total], [Aliases], [Stored MB]),
  [On], [0.096], [0.106], [0.202], [576], [0.169],
  [Off], [0.904], [0.053], [0.958], [0], [1.387],
)

These are reduction-plus-expansion seconds, excluding the approximately
2.4 ms trace-collection step. Use `--symbol-indices` with the quark driver to
reproduce this diagnostic. Renaming is an experiment, not a notebook workaround
or a proposed requirement on users' dummy labels. The ordinary generated input
should receive the same eligible optimization after its invariants are checked.

Consequently, the earlier 0.095-second standalone result already used aliases.
The slow-path comparison alone does not determine their value on the optimized
path. For the matched fast-path input, retention saves copying and processing
repeated factorized subexpressions, reducing the stored payload by about eight
and the reduction-plus-expansion time by about five in this observation.
The native profile below now distinguishes trace construction from subsequent
observation and notation work; the original end-to-end timings alone did not.

== Native profile of the matched quark reduction

Fresh reductions of the same admitted input were measured in sequential
processes, with retention order reversed in the second pair. Forty alias runs
average 98.3 ms; 24 no-alias runs average 801.8 ms, a factor of 8.2. These repeat
measurements exclude profiling, generation, admission, expansion, completed
reruns and result destruction. The original 95.8/904.5 ms figures above were
single fresh-process diagnostic measurements.

Separate native CPU profiles use the unchanged binaries, `cpu-clock:u` at a
requested 499 Hz, 65,528-byte DWARF stack captures, and a 1,024-frame unwind
limit. FIFO controls enable recording only around the reduction loop. They
contain 3,057 alias samples and 6,519 no-alias samples, with no lost samples.
Each sample is assigned to the nearest relevant frame walking from the leaf
toward the caller; categories below are disjoint. Observation work nested inside
dot conversion belongs to observations, not both categories.

#table(
  columns: (2fr, 1fr, 1fr, 1fr),
  align: (left, right, right, right),
  inset: 5pt,
  table.header([Work], [No-alias share], [Aliases, est. ms], [No aliases, est. ms]),
  [Observation maintenance], [81.4%], [57.5], [652.9],
  [Dot conversion, excluding observations], [14.7%], [16.2], [117.7],
  [Trace kernel, including its expression construction], [1.0%], [8.6], [7.7],
  [Other work], [2.9%], [16.0], [23.5],
  [Total], [100%], [98.3], [801.8],
)

The milliseconds are estimates: independently measured mean reduction time
multiplied by sampled CPU share. They are not exact instrumented phase clocks.
The small trace-kernel difference is not evidence that either representation is
faster there. Observation maintenance explains approximately 85% of the additional
elapsed time; dot conversion explains approximately 14%. Actual gamma evaluation
does not explain the factor of eight.

The expensive owner is `tensor/simplification/observation.rs`:

- `scan_regions` copies every visited expression node into owned storage,
  including overlapping sum/product subtrees. Functions are also copied into
  leaf sets.
- `finish_observation` merges descendant sets into each ancestor and retains
  observations keyed by the complete arithmetic subtree. Hashing these keys
  reads their full encoded bytes.
- `updated` clones the old scope map and only reuses matching top-level regions.
  A changed region is scanned recursively even when some nested subexpressions
  survive unchanged.
- `RegionObservation::observe` repeatedly looks up the fixed `scalar_alias`
  symbol. Native stacks place symbol-building and built-in lookup work under
  this observer.

Thus the large, still-factorized no-alias result incurs repeated copying,
hashing, allocation and merging while constructing planning observations.
The aliases mostly help this implementation by keeping those observed bodies
small. This corrects the earlier emphasis on copying inside the trace recurrence:
the profile identifies processing of the generated result as the dominant cost.
It does not show that gamma algebra intrinsically needs aliasing for speed.

No new phase instrumentation, source patch, extension rebuild or notebook reload
was needed. Both variants return `Complete`, retain the expected 576/0 aliases,
and their final expanded expressions match the saved 525-term polynomial exactly.
Each process checks its warmup and final timed result outside the profile window.

#link("profile-results.json")[`profile-results.json`] records timings, sample
counts, exclusive/inclusive attribution, stack coverage, binary hashes and raw
artifact paths. The Python binding appears in 99.5% of alias stacks and 96.7%
of no-alias stacks. A preliminary 16 KiB capture truncated many outer frames,
but its observation/dot/trace shares agree within 0.2 percentage points with
the deeper capture. Inlined source-level attribution is unavailable in these
release binaries.

A later health audit found kernel throttle notifications in these 499 Hz
captures, despite zero lost samples. Their requested cadence was not fully
achieved. See #link("trusted-observations.typ")[the producer-observation followup]
for a lower-rate sensitivity check, the implemented change and matched new
measurements; the phase numbers remain estimates.

Use #link("profile.py")[`profile.py`] with `--mode`, `--iterations`,
`--hep-examples`, `--compare`, and `--save`. Select the isolated no-alias package
with `--python-root`. Omit `--control` and `--ack` for ordinary wall timings.
For profiling, create two FIFOs and supply their paths to both the script and
`perf record --delay=-1 --control=fifo:CONTROL,ACK`. Decode using:

```sh
perf script --max-stack 1024 -i capture.data \
  -F comm,pid,tid,time,event,period,ip,sym,dso > capture.stacks
python profile_stacks.py capture.stacks
```

The #link("profile_stacks.py")[sample classifier] retains inclusive counts
separately; they must not be summed to estimate total runtime.

== A concrete shared subtrace

Both retention modes already memoize the same subword recurrence. Aliases do
not explain a reduction in the number of recurrence evaluations in this matched
comparison. They preserve sharing when those cached answers are assembled into
larger expressions. Ordinary Symbolica sums and products embed their children's
byte representations; retrieving a cached polynomial does not make later
occurrences share its storage.

The audit extracts the actual alias definitions from the fast quark result.
Writing $d_(i j) = q_i dot q_j$, definition 53 is

$ A_53 = 4 (d_(2 3) d_(9 10) - d_(2 9) d_(3 10) + d_(2 10) d_(3 9)). $

This is the four-gamma trace for momenta $(q_2, q_3, q_9, q_10)$. It has 27
distinct immediate parent definitions. Three actual parents, numbered 230,
231 and 232, have the common form

$ B(p) = d_(p 9) A_43 - d_(p 10) A_46 + d_(p 2) A_49
  - d_(p 3) A_52 + d_(p 8) A_53, quad p in {5, 4, 7}. $

The other four handles denote the corresponding four-gamma traces. Their full
definitions and the parent list are recorded in
#link("walkthrough.json")[`walkthrough.json`]. Thus these parents store short
references to the same sum rather than embedding it again. The body of $A_53$
occupies 278 bytes and its handle 9 bytes. Counting only one occurrence per
immediate parent gives $27 times 278 = 7506$ bytes of embedded bodies versus
$278 + 27 times 9 = 521$ bytes with sharing. This illustrates payload copying;
it is not an allocator/RSS model and excludes headers and registry bookkeeping.
Other shared bodies compose similarly, explaining why shallow nesting can still
save substantial work.

The audit also checks the full 16-gamma word independently at $D=4$. It sums
all $4^4=256$ assignments to the four contracted Lorentz indices, multiplying
vectors in the 16-component Clifford basis with exact integer coefficients and
metric $(+, -, -, -)$. It uses neither the Symbolica trace recurrence nor its
tensor contraction/admission implementation. Three deterministic momentum
assignments give $-209664$, $125952$, and $46848$, each exactly matching the
Symbolica result after removing the common scalar prefactor. These are independent
component checks; they supplement, rather than replace, exact symbolic equality
of the two retention modes for symbolic $D$.

Run #link("walkthrough.py")[`walkthrough.py`] with the same arguments as the
quark driver, adding `--symbol-indices` and using the no-alias binary result
as `--compare`. It checks exact equality, extracts the sharing example and
writes the independent checks to `--save` with suffix `.audit.json`.

== Why there are two execution routes

The terminal evaluator handles an eligible whole gamma word in one invocation.
It contracts repeated indices inside that word and caches subwords before
emitting scalar expressions. General rewriting also handles open chains,
projectors, special matrices and callback-sensitive contexts. It applies local
identities, publishes intermediate traces/metrics and lets subsequent scheduled
work finish them. It can eventually call the same trace kernel, but after a
different sequence of intermediate construction and contraction.

The existence of general identities does not justify excluding an already
admitted abstract index merely because its label is a function. In this case the
`AtomView::Var` gate couples the algorithm to index spelling and the presence of
a scalar spectator. Eligibility should use established slot identity, dimension,
scope and callback facts. This investigation does not change that production
gate. It also does not justify a separate public aliased tensor type: the measured
benefit is internal sharing, which can be represented behind one tensor API.

Change `uqrmxmwx` (`90f1a6e752161c243f0fe51d3294b21602fcce08`) also used
aliases, as did its baseline. Its recorded improvements came from retaining
validation/interface facts, inspecting one representative term of homogeneous
trace sums, carrying intrinsic-normalization proofs and avoiding unnecessary
collection of standalone traces. Its warm native Rust timings explicitly exclude
final alias resolution. See
#link("../trusted-gamma-algebra/report.typ")[the historical report] and
#link("../trusted-gamma-algebra/driver.rs")[its driver]. Those results are not an
alias-on versus alias-off comparison and are not fully materialized Python
ladder timings.

== Reproduction and scope

#link("environment.json")[`environment.json`] records revisions, the release
build command, dependency override, CPU affinity and extension hashes.
#link("no-generated-aliases.patch")[`no-generated-aliases.patch`] is a diagnostic
patch, not a supported public API. Apply it only to an isolated build: it
intentionally rejects any operation that still needs an alias handle, including
explicit user alias construction. It does not establish support for every
callback, dummy-scope or tensor-collection case.

Use #link("ladder.py")[`ladder.py`] with the unchanged extension first:

```sh
python ladder.py gluon --hep-examples /path/to/community/examples/hep \
  --save /tmp/gluon-baseline.bin
python ladder.py quark --hep-examples /path/to/community/examples/hep \
  --save /tmp/quark-baseline.bin --max-steps-per-domain 128
```

Build the diagnostic extension with the same release configuration and use its
Python package in a separate environment. Repeat the commands with distinct
`--save` paths, `--no-aliases`, and `--compare` pointing to the corresponding
baseline file. Keep the quark budget at 128 in both builds. The driver requires
complete results and asserts exact equality. It writes timing and size records
alongside each binary expression.

The experimental source edits have been reverted; the notebook extension was
not replaced. The separate scalar-prefactor scheduling regression is unchanged
in both builds. Blanket removal of aliases is not supported by these results:
sharing helps both the gluon contraction and the fast quark path; the quark
scalar-context eligibility check needs correction. Fix the scalar-context dispatch discrepancy before
choosing a new default or exposing a permanent retention option.
