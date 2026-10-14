= Explicit native cache preparation

This interface is opt-in. Constructing Vakint and evaluating without an attached
cache keep the existing behavior; no cache files are read, written or downloaded
automatically. Cached payloads are exact, unit-mass, unnormalized reductions and
their complete native nonzero conditions. They are not numerical master values,
epsilon expansions, or a claim of additional reduction coverage.

== Current integration status

The current dependency is RustRed
`78cbc5c1c6e98d057a20390534cae25184e7ab99`. Preparation now captures newly
computed admitted intermediate reductions as well as selected roots. Each
intermediate carries its own exact nonzero conditions, derived from bounded
request-local dependency provenance; it is not assigned the entire parent
guard list and is never exported from an unguarded coefficient memo.
Ordinary evaluation remains read-only against attached packs.

Fresh native one-, two- and three-loop default primers completed all 5, 77 and
2,320 selected roots with the existing eight-chunk lookup allowance. Preparation
runner times were 0.00457, 0.00709 and 1.55330 seconds; including native program
setup, the corresponding measurements were 0.04844, 0.04599 and 2.28218 seconds.
These are shared-host functional timings, not matched speedup claims. The first
two primers already contain every captured descendant. An exact metadata census
of native verified committed snapshots confirms 5, 77 and 2,352 distinct keys,
respectively: the three-loop pack adds 32 keys outside its 2,320-root primer.
No coefficient import or custom binary parser is used for this census.

The enriched four-parent 38-root gate completed in 94.203 seconds of preparation,
101.495 seconds of native-process wall time, with a 3.472 GB sampled process-tree
peak. It stored 134,821 record references in 34,303,404 charged bytes. The exact
committed-key census confirms 134,821 distinct family-qualified keys, including
134,783 outside the selected 38-root workload. Six selected
roots reused captured intermediate data. Two requests exhausted their optional
provenance-metadata allowance and continued normally; this is not a claim that
every internal node was captured. All 38 selected records and eight non-selected
records per parent passed fresh native exact coefficient and guard-set comparisons.
The subsequent full enriched four-loop primer also completed all 164,352
selected roots. Its runner took 552.794 seconds (9 min 13 s), with 564.834
seconds including native setup and a sampled process-tree peak of 5.987 GB.
The pack uses 90,790,870 charged bytes. An exact committed-key census confirms
286,527 distinct family-qualified keys: the selected 164,352 plus 122,175
additional intermediate keys. One optional provenance-metadata exhaustion was
reported, without a selected-root or checkpoint failure. Per-parent run counts
are 1,887 / 778 / 922 / 2,493 for H / FG / BMW / X, within the unchanged 4,096
limit. This is the requested finite warmup, not additional integral coverage.

Fresh native verification of the full pack passed 38 selected records and 32
non-primer samples, including exact coefficient and guard-set agreement. A
same-engine comparison then loaded the old full root-only pack and the enriched
full pack in separate fresh processes for each held-out input:

#table(columns: 4,
  [Parent], [Root-only full pack], [Enriched full pack], [Native rules, old → new],
  [H], [3.29911 s], [3.04448 s], [18,538 → 1,850],
  [FG], [0.22993 s], [0.05056 s], [1,313 → 15],
  [BMW], [1.36707 s], [0.70897 s], [9,536 → 251],
  [X], [17.76424 s], [12.39088 s], [68,124 → 6,297],
)

These are first applications after loading, excluding package load. Every
input remains absent as an entry and has positive intermediate backing hits;
exact ordinary-reference and RAM-repeat parity pass. These are single matched
samples on a shared host, with concurrent exact-record checks on another core
and independent review of bindings and correctness receipts, not statistical
speedup guarantees. The older 3 h 22 min preparation was not
rerun on this engine, so its ratio to 9 min 13 s is not a controlled speedup.

Four fresh-process public held-out calls used the enriched 38-root pack, followed
by a same-process RAM repeat and an ordinary reference. All outputs agreed exactly.
Each first call had zero entry hits and positive descendant hits: these genuinely
unseen roots used intermediate shortcuts.

#table(columns: 4,
  [Parent], [Loaded first call], [RAM repeat], [Ordinary reference],
  [H], [3.88129 s], [0.05320 s], [3.78110 s],
  [FG], [0.28171 s], [0.02093 s], [0.30175 s],
  [BMW], [2.32338 s], [0.02020 s], [2.22645 s],
  [X], [25.25758 s], [0.09448 s], [21.00561 s],
)

The ordinary reference follows cached execution in the same process, so shared
initialization is already warm; these are not paired cold-process speedups.
Application uses unchanged defaults: 1 GiB aggregate retained coefficients,
256 MiB per four-loop namespace. Preparation instead uses decimal 16 GB
aggregate, 4 GB per namespace. X incurred 12,521 application evictions and
2.489 GB of decoded chunk traffic. Fewer native rules therefore do not imply a
faster first call. Enlarging the preparation batch from 64 to 256 reduced this
gate's preparation to 84.848 seconds and its pack to 26,089,874 bytes, but made
every first held-out call slower (H 5.59869, FG 0.40846, BMW 3.41372, X 32.49444
seconds). The default remains 64; no application allowance was raised.

New packs write index schema 03 with compact per-chunk key fingerprints. A
definite negative avoids reading an unrelated chunk; a possible match still
undergoes full native exact-key and payload validation. This fixes a measured
three-loop refusal caused by ten overlapping runs against the unchanged
eight-chunk lookup grant. Native record and bundle schemas stay 02, and older
schema-02 indexes remain readable. Old root-only preparation checkpoints are
not resumable as intermediate-complete work: regenerate into a fresh directory.
Their generated caches remain usable for ordinary evaluation.

The focused Rust cache/application and native public gates are complete at this
revision. A freshly rebuilt private coherent Symbolica 3.0.1 host also passes
all eleven Python cache smoke tests and 32 additional public API checks,
including lazy loading and exact evaluation of the complete enriched four-loop
pack. The host uses the explicitly pinned historical Community baseline; this
is not a Community-main/PyPI release or a full-workspace CI claim. Earlier Python
receipts below retain their older pinned dependency. No numerical masters or
IBP rules are generated by these preparation or validation steps.

== Earlier root-only four-loop preparation and reader gates

As of October 14, the complete default four-loop primer has prepared and durably
stored all 164,352 requested roots. Its native completion report records
12,171.30 seconds (3 h 22 min 51 s) and 72,253,397 physical/charged bytes
(72.253 MB). This is the earlier selected-root-only pack, retained as a baseline.
No new IBP generation was needed. These are native
completion measurements; the missing final supervisor receipt does not support
a full-run peak-RSS claim.

RustRed `0974effab7bde243e5360c96ae7296639d5018e5` adds candidate descendant
reuse: an unseen starting input can reuse a matching completed guarded record
for an intermediate integral. The prepared schema-02 pack and its exact
artifact/terminal-convention binding remain unchanged. The first 38 records
(431 coefficients and 60,007 guards) agree with old-pack and fresh native
controls. Four unseen power-four inputs also agree exactly in coefficients and
guards; the updated X check passes with LRU at the unchanged 256 MiB namespace
allowance.

Measured elapsed-time effects are mixed: fewer rule applications improved the
sampled FG/BMW calls, while H/X became slower because repeated chunk reads and
decompression outweighed the saved reduction work. That reader milestone,
RustRed `301beec66e6a2fcdc12ffe2c665181206615550c`, retains bounded exact chunk-key
metadata to avoid repeated known-absent reads. Native matched checks preserve
all terms and guards while reducing H chunk reads from 12,924 to 3,293 and X
from 47,593 to 9,721. These I/O counts establish the mechanism; no general
elapsed-time speedup claim follows from them, descendant hits or exact repeats.

The same Rust public interface was also measured on those four held-out
power-four requests, with LRU, unchanged allowances and the same prepared pack.
Setup is excluded. The comparison isolates the reader change; both columns
already use candidate descendant reuse.

#table(columns: 3,
  [Parent], [Range-only reader], [Bounded exact-key reader],
  [H], [4.72166 s], [3.34233 s],
  [FG], [0.28129 s], [0.22426 s],
  [BMW], [2.02304 s], [1.40644 s],
  [X], [30.63170 s], [18.15692 s],
)

These are single fresh-process samples on a shared host, not statistical
guarantees. Exact outputs, native rule counts and descendant hits are unchanged;
the eliminated reads match `key_index_skips`, with 0.69--5.78 MB of additional
metadata in these public probes.

The primer is a finite operational selection: all physical parent pinches with
at most two extra positive powers and numerator rank at most one. Completing
it does not establish unrestricted integral coverage or candidate closure.
The earlier lower-loop, collector and Python receipts below retain their
original scope. Explicit import/download APIs and their completed native/private
Python-runtime gates are described below. Five-loop cache preparation is deferred
until the updated artifact is available.

== Explicit import and download

On Linux and macOS, `import_cache` and `download_cache` publish an independently
validated native snapshot into a new directory. Both require the SHA-256 of a
trusted USTAR gzip archive, obtained independently of the archive bytes. Its
single root is `cache/`, containing the native `bundle.rrindex` and `bindings/`
layout. Only regular files and directories are accepted; links, sparse files,
PAX/GNU extension headers, duplicate paths and traversal paths are refused.

```python
report = engine.import_cache(
    4, "vakint-native-cache-4l-ustar.tar.gz", "local-four-loop-cache",
    sha256=trusted_archive_sha256,
)
engine.load_cache(4, report["output"])
```

Alternatively, download from a supplied final HTTP(S) URL; redirects are refused:

```python
report = engine.download_cache(
    4, trusted_final_url, "downloaded-four-loop-cache",
    sha256=trusted_archive_sha256, timeout_seconds=3600,
)
engine.load_cache(4, "downloaded-four-loop-cache")
```

The local October 14 four-loop transport archive is 63,548,634 bytes and has
SHA-256 `de037d50a105c31d6126563525dc951436e252037b47ad6f763951b1f66a6276`.
The one-through-four-loop archives use the names
`vakint-native-cache-1l-ustar.tar.gz` through
`vakint-native-cache-4l-ustar.tar.gz`. These are local delivery artifacts; no
public hosting URL is implied. Use the corresponding independently supplied
checksum for each file, and load the importer’s output directory directly,
without appending another `cache/` component.

Those named archives are the earlier root-only delivery. The new enriched
delivery uses `vakint-native-cache-1l-intermediates-index03.tar.gz` through
`vakint-native-cache-4l-intermediates-index03.tar.gz`. The four-loop archive is
89,680,390 bytes, SHA-256
`f0b718fe01e6801eb9d3500959fc0087d46a50ec64fd65eced90af18e4942ea2`.
It contains the complete finite preparation with 286,527 distinct keys, including
122,175 outside the selected primer. These are local, untracked data assets,
not a hosted release. The separately delivered manifest supplies every archive's
checksum, source-pack identity and exact-key census; do not infer a download URL.

To package an existing native bundle whose directory is named `cache`, choose
a new output archive path and use standard USTAR plus deterministic gzip:

```sh
tar --format=ustar --sort=name --owner=0 --group=0 --numeric-owner \
    --mtime=@0 --mode='u=rwX,go=rX' --use-compress-program='gzip -n -6' \
    --create --file=vakint-native-cache-4l-ustar.tar.gz \
    --directory=/path/to/prepared-parent cache
sha256sum vakint-native-cache-4l-ustar.tar.gz
```

This only packages existing bytes. Distribute the archive through the chosen
hosting service and its checksum through a trusted release channel; `download_cache`
does not discover a host or obtain its own trust anchor.

Transport does not attach a cache, evaluate an integral, regenerate rules or
compute numerical masters. Normal evaluation never downloads. A byte-verified
archive is not an independent mathematical certificate: native generated-cache
provenance remains the trust boundary, and exact shipped program/terminal
bindings must match before publication. Transport accepts the native index
versions understood by the pinned RustRed dependency; it does not change record
payloads. The named historical archives above contain schema-02 indexes.

The Rust methods take `RustRedCacheTransferOptions::new([u8; 32])`. Defaults allow
16 billion archive bytes, 16 billion extracted file bytes, 100,000 entries and
48 billion concurrent temporary bytes. Python exposes `max_archive_bytes`,
`max_extracted_bytes`, `max_entries`, `max_temporary_bytes` and `disk_bytes`
overrides. Bounds count logical bytes, not filesystem allocated blocks or
metadata, process RSS or reduction RAM. Tar headers and padding have a separate
derived decompression cap. HTTP requests ask for identity encoding; the archive
checksum and byte count apply after any HTTP content decoding by the transport.

The implementation streams checksum verification before extraction. RustRed then
checks exact namespace bindings and copies every referenced native chunk with
its existing bounded, checksummed 64 KiB copy path. The staging allowance covers
the archive, extracted payload and verified copy simultaneously. Only the
verified native snapshot is published with an atomic no-replace rename; existing
destinations always refuse, including empty directories. Temporary files are
cleaned on failure. A directory-sync failure after the rename is reported as an
error and may leave the complete published destination; it is never removed as
rollback. Preparation checkpoints and unrelated archive files do not enter the
published snapshot, which is loaded separately using `load_cache`.

== Earlier four-loop collector gate

The initial local four-loop gate was incomplete: it checkpointed seven roots before
the native retained-term census requested 16,000,086 terms against the
unchanged 16,000,000-term ceiling. The next round-robin request was the X
parent key `(1,1,1,1,1,1,1,1,2,0)`: nine distinct positive denominator slots,
positive-power sum 10 and numerator rank zero. This is an operational refusal,
not an out-of-scope or mathematical reduction result. Fresh-process public
four-loop hit, warm, held-out and numerical-table comparisons were not reached.
The partial pack and failing test receipt are retained; no assertion or
allowance was relaxed.

That functional local-patch run used frozen executable
`126871df24f6bf2a39b84ff0713d9ec490d5932b834e51f47c30f89e3e8e482b`.
Native program loading took 6.681 seconds and preparation 19.757 seconds;
the supervised process-tree peak was 3.103 GB, without a memory-guard stop.
Preparation had a decimal 16-billion-byte aggregate coefficient allowance,
partitioned equally over four sessions. Its separate native term ceiling was
not increased. The planned fresh reader uses the ordinary 1 GiB aggregate
allowance (256 MiB per namespace), not the preparation allowance.

== Earlier passing four-loop first-38 gate

The matched diagnostic established that ordinary coefficient reduction at the
failed X point used about 1.95 million retained terms, while duplicating complete
guard vectors in every guarded subtree exhausted 16 million terms. The repaired
native path collects one bounded request guard union only for provably fresh
candidate sessions. It preserves the legacy guarded API, original poles and
shared memory accounting; no configured allowance was raised.

The unchanged 38-root test then passed with frozen executable
`0789c7e6ec22bd82c9d2055672550948919301827306ec07b23a243e09673473`.
Preparation completed all 38 roots in 173.084 seconds and stored 397,885 bytes.
Initial program loading took 7.507 seconds; the fresh child separately loaded
programs in 7.145 seconds and attached the indexes in 0.001 seconds, with no
coefficient reads during attachment. The full supervised run took 278.452
seconds with a 4.468 GB sampled process-tree peak and no guard stop.

The gate checked 76 mass-squared-1/4 prepared public evaluations and RAM
repeats, four independent exact coefficient/ordered-guard records, eight
denominator-cancellation/pinch comparisons, four held-out cubed-line roots,
four substitutions using existing numerical master tables, and clear/reload.
The prepared phase recorded 38 backing hits, 114 RAM hits and zero IBP rule
applications. Ordinary shared reducers were unchanged by preparation and cached
calls. Later cancellation and held-out requests performed real native work;
they are not counted as prepared hits. That earlier gate did not exercise
candidate descendant-cache reuse.

#table(columns: 3,
  [Held-out family], [First public call], [RAM repeat],
  [H], [3.128 s], [0.048 s],
  [FG], [0.271 s], [0.018 s],
  [BMW], [1.623 s], [0.016 s],
  [X], [15.897 s], [0.083 s],
)

These are same-input functional measurements on a shared host, not a general
performance claim. Preparation uses decimal 16 GB aggregate RAM; fresh attachment
uses the unchanged native 1 GiB aggregate default. Passing this finite sample
did not by itself certify candidate closure or completion of the
164,352-request primer. The later complete preparation is recorded above.

That RustRed dependency used mandatory gzip chunks and cache-pack/bundle
schema 02; current index evolution is described above. The earlier `59e45868` gate below used schema 01: its development
packs must be prepared again. Native artifact and record schemas are unchanged.
No stable on-disk migration contract is promised here, and an incompatible pack
is refused rather than guessed at.

== Initial one-loop gate (schema 01)

The optimized actual-Git test prepared native roots 0 through 4, then spawned a
fresh process. Public symbolic I4 evaluation at mass squared 4 matched the
uncached reference. The loaded engine recorded one backing hit and zero rule
applications on that call; a repeated call hit RAM. Held-out I5 ran native
reduction and matched its reference. Clearing RAM retained the backing, and an
explicit one-byte RAM policy still refused the result after default reload.
No numerical master values or new rules were computed.

The libtest result was one pass in 0.12 seconds. GNU time reported 0.11 seconds
wall time (centisecond resolution), 0.10 seconds user and 0.01 seconds system,
with maximum resident size 30,924 KiB. The supervisor's 2.175 seconds includes
launch, sampling and process-group drain; its sampled tree peak 18.96 MB can
miss a short-lived peak. These are correctness-gate measurements, not a paired
speed comparison. The reference and held-out calls perform native work: the
zero-rule count applies only to the loaded engine's I4 hit. Per-stage timings,
actual disk bytes and descendant-cache reuse are not established by this gate.

The frozen executable SHA256 is
`531a731b069b931a4d630be32afe19d3e0e5f0360d97c8078e35012a51f165c7`.
The assertion-only test source is
`065a40827c4abb470a36c3c08bbd01b0c4caed29c569ffcaf9629534656526ce`;
later uncompiled timing-print edits were removed from this checkpoint.

== Earlier one-through-three-loop gates (schema 02)

One frozen optimized executable, SHA256
`ecfc2551dac3895d4c8ccb18f21023641d150f2db47474b99dc567ff9d16304e`,
tested the actual Git dependency `a0649fed8f62f7a9fb7c051aa8bd7fc6938d6399`.
The one-loop fixture again prepares powers 0 through 4. The two-/three-loop
fixtures use all physical pinches with one extra dot and numerator rank zero.
Their complete finite workloads contain 20 and 256 roots, respectively.

#table(columns: 5,
  [Loops], [Roots], [Preparation], [Fresh load], [Disk bytes],
  [1], [5], [0.008873 s], [0.000522 s], [3,251],
  [2], [20], [0.010092 s], [0.003162 s], [3,927],
  [3], [256], [1.101400 s], [0.776528 s], [12,167],
)

Preparation timings exclude `Vakint::new` but include program loading and
workload planning. Fresh load includes native artifact loading in the child
process; it performs zero coefficient-chunk reads. These are functional
shared-host measurements, not controlled speedup comparisons.

The two-/three-loop children compare prepared dotted full-support roots and
known pinches at mass squared 1 and 4 with independent uncached evaluations;
the three-loop pinch itself is undotted. Prepared evaluations and their warm
repeats use zero native rule applications. A separate typed-record comparison
checks exact unit-mass coefficients and the full native guard set, ignoring
guard order only. Held-out cubed-line roots still perform native reductions:
the two-loop call uses 5 rules and one descendant backing hit; the three-loop
call uses 62 rules, one descendant RAM hit and 32 descendant backing hits.
The one-loop held-out I5 call now verifies reuse of resident I4 with one rule.
All values agree exactly; no numerical master data or new rules are produced.

The three parent tests pass in 0.15, 0.13 and 2.21 seconds. GNU time wall values
are 0.16, 0.13 and 2.23 seconds; their per-process peak RSS values are 12,324,
15,432 and 245,876 KiB. The three-loop supervisor sampled 516.7 MB across the
still-live parent and fresh child together, which is a different measurement.
This establishes representative cache correctness, not exhaustive coverage of
every admitted numerator or a claim about four-loop cache readiness.

== Default finite one-through-three-loop primers

The public-API driver subsequently completed the default relative-two-dot,
rank-one primers under explicit decimal 1 GB RAM and disk allowances, using
the same frozen native libraries as the passing four-loop gate. Its executable
SHA256 is `55cd4b5a690cdecd1fb97ec7cc8f63bf1363c53546b013a45632867a2f31bc81`.

#table(columns: 4,
  [Loops], [Completed roots], [Preparation call], [Charged disk bytes],
  [1], [5], [0.00922 s], [3,201],
  [2], [77], [0.01248 s], [4,514],
  [3], [2,320], [5.95493 s], [75,240],
)

Each run spawned a fresh child, verified index-only attachment, checked
prepared full-support and pinch witnesses at mass squared 1 and 4 with zero
IBP rules, and compared held-out fourth-power roots with ordinary evaluation.
The three-loop pinch witness is undotted. These driver comparisons are exact
symbolic-value checks; separate typed-record guard parity is established by
the preceding tests, not inferred from these values. Clear-to-disk reuse also
passed. Timings include program loading and planning but exclude construction
of Vakint; they are shared-host functional measurements, not controlled speedups
or a physics-coverage bound.

== Preparation and attachment

`Vakint::prepare_cache` takes a loop count, an output directory,
`RustRedCachePreparationOptions`, a cooperative cancellation flag and a progress
callback. Inspect the returned status: completion, time-budget stop,
cancellation, storage-budget stop, native resource refusal and failure are
distinct. A partial report does not certify a completed workload.

`options.storage.batch_records` is a preparation layout target, defaulting to
64 staged completed records across namespaces, including captured intermediates.
Python exposes the same optional
`batch_records` keyword. Positive values no larger than the independent native
chunk record limit are accepted; zero or an oversized target is rejected before
filesystem creation. Earlier byte/term or checkpoint limits may still flush a
smaller batch. This target changes neither the selected roots nor their native
admission, exact reduction, guards, or RAM/work allowances.

The reader's default hard limit remains 4096 records per chunk. Existing larger
chunks are therefore still loadable when the preparation target is reduced to 64.
Preparation resume additionally requires the current intermediate-capture profile.
Changing that target does not change a workload's identity,
rewrite completed chunks, or imply that an already loaded pack is repartitioned.
The batch-size choice does not change the on-disk schema. This separates a read-through cost choice from
the decoder's safety ceiling; it is not a new cache algorithm. The layout
implementation was tested at RustRed `d44e3beea9de1370a9cdd6fff08cbfbb486c1c4e`.
Its 53 native application tests passed, including exact record/ordered-guard/work
parity across batch sizes and reuse of a larger existing chunk.

`RustRedCacheWorkload::Explicit` supplies native integral keys with an explicit
namespace index; the exact program bindings supply their family identity.
`Finite` uses RustRed's finite A/R domain iterator. The default `Primer` selects
all physical parent pinches with at most two extra positive powers above the
surviving unit-power lines and numerator rank at most one. `RelativePrimer`
allows those two operational allowances to be changed explicitly. All counting,
pinch enumeration and power compositions are owned by RustRed; Vakint supplies
the exact program bindings and physical-slot masks. The default counts are
5, 77, 2320 and 164,352 roots for the shipped one-through-four-loop parents. These
are finite operational choices, not physical power-counting envelopes. Every root,
including one already stored during an earlier run, crosses the active native
program's admission before reuse. The planner errors on its root/request cap;
it never silently truncates the selected workload.
Execution uses deterministic namespace round-robin order, not a cost-ranked or
difficulty-ranked selection. A smaller chunk target does not select easier
integrals or enlarge the admitted workload.

Preparation does not attach its output. `Vakint::load_cache` explicitly attaches
a trusted-generated read-only bundle. Disk misses may fill the configured native
RAM cache but do not modify the pack. `clear_cache` clears RAM only, retaining
the backing and cumulative counters. `cache_statistics` distinguishes native
RAM hits, backing hits, actual reductions, eviction and bounded I/O activity.

The current local four-loop implementation shares immutable native programs
but gives each attachment isolated mutable sessions. It does not clear or warm
the ordinary static reducers during preparation or cached evaluation. Calling
statistics, clear or save without an attachment does not load the four programs.
Namespace details retain exact bindings; aggregate telemetry adds current
retained payload and work counters with saturation, never estimates process RSS.

== Explicit snapshots and optional reads

`Vakint::save_cache` requires an explicit `RustRedCacheSaveMode`.
`ResidentRoots` exports only completed records currently in RAM; never-loaded
and evicted backing records are deliberately excluded. `FullSnapshot` instead
copies every captured immutable base through bounded, checksum-verified byte
transfers without decoding coefficients, then appends the resident overlay.
It requires a new destination. Source paths can subsequently move without
invalidating a successfully completed portable snapshot.

Inspect the report, rather than treating a return value as completion.
Callback-accepted records, durably checkpointed resident records, and the
base-copy phase are separate. A successful base copy followed by an overlay
failure remains a partial snapshot. Counts of preserved base record references
include shadowed entries; they are not a distinct-key count. All namespaces
share one configured disk quota. Resident-only resume can add overlapping runs
and consumes quota; no implicit compaction or partial-base-copy resume occurs.

`RustRedCacheSaveOptions` reuses preparation's disk defaults without changing
the running cache's RAM, eviction or arithmetic policy. Load options independently
accept an optional lookup-budget policy; omitting it on reload preserves the
active policy. A new optional attachment defaults to cold native reduction only
when the typed prospective lookup-work allowance is exhausted, recording a
separate fallback counter. `Strict` reports that exhaustion. Corruption,
incompatible bindings, malformed records, failed native admission and arithmetic
limits never become misses. Preparation's own completed-root verification stays
strict; its optional descendant reads can fall back under the same explicit rule.

The small #link("examples/native_cache.rs")[Rust cache example] demonstrates
one-loop preparation, exact loaded/warm parity, a held-out reduction and a
full portable snapshot. It creates only caller-specified new directories and
does not evaluate new numerical masters. Its optimized local-patch runtime gate
passes. The first delivery selected published RustRed revision
`6a70dccfac11dbbaa9ffe3facfa5710b4f33a65f`; its actual-Git source check passed,
including the Python module and stub-generation feature combination. The subsequent
layout checkpoint selected the separately tested `d44e3bee` revision, before
the current descendant-reuse engine revision.
The default five-root one-loop primer used 3,201 disk bytes. A full snapshot
preserved five base references plus two resident records in 1,697 bytes, without
copying preparation-job metadata. Restored held-out I5 matched exactly with zero
native rules. GNU time measured 0.10 seconds total and 12,288 KiB per-process
maximum RSS; this small workflow is not a performance comparison.

== Python and HepKit steering

The current local `Vakint` Python wrapper exposes the same `prepare_cache`,
`load_cache`, `save_cache`, `clear_cache` and `cache_statistics` operations for
one through four loops, plus explicit `import_cache` and `download_cache` on
supported native platforms. The transfer facade is source-checked and its native
implementation passes eight focused tests and an independent repeat. The final
private host passes all eleven Python smoke tests, including explicit import,
and the nine public transfer cases described next. The earlier host receipts
below retain their original scope. Use a coherent Symbolica /
HepKit extension build, never a separately loaded second Symbolica kernel.
#link("examples/native_cache.py")[The Python example] uses only a tiny one-loop
workload; #link("tests/python/native_cache_smoke.py")[the host smoke tests]
cover callbacks, interrupts, GIL release, snapshots and partial native overrides.

=== Earlier October 14 transfer/runtime gate

That completed gate selected actual Git RustRed
`301beec66e6a2fcdc12ffe2c665181206615550c` consistently across the core,
application and FeynKit bridge. The private Community host contains one shared
Symbolica 3.0.1 kernel. Its wheel SHA-256 is
`d333ed26002bbb8e88465a639d292086ee017139abd0bfbff715bcae43565d20`;
the installed module SHA-256 is
`22b3eada65bc29fa2580effa3de409f4c57fe7928a2ba93444ec898f53a4a4a7`.
The supervised build completed in 1,418.49 seconds with a sampled process-tree
peak of 14.962 GB and unchanged captured source inputs.

All eleven Python smoke tests passed in 0.889 seconds. Nine separate public
cases imported and loaded the unchanged one-through-four-loop transport
archives, preserving 5 / 77 / 2,320 / 164,352 stored root references and exact
native chunk bytes and bindings. Prepared outputs and lower-loop held-out
controls agreed exactly. Four fresh four-loop held-out requests agreed with
both source-pack evaluation and native binary Atom goldens from the frozen
older host; descendant hits and exact-key read skips were positive. An explicit
loopback HTTP download also passed. Those nine public cases did not run preparation, modify
the source packs or archives, attach a cache implicitly, or change the frozen
older installation. These are functional controls, not general timing claims.

This private wheel requires the validated OpenSSL 3.6 loader environment used
by its test runner. An initial bare Nix Python invocation selected an incompatible
OpenSSL 3.5 library and failed before importing the module or running tests;
the failure receipt was retained. Successful tests used the matching loader
environment with the same wheel. These results do not establish portable-wheel
or default-shell import, a Community/PyPI release, or whole-workspace CI.

The final native transfer gate and its frozen-binary repeat each pass all eight
tests, including checksums/bindings, malformed archives, independent quotas,
atomic destination refusal, standard USTAR padding and bounded loopback HTTP.
Python-feature source checking and Clippy complete; existing warnings and the
nonfunctional explicit tar-borrow-ending warning remain visible. Hakari verifies
without changes, and one `just ci-update` run reproduces the existing graph and
configuration with CI disabled. No full CI suite or upload is claimed.

=== Earlier private host gates

The initial private complete Community wheel, SHA256
`8ddfee95fcf7cd4e28658f0614e3b20291f1a92a15aefdf168e88604f645676d`,
uses actual Git RustRed `6a70dccf`, Symbolica 3.0.1 and one shared native kernel.
All eight smoke tests pass in 0.476 seconds. The first attempt passed seven;
the remaining fixture incorrectly supplied a string instead of the existing
constructor's required Symbolica epsilon variable. Correcting only that input
made all eight pass, without changing native code or weakening assertions.
The public Python example prepared five roots (3,212 charged bytes), made a full
snapshot (1,716 bytes), and restored held-out I5 exactly with zero native rules.
The environment and wheel are private validation artifacts, not a Community
release or modification of a shared installation. No new master values were
computed. That earlier wheel gate did not establish completion of the full
default four-loop primer.

The incremental private host rebuilt against actual Git `d44e3bee`, wheel SHA256
`430276d333ec659e24fadcf8be7df152ebcaf9bf9f8d4ec30285074490aac380`,
passes all ten Python cache tests in 0.762 seconds. The added controls verify the default
64-record target, an explicit two-record target producing three chunks for five
roots, invalid targets rejected before output creation, and larger-chunk resume.
In a separate cross-version check, the preceding `6a70dccf` wheel produced one
65-record chunk with no batch-target field in its checkpoint options. The new
wheel loaded that unchanged pack and resumed a copy under the 64-record target:
65 records were reused and two were newly completed. The source pack's hashes
were unchanged. The public Python example also passed on the new host, using
3,231 bytes for preparation and 1,716 bytes for its full snapshot. Installed
native-module hashes were stable across all these checks. This is private
runtime verification, not a Community publication or a full four-loop warmup.

Preparation accepts `ram_bytes`, `disk_bytes`, `cache_policy` (`evict-lru`,
`evict-aged-lfu` or `fail-on-limit`), optional native term/integral/rule ceilings,
time/checkpoint/root allowances, the independent `batch_records` layout target,
and `resume`. The workload is either the native
relative primer controlled by `extra_dots` / `numerator_rank`, or an explicit
sequence of `(namespace, complete_integer_powers)` pairs; the two modes cannot
be mixed. No network access or numerical master computation is introduced.

The GIL is released during native work. Optional `progress(report)` and
`cancel()` callables run under the GIL only between native operations. A False
progress result or True cancellation predicate requests cooperative cancellation.
Progress callbacks are throttled (one second by default), while Python signals
and cancellation are checked at the runner's normal root boundaries. One native
operation can overrun the time allowance; this is not asynchronous interruption
inside the algebra engine.

Python callback errors and KeyboardInterrupt are re-raised after the native
runner finishes checkpointing or records a checkpoint failure. When the Python
exception supports an added attribute, `cache_report` contains that final report.
Always inspect `checkpoint_failure`; the latest checkpoint is not guaranteed
durable merely because a Python exception was preserved. Ordinary partial
resource/time/cancellation results remain report dictionaries, not false success.

Python load controls apply only explicitly supplied native overrides to the
current RAM/work limits. Omitted eviction and lookup policy also stay active.
Reader storage/I/O limits are different: the Rust load API uses the supplied
whole storage policy, or its defaults, on every load. In particular, omitted
Python `disk_bytes` selects the default reader disk allowance rather than
preserving a prior custom I/O policy. Full snapshots require
`mode="full-snapshot"`; a RAM-only export requires `mode="resident-roots"`.

The binding combines the full family identity, physical arity, cache semantic
layer and digest of the exact immutable bytes used to load the program. Neither
the path nor loop count is a compatibility identity. The one-through-three-loop
loader retains its SHA256 digest together with the artifact it decoded; it
does not reencode Symbolica state or reread a possibly replaced path later.
Content hashes detect corruption and incompatibility, not maliciously forged
mathematical identities: native generated-cache provenance remains required.

== Operational allowances

LRU is the default; aged LFU and fail-on-limit are explicit alternatives. The
ordinary evaluation defaults are unchanged. Preparation's decimal disk targets
are one billion bytes separately for one, two and three loops, and sixteen
billion bytes for the aggregate four-family bundle. These are maxima,
not minimum data quantities to manufacture. Its default RAM allowance is the
lesser of the disk target and thirty-two billion bytes. Byte counts use checked
portable integer conversions. An existing attachment keeps its active native
RAM limit unless load options explicitly override it.

RAM coefficient payload, bounded index/chunk/decode buffers and process RSS are
different quantities. Reports include the selected storage and native work
policies; increasing storage does not silently raise unrelated polynomial,
coalescing or rule-application limits. A genuine active-working-set refusal is
not converted into a cache miss or an out-of-scope result.

The default preparation time allowance is 36,000 seconds and is cooperative:
one noninterruptible native operation may overrun it. Completed root records are
published through bounded atomic child checkpoints. The binary progress record
shares the bundle's disk quota, including temporary replacement files and
interrupted tails. Resume fixes the artifact and workload identities but permits
an explicit change to operational budgets. It reuses completed guarded roots;
it does not restore an in-flight native stack. Both closed one-through-three-loop
programs and candidate four-loop programs can reuse completed guarded records
for matching descendants while reducing an unseen starting root. This lookup
belongs to RustRed’s native traversal, with source/applicability checks and all
imported guards preserved. Entry-only rank limits do not clip native descendants.
Preparation also streams newly computed admitted intermediate identities through
the same bounded writer. Their guards come from exact per-node provenance, not
coefficient-only memo entries. Native entry admission still limits which
intermediates can populate the external guarded-root namespace; it does not clip
the actual reduction. The existing directory-memory allowance bounds optional
provenance metadata per request. Exhaustion disables further capture for that
request, not its valid reduction; reports count these events. Oversized optional
records and inconclusive optional lookup grants are reported separately. Selected
roots become complete only when their own records are durable or validated from
an existing committed snapshot. A failed parent can therefore leave useful
durable descendants without falsely advancing selected-root completion.

`newly_stored_intermediates` records capture provenance at publication time. Some
of those keys may later be selected themselves, so it is not a count of keys
outside the complete primer. Ordinary `save_cache` exports resident completed
records or the explicitly requested immutable snapshot; it does not construct
missing provenance for arbitrary transient coefficient-only memo entries.
LRU/LFU governs
the outer completed-record cache. The transient candidate inner DAG still has
no LRU/LFU eviction of its own intermediate nodes.

Selected-root and captured-write progress counters are live and durable. During
an active preparation, the detailed native statistics, charged-disk census and
reader-directory bytes in progress callbacks are currently setup/final snapshots,
not per-callback live measurements. Quotas remain enforced by the writer, and
the final report refreshes these values. Use independent process/RSS and native
directory inspection for live resource headroom; a stale progress census must
not be interpreted as actual current disk use.

The current local opt-in request cap is one million, allowing the native
four-parent relative-two-dot/rank-one plan of 164,352 requests. An explicit
100,000 cap still refuses that plan before execution; no workload is silently
clipped. Both the earlier root-only workload and the new enriched full four-loop
workload have completed, with separately retained packs and measurements.
The first-38 and four unseen-input exact checks establish their stated finite
validation scope; neither preparation completion nor those witnesses certify
all possible reduction paths or an unrestricted integral domain.
