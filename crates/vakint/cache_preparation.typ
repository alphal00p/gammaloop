= Explicit native cache preparation

This interface is opt-in. Constructing Vakint and evaluating without an attached
cache keep the existing behavior; no cache files are read, written or downloaded
automatically. Cached payloads are exact, unit-mass, unnormalized reductions and
their complete native nonzero conditions. They are not numerical master values,
epsilon expansions, or a claim of additional reduction coverage.

== Current integration status

The one-through-three-loop sealed-program adapter and Rust preparation/load
surface are implemented and source-checked. The real one-loop fresh-process
acceptance test passes on RustRed `59e45868`; broader two-/three-loop cache
comparisons are pending. Four-loop candidate-session integration, the Python
surface and explicit download support are follow-up work; do not treat this
development interface as completed one-through-four-loop cache delivery.
Five-loop cache preparation is deferred until the updated artifact is available.

This checkpoint uses the first development cache-pack format. The forthcoming
mandatory compressed format deliberately changes its schema; previously made
development packs must be prepared again. No stable on-disk migration contract
is promised here, and an incompatible pack is refused rather than guessed at.

== Measured one-loop gate

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

== Preparation and attachment

`Vakint::prepare_cache` takes a loop count, an output directory,
`RustRedCachePreparationOptions`, a cooperative cancellation flag and a progress
callback. Inspect the returned status: completion, time-budget stop,
cancellation, storage-budget stop, native resource refusal and failure are
distinct. A partial report does not certify a completed workload.

`RustRedCacheWorkload::Explicit` supplies native integral keys with an explicit
namespace index; the exact program bindings supply their family identity.
`Finite` uses RustRed's finite A/R domain iterator. The small default `Primer`
is an operational choice, not a physical power-counting envelope. Every root,
including one already stored during an earlier run, crosses the active native
program's admission before reuse. The planner errors on its root/request cap;
it never silently truncates the selected workload.

Preparation does not attach its output. `Vakint::load_cache` explicitly attaches
a trusted-generated read-only bundle. Disk misses may fill the configured native
RAM cache but do not modify the pack. `clear_cache` clears RAM only, retaining
the backing and cumulative counters. `cache_statistics` distinguishes native
RAM hits, backing hits, actual reductions, eviction and bounded I/O activity.

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
billion bytes for the eventual aggregate four-family bundle. These are maxima,
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
it does not restore an in-flight native stack or promise lower-loop descendant
reuse inside an unfinished root.
