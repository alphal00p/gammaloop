= Explicit native cache preparation

This interface is opt-in. Constructing Vakint and evaluating without an attached
cache keep the existing behavior; no cache files are read, written or downloaded
automatically. Cached payloads are exact, unit-mass, unnormalized reductions and
their complete native nonzero conditions. They are not numerical master values,
epsilon expansions, or a claim of additional reduction coverage.

== Current integration status

The one-through-three-loop sealed-program adapter and Rust preparation/load
surface pass real fresh-process gates on RustRed `a0649fed`, including exact
mass restoration, guard preservation and held-out native reductions. Four-loop
candidate-session integration, the Python
surface and explicit download support are follow-up work; do not treat this
development interface as completed one-through-four-loop cache delivery.
Five-loop cache preparation is deferred until the updated artifact is available.

The current RustRed dependency uses mandatory gzip chunks and cache-pack/bundle
schema 02. The earlier `59e45868` gate below used schema 01: its development
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

== Current one-through-three-loop gates (schema 02)

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

== Preparation and attachment

`Vakint::prepare_cache` takes a loop count, an output directory,
`RustRedCachePreparationOptions`, a cooperative cancellation flag and a progress
callback. Inspect the returned status: completion, time-budget stop,
cancellation, storage-budget stop, native resource refusal and failure are
distinct. A partial report does not certify a completed workload.

`RustRedCacheWorkload::Explicit` supplies native integral keys with an explicit
namespace index; the exact program bindings supply their family identity.
`Finite` uses RustRed's finite A/R domain iterator. The default `Primer` selects
all physical parent pinches with at most two extra positive powers above the
surviving unit-power lines and numerator rank at most one. `RelativePrimer`
allows those two operational allowances to be changed explicitly. All counting,
pinch enumeration and power compositions are owned by RustRed; Vakint supplies
the exact program bindings and physical-slot masks. The default counts are
5, 77 and 2320 roots for the shipped one-, two- and three-loop parents. These
are finite operational choices, not physical power-counting envelopes. Every root,
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
