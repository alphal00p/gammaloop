= Three one-diagram QED numerator oracle migration

== Scope

This bounded migration updates exactly three numerator digest literals in
`crates/feynkit-generator/tests/generation_oracle/main.rs`. It is a separate
descendant of the independently approved certificate at `16bfd2c8`; that
certificate and all earlier diagnostic rounds remain frozen. No production
algorithm, request, option, filter, class mode, naming, lock, kernel, project
manifest or NixCI setting changed. The other 78 failing rows and any broader
physical claim remain outside scope.

#table(
  columns: (auto, auto, auto),
  table.header([Case], [Old numerator digest], [New numerator digest]),
  [`xs_qed_eedd_L1`], [`f5a21f1a84870bc7`], [`d3ad2f3a0802b6ab`],
  [`xs_qed_eemm_L1`], [`8b7b40bba961420d`], [`414e2a7a954beae8`],
  [`xs_qed_eeee_1pi_L1`], [`0db0b7cc52558270`], [`1fa84d9ad4fdf0b0`],
)

The old and new sources have the same line count and differ on only lines 331,
337 and 347. Each changed line is identified by its exact named row and differs
only in the `numerators` hexadecimal literal. Counts, inverse symmetry,
rational weight, symmetry histogram, content digest, name digest and ID digest
remain byte-identical. No record, update or force-pass mode was used.

== Producers and actual pool runs

Before rebuilding, the old-literal current producer was copied to
`target/oracle-numerator-migration/old/generation_oracle`, SHA-256
`c2ef6aea...be84`, together with source SHA-256 `499b49e4...8aa8`. It is the
same qualified current producer used by the certificate. With explicit
`FEYNKIT_ORACLE_THREADS` values 1, 8 and 32, it executes three named tests per
run and fails all nine executions at the unchanged old expectations. Every
failure reports only the numerator digest difference and the same current
values shown above.

The new source SHA-256 is `987e9c27...049b`. A fresh locked native-default
`dev-optim` build produced ELF SHA-256 `8cde4358...96a`. That retained binary
executes three nonzero named tests separately at pools 1, 8 and 32; all nine
executions pass. The owner reads `FEYNKIT_ORACLE_THREADS`, iterates the requested
counts, and inside each count runs `Classes::None` followed by
`Classes::Auto`. Distinct `t1`, `t8` and `t32` inventories prove that pools 8
and 32 were actually run rather than inferred from pool 1.

For every case and pool, None, Auto and their `sources.json` companions are
byte-identical. Old/new and cross-pool payloads and census text are also
byte-identical. An independent Summary digest replay reproduces content,
numerator, name and ID digests for all 18 route/pool/case records. Passing the
new full `Summary` assertion additionally checks unchanged counts, rational
weights and symmetry histograms. This is fresh pool-8/32 regression evidence,
not an extension of the pool-1 symbolic equivalence certificate.

Both run scripts explicitly unset record, diagnostic and common update/force
variables, print only nonsecret environment and argv, and use the private
license without storing its value. The old and new ELFs, sources, pool logs,
dumps and census files remain in the task-owned target tree. The migration
#link("evidence/generation-oracle-numerator-migration/receipt.json")[receipt]
and #link("evidence/generation-oracle-numerator-migration/validate.mjs")[validator]
bind their hashes and exact commands.

== Certificate reuse and limits

The approved round-four combined review has verdict `APPROVE_CERTIFICATE` and
SHA-256 `1202230f...e4b30`. The frozen round-four capture and semantic validator
was replayed against its original inputs before changing the active registry;
it still reports six capture rows, 92 typed axes, six transports and all
approved controls. The migrated expectations do not feed that checker and are
not treated as an independent oracle.

The approval is fixed-case equivalence for these three one-diagram QED rows at
pool 1. Fresh pools 8 and 32 establish only regression consistency for the same
requests and producer. Nonblocking certificate notes remain qualified: the
semantic theorem acts on parsed JSON while raw-byte identity binds duplicate
keys; some leaf values remain capture-bound; and no whitespace no-op or process
exit status is counted as broader mathematical evidence.

Rustfmt with explicit edition 2024 and scoped Clippy pass. An initial direct
Rustfmt attempt without the edition failed while parsing the Rust-2024 let-chain
in the unchanged census module; it receives no formatting credit and is
retained as a failed attempt. The registered Typst source, PDF and semantic HTML
bundle are hash-bound. The full documentation builder remains blocked by the
inherited eight tracked mark-geometry HTML files plus two retained review HTML
files; none were removed.

No integration, main edit, publication, push, fetch, pull request, upload or CI
activation occurred. This migration stops for a new independent Opus migration
and test review, with Astra fallback authorized.
