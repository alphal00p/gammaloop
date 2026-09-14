= BNL Scalar Alias Capture

#quote(block: true)[
#strong[Status:] Archived raw-output evidence

This record describes an undated diagnostic capture present at evidence revision
`4fdbf430b29edd24b9e1292c07ab24a151427dc7`. It is not a benchmark and does not establish current
runtime behavior. The four complete expression payloads occupy a 1.9 MB non-prose data artifact,
so the developer page remains navigable.
]

== Provenance

- Source atom:
  #link("../../examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_unfiltered_pre_network.sym")[`bnl_integrated_evaluator_atom_unfiltered_pre_network.sym`]
- Rendering: `log_print(Some(120))`
- Scalar-alias threshold: `4096` bytes
- Captured aliases: four, each with one term
- Artifact extent: `9481` lines, `1977598` bytes, including `230032` ANSI escape bytes
- Source-atom SHA-256: `01de1b29c85545ae1234c1c68674825e9950a2bcd075871a5479710fd816570c`
- Migrated artifact SHA-256: `63e22c639d4021e247ed40cb05247fa135eb02f27937aef913b48527402cf684`
- Pre-migration Markdown SHA-256: `f919b62ff11fc5afa8d7224500fce65b923ffbe667ac508480549f7404dfe4f9`
- Raw artifact:
  #link("../../examples/cli/BNL/profiling/bnl_scalar_alias_captures.ansi.txt")[`bnl_scalar_alias_captures.ansi.txt`]

The migration removed only the Markdown wrapper and added machine-readable capture headers; the
expression payloads retain their original ANSI SGR color sequences. View the artifact with a
terminal that understands ANSI escapes or strip them in a disposable copy. Do not treat the
printed byte counts as stable across Symbolica versions or formatting changes.

== Captured aliases

#table(
  columns: (1fr, 1fr, 1fr),
  table.header([*Alias*], [*Terms*], [*Reported bytes*]),
  [`scalar(0)`], [`1`], [`5272`],
  [`scalar(2)`], [`1`], [`59391`],
  [`scalar(3)`], [`1`], [`386472`],
  [`scalar(8)`], [`1`], [`953984`],
)

The size distribution explains why aliasing materially changes the evaluator-network input: the
largest captured scalar dominates the raw printed payload. This record preserves the observation
without embedding thousands of terminal-output lines in authored documentation.

For the related threshold-8 timing summary over the same input atom, see
#link("bnl-evaluator-atom-capture.typ")[BNL evaluator-atom capture].
