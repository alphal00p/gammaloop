= Factorized numerator restoration validation

This record validates the nine-change raised-energy stack after its Symbolica 3
rebase. Measurements were taken on 2026-09-17 at
`ac2c553c36309e02529574812e70d4bdaa37eba8`. The accompanying documentation update
does not change its production sources, dependencies, or benchmark inputs.

== Representation and correctness

Branch numerators remain shared functions of affine energy arguments. A branch
supplies the coefficients of internal on-shell energies, external energies,
the integer sampling scale, and a constant. The complete source energy map
defines branch identity; a sign alone does not. Arguments can span the complete
graph, including energies external to a UV component. This accommodates future
affine residue maps without implementing LTD.

The existing `NumeratorAtomExt::series_preserving_factors` owns factor-preserving
Taylor expansion. It separates factors independent of the expansion variable,
retains shared numerator coefficient families, and delegates the required Laurent
series algebra to Symbolica. Energy arguments that depend on the expansion
variable remain visible to differentiation. Both UV routes use this machinery.

`Integrands` carries cut expressions and shared definitions separately. Evaluator
preparation contracts each numerator body once, keeps its branch arguments on
component calls, preserves sparse components before outer contraction, and
registers reachable functions and scalar aliases. Explicit orientation sums
retain products and exact cancellation. Independent UV component sums cannot
always be separated: derivatives can couple their energy arguments. The shared
bodies survive that combined sum without copying and expanding each numerator.

Regression coverage includes independent scalar/tensor factors, Laurent and
zero series, full affine rows with independent sampling nodes, exact summed
cancellation, inactive singular definitions, alias collisions, sparse component
lowering, and standalone archive export/reload. The owning tests live beside
`numerator/symbolica_ext.rs`, `uv/approx/direct_3d/branches.rs`, and
`integrands/process/evaluators.rs` in `crates/gammalooprs/src/`.

== Comparison with main

Live main was verified as `312cf6aaf4cd1414f6742ad87ce39922c497a0b6`.
The preserved native main binary was built from
`226b7bea45398023f767fb4e53a1b5c6c3a43a2d`, through its unchanged child checkout
`6fbfe7cd89839ca6b9f10bdd1e5bbc9eba1152f1`. The eight files changed between that
reference and live main concern CI and contributor guidance; Rust sources,
assets, Cargo manifests, lockfile, and toolchain are identical. The archive lists
those paths. This is reuse of a verified equivalent native binary.

Both binaries use Rust 1.98.1, Symbolica 3.0.0, and the `dev-optim` profile. The
candidate uses `intermediate_cost`; main uses its available `sparse_atom_aware`.
This comparison intentionally uses the established preferred executor on each
side. Compilation is excluded. Each entry is one fresh run, with no retries.

#table(
  columns: (2fr, 1fr, 1fr, 1fr),
  table.header([Native generation], [Main], [Candidate], [Main/candidate]),
  [GL1, UV disabled], [0.284 s], [0.092 s], [3.09],
  [GL0 box, logarithmic UV only], [4.278 s], [1.915 s], [2.23],
)

The no-UV graph has a numerator multilinear in its edge energies. The box has
superficial degree zero and no proper loop subgraphs, so its local UV subtraction
requires no positive-order Taylor derivatives. These are eligible numerical
comparisons within main's supported numerator powers. The frozen graph and
Taylor-order certificates are included in the archive.

Each side performed two identical inspections with finite, nonzero, repeatable
results. Persisted runtime settings and model parameters agree. Applying the
previously established phase conversion `(-1)^L` gives relative differences
`8.03e-14` for GL1 and `2.37e-16` for the box. The archive retains raw complex
values, rather than only the normalized differences.

Whole-command wall times, which also include startup, inspection, and exports,
were 1.259/1.274 s for main/candidate GL1 and 5.144/2.977 s for the box. Thus the
short no-UV command is approximately at parity. Sampled peak process-tree memory
was 49.2/48.8 MB and 89.4/139.7 MB respectively; the box uses more memory.

The native preparation bucket decreased from 0.198 to 0.043 s for GL1 and from
3.975 to 0.913 s for the box. Symbolica evaluator construction itself increased
from 0.012 to 0.018 s and from 0.197 to 0.922 s. The improvement is in total native
generation; it does not establish that every stage is faster.

== Existing derivative-UV workloads

The established GL1 derivative workloads also completed with the candidate's
`intermediate_cost` setting:

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  table.header([Local UV route], [Native generation], [Command wall], [Sampled peak]),
  [Direct 3D], [100.382 s], [101.568 s], [2.092 GB],
  [Projected 4D], [10.655 s], [11.863 s], [1.351 GB],
)

These cards measure generation only. Their timings are not an equal-physics
comparison with main: derivative counterterms require higher numerator energy
powers that main does not support, and main has no equivalent projected 4D route.
Numerical correctness is covered by the local UV regression suite.

== Evidence and reproduction

The #link("raised-energy-cff-restoration-evidence.tar.gz")[evidence archive] and
its #link("raised-energy-cff-restoration-evidence.json")[SHA-256 manifest] contain
the six run receipts, native generation reports, inspection events, effective
cards, model parameters, persisted settings, graph exports, and frozen input
bundles. Source, binary, and input hashes were checked before and after every
run. Full native binaries and saved binary states are excluded.

Build a matching checkout in its licensed development environment:

```sh
NO_SYMBOLICA_OEM_LICENSE=1 cargo build --locked -p gammaloop-api \
  --bin gammaloop --profile dev-optim
```

Keep `SYMBOLICA_LICENSE` in the environment. Preserve the executable with a
`BINARY.build.json` receipt containing its SHA-256 and the checkout commit,
`Cargo.toml` SHA-256, and `Cargo.lock` SHA-256; `current-build.json` provides the
exact schema used. The runner rejects a receipt that does not match the source
and executable. After extracting the archive outside the checkout, run:

```sh
TASK_BUNDLE=/path/to/extracted/inputs/eligible-controls
python3 "$TASK_BUNDLE/tools/run-card.py" \
  /path/to/matching/checkout /path/to/preserved/gammaloop \
  "$TASK_BUNDLE/cards/candidate/g_gl1_no_uv.toml" \
  /path/to/fresh/external/output 900
```

Use the adjacent box card and the `cards/main/` controls for the other eligible
runs. The two UV cards are under `inputs/uv-routes/cards/candidate/`, named
`g_gl1_local3d_intermediate.toml` and `g_gl1_local4d_intermediate.toml`, with their
own identical runner interface. Each runner retains the established 30 GB guard,
900-second timeout, one generation core, eight Rayon threads, and 128 MiB worker
stack. Archived README preparation notes describe their original campaign;
the measured revisions and completed results are those recorded here and in
`comparison.json`.

On the measured revision, `just ci-checks-and-upload` passed: 2,784 selected
tests passed, none failed, and 301 were skipped. Formatting, Clippy, doctests,
and workspace/CI graph checks passed. Check derivations and the successful
check/upload receipt are retained under `validation/`. The required remote
NixCI deployment also passed on that revision and reused the uploaded outputs.
The final documentation update follows the same pre-push check/upload workflow.

These are single observations, not timing distributions or evaluator-throughput
measurements. Factor preservation is tested structurally, but this campaign does
not isolate the speedup due to the series method or parametric representation
from the other changes in the stack. Older run cards and runtime measurements
remain historical evidence; they were not rerun in this six-run validation.
