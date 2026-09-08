# IR-safe threshold update and rebase plan

## Goal

Rebase `ir_safe_threshold_subtraction` onto
`origin/codex/raised_energy_cff_wip_optimized_phase_fix`, retain the target
branch's raised-energy CFF, three-dimensional-representation, archive, and CLI
architecture, and correct IR-safe threshold grouping.

The runtime must solve only compatible thresholds together. For every solve
group, maximal overlap discovery, E-surface classification, SOCP construction,
center assignment, radial-root construction, and threshold multichanneling are
performed independently. Groups are then recombined only when forming the final
counterterm value.

## Group identity and directives

For a fully resolved threshold subspace, define its canonical solve-space
signature as the ordered formal loop variables and their defining edges in the
parent LMB used by the SOCP:

```
solve_signature = ((k_1, e_1), (k_2, e_2), ...)
```

The graph/cograph LMB used to express external sampled data is not part of this
signature. Thresholds with the same formal loop variables but different
defining edges have different signatures and cannot share a solve group.

Automatic grouping uses the solve signature. An explicit `group_id` is an
additional subdivision label:

- equal signature and equal group ID: one group;
- equal signature and different group IDs: separate groups;
- different signatures and different group IDs: separate groups;
- different signatures and the same explicit group ID: a clean resolution error,
  never an implicit merge.

The error must be deterministic and exhaustive. It will report the graph, side,
explicit ID, threshold variant names/IDs, physical cut and cut-group IDs,
threshold/cut edge sets, requested subspace and parent LMB, and the conflicting
canonical `(loop variable, defining edge)` signatures.

Group IDs are serialized on threshold-counterterm variants and propagated
through parsing, validation, resolution, materialization, runtime metadata,
archives, display, and Python APIs. Cross-section IDs are namespaced by side;
amplitude IDs are graph-local. Omitted IDs use automatic grouping and never join
an explicitly tagged group. IDs start at zero and are validated consistently in
each namespace.

When no threshold metadata is present in a loaded DOT graph, preserve the
legacy default: all left thresholds of a cut use the maximal available left
loop subspace and therefore form one left group; all right thresholds use the
maximal available right loop subspace and therefore form one right group.

## Cross-section kinematics and shared groups

Cross-section evaluation must not construct a shared SOCP before the cut
kinematics exist. The graph evaluation will be split into explicit phases:

1. For every active Cutkosky cut group, solve its kinematics, including all
   required `t^\star` rescalings and derivative/residue data, and retain the
   resulting representative sample and external data.
2. Classify the threshold instances using each cut's solved kinematics and
   cut-local reversed-edge rules.
3. Group compatible instances across cut groups within the same side. Each
   group carries the per-instance solved kinematic/external data needed to
   construct its E-surface constraints, while sharing only the compatible
   solve-space coordinates and center problem.
4. Build one group-local maximal-overlap/SOCP structure and one threshold
   multichanneling structure for each group. A group with several overlap
   centers gets its own center channels and PDFs.
5. Remap each shared group result back to the local cut-group threshold IDs, then
   run the existing residue, evaluator, CFF, and iterated-counterterm logic with
   its original CutGroupId indexing.

Left and right sides remain independent. Thresholds in different solve-space
signatures remain independent even if they have the same numeric group ID.

## Rebase procedure

- Fetch and verify the remote target and common ancestor.
- Preserve a recovery ref for the pre-rebase branch.
- Replay all existing branch commits with `git rebase --onto` from the common
  ancestor onto the target tip.
- Resolve shared files by preserving target CFF/3D-representation/archive and
  command-template changes, then restoring threshold directive, metadata,
  multiplier, and runtime behavior around them.
- Keep target-only files and dependencies intact.
- Run formatting and compilation after the core threshold replay and after the
  complete commit series.

## Implementation areas

- Add `group_id` to the graph-level threshold variant schema and all resolved
  threshold draft/metadata structures.
- Add a canonical solve-signature formatter/comparator on `SubspaceData` using
  formal loop variables and defining edges, with diagnostics suitable for
  incompatible explicit IDs.
- Validate explicit group-ID compatibility after all parent-LMB and subspace
  resolution/rebasing is complete, before threshold generation or runtime setup.
- Partition amplitude variants before existence classification and overlap
  construction; give every group independent overlap, SOCP, centers, radial
  roots, multiplier workspace, and multichanneling.
- Add cross-section group metadata/cache construction after all cut kinematics
  are solved, with per-instance external data and local-ID remapping.
- Keep evaluator/archive indexing and legacy no-directive paths stable except
  where version markers must change for new serialized metadata.

## Tests and physics validation

Add focused tests for schema round trips, automatic signatures, same-cycle /
different-edge rejection, explicit-ID incompatibility diagnostics, omitted IDs,
side-local namespaces, and independent multichanneling.

Extend the existing GL297 and GL638 tth@NNLO fixtures and command cards. For
both graphs, verify a soft-limit cancellation in which threshold counterterms
from different cuts share one compatible solve signature and are solved in one
common group, while thresholds with different signatures remain separate.
Verify the cancellation after the per-cut `t^\star` kinematics and external
data have been prepared. Retain the existing forced-cut, full-cut,
serialization, no-directive-equivalence, raised-energy CFF, archive, and overlap
regressions.

Run `cargo fmt`, targeted `cargo check`, focused `cargo nextest` tests, broader
workspace checks, and clippy.
