# IR-safe threshold grouping and subtraction plan

This document records the intended IR-safe threshold behavior and the implementation
boundaries used while rebasing the threshold work.

## Group identity

Threshold counterterms are grouped by the complete loop-variable subspace on which
they are defined. The identity includes the formal loop cycles and the selected graph
edges that realize those cycles. Thresholds from different Cutkosky cuts may share a
group when this identity is equal; thresholds with different subspaces must never
share one, even when their formal cycle notation looks similar.

When no threshold metadata is supplied in a loaded dot graph, each side of a cut
defaults to its maximal available loop subspace. Consequently all left thresholds of
that cut and all right thresholds of that cut are in the same automatically derived
group.

An optional user `group_id` starts at zero and can split an automatically compatible
group into smaller groups. It can never merge incompatible subspaces. Reusing one
explicit ID for incompatible thresholds is an error. The diagnostic must identify
every conflicting threshold, cut, amplitude or cross-section graph, its loop-cycle /
edge realization, the derived subspace, and the requested group ID.

## Per-group solve pipeline

For every group independently:

1. Resolve and validate all participating cuts and their threshold metadata.
2. Solve the kinematics for every participating cut, including each cut's
   `t^star` rescaling, before constructing any shared SOCP input.
3. Build the complete maximal-overlap structure over the group's e-surfaces.
4. Construct and solve one SOCP problem for that group to obtain its common center.
5. Assign the resulting center only to thresholds in that group.
6. Perform the group's MC-over-e-surfaces and integration setup using the same group
   boundary; do not share centers or overlap state with another group.

This order is required when thresholds from different cuts occur in one group: the
SOCP external data are valid only after all cut-specific kinematics and rescalings
have been solved.

If a group's maximal-overlap structure has two or more centers, multichanneling over
the threshold must also be built and run separately for that group. Channels and
weights must never combine thresholds from distinct groups.

## Defaults and compatibility

The automatic grouping path is used when `group_id` is absent. Explicit IDs are
validated against the automatic subspace identity before any solve starts. Existing
graphs without metadata retain the maximal-space defaults above, preserving the
historical behavior while making the grouping explicit and deterministic.

## Verification

Validation must cover the two established `tth@NNLO` graphs (GL297 and GL638):

- thresholds from different cuts that share a subspace are solved in one group;
- thresholds in different subspaces receive independent SOCP centers;
- explicit incompatible IDs fail with exhaustive diagnostics;
- nontrivial overlap groups build independent multichannel structures; and
- soft-limit evaluations show the expected IR cancellation across thresholds from
  different cuts when they belong to the same group.

The filtered state-loading interface implemented alongside this work is orthogonal:
it selects saved processes and integrands before integrand deserialization, forces a
read-only, no-save session, and keeps process IDs dense within the selected view.
