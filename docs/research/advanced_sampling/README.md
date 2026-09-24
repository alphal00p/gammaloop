# Advanced sampling research archive

This archive contains the research supporting
[ADVANCED_SAMPLING_PLAN.typ](../../../ADVANCED_SAMPLING_PLAN.typ), subsequent
implementation audits, and explicitly labelled numerical evidence. The plan
records current implementation status; a proposal or arithmetic check alone is
not evidence that the corresponding production feature works.

## Reading order and precedence

Agents implementing or auditing this work should consult the plan and relevant
research again at each implementation slice, before changing an invariant, and
when investigating a mismatch. Do not infer current requirements from an older
example in isolation. Where proposal details differ, use this precedence:

1. [ADVANCED_SAMPLING_PLAN.typ](../../../ADVANCED_SAMPLING_PLAN.typ): consolidated
   requirements, milestones and acceptance criteria.
2. [SOPER_AND_SAMPLING_API.typ](SOPER_AND_SAMPLING_API.typ): latest agreed string
   language, automatic presets, conditional existence, eager/dual derivatives,
   exact/proxy partitions and saved-state acceptance harness.
3. [API_AND_GENERICITY_ADDENDUM.typ](API_AND_GENERICITY_ADDENDUM.typ): genericity,
   amplitudes and physical-cut/left/right composition.
4. [REPORT.typ](REPORT.typ): original mathematical study, source audit and GL638
   derivations. Its original API syntax and staged scope are historical where
   the documents above supersede them.

The focused agent studies preserve independent reasoning and review:

- [LU_H_MATCHED_SAMPLING.typ](LU_H_MATCHED_SAMPLING.typ): current LU-scale measure,
  invertible h-matched proposals, raised-cut derivative coverage and finite-variance
  reference requirements. This extends the earlier cut-channel proposals.
- [AMPLITUDE_BENCHMARK_CANDIDATES.typ](AMPLITUDE_BENCHMARK_CANDIDATES.typ): generated
  UV-finite massive two-loop kite and double-box candidates, verified threshold
  intersections, production coverage gaps and correctness/variance gates.
- [DOUBLE_BOX_GEOMETRY_CHECK.typ](DOUBLE_BOX_GEOMETRY_CHECK.typ): independent
  routing/rank/center checks and an actual same-orientation A/B intersection.
- [WARMUP_CACHE_DESIGN.typ](WARMUP_CACHE_DESIGN.typ): implemented cache ownership,
  invalidation and worker-local eager buffers, with validation scope.
- [SAMPLING_PRECISION_RESCUE.typ](SAMPLING_PRECISION_RESCUE.typ): native map/density
  migration and original-source rescue required before strong-focus long runs.
- [CROSS_SECTION_DELIVERY_PATH.typ](CROSS_SECTION_DELIVERY_PATH.typ): conditional
  host/block proposal, parallel cut-h/side/star work and early GL638 milestones.
- [AFFINE_STAR_IMPLEMENTATION_AUDIT.typ](AFFINE_STAR_IMPLEMENTATION_AUDIT.typ):
  unimplemented X4 dependency certificate, shared physical centers/projection
  scale, canonical center branches and minimum validation gates.
- [PINCHED_COLLINEAR_SAMPLING.typ](PINCHED_COLLINEAR_SAMPLING.typ): future light-cone
  and prolate charts, normalized transverse enhancement, soft endpoints and
  architectural requirements for cut-dependent pinched geometry.
- [api.typ](api.typ): existing code owners and estimator boundaries; API sketches
  are historical.
- [gl638.typ](gl638.typ): GL638 geometry, projected targets and numerical scope.
- [literature.typ](literature.typ): primary literature and mathematical derivations.
- [SOPER_ANGULAR_REVIEW.typ](SOPER_ANGULAR_REVIEW.typ): angular and correlated density
  construction.
- [TWO_NORMAL_PROPOSAL.typ](TWO_NORMAL_PROPOSAL.typ): exact H/Z chart and A-star
  pullback derivation, with explicit regularity and branch conditions.

## Numerical evidence and portability

The small records included here are:

| File | Meaning |
|---|---|
| [lu_h_profile_checks.py](lu_h_profile_checks.py), [output](lu_h_profile_checks.json) | Portable one-dimensional LU h-profile variance and peak-weight comparisons; not a GL638 integration |
| [geometry_checks.json](geometry_checks.json) | Ellipsoid forward/inverse, independent Jacobian, volume and moment checks |
| [gl638_geometry_checks.json](gl638_geometry_checks.json) | Recorded GL638 rank, LU radial invariance and A-star affine arithmetic |
| [normal_plane_moments.csv](normal_plane_moments.csv) | Analytically normalized model comparisons for a codimension-two singularity |
| [hard_hz_points.json](hard_hz_points.json) | Frozen inputs and derived geometry for the H/Z approach |
| [a_star_pullback_points.json](a_star_pullback_points.json) | Frozen actual center, complements and A-star points |
| [AMPLITUDE_CHECKPOINT06D_SMOKES.json](AMPLITUDE_CHECKPOINT06D_SMOKES.json) | Six all-orientation kite/box smoke runs at a fixed checkpoint; no variance claim |
| [GL638_CHECKPOINT06D_BASELINE.json](GL638_CHECKPOINT06D_BASELINE.json) | Fresh full-UV baseline, all 936 orientations, six cuts and native registry/reload checks; no sampling improvement claim |
| [AMPLITUDE_X1_KITE_MATRIX.typ](AMPLITUDE_X1_KITE_MATRIX.typ), [data](AMPLITUDE_X1_KITE_MATRIX.json) | All-18-orientation normalized-reference and three-seed physical comparison; mixed channels improve the real-component pilot, surface alone worsens it |
| [AMPLITUDE_X1_DOUBLE_BOX_MATRIX.typ](AMPLITUDE_X1_DOUBLE_BOX_MATRIX.typ), [data](AMPLITUDE_X1_DOUBLE_BOX_MATRIX.json) | All-98-orientation reference and three-seed physical comparison; gains depend on the selected complex component |
| [GL638_X1_PILOT.typ](GL638_X1_PILOT.typ), [data](GL638_X1_PILOT.json) | All-936-orientation, full-UV, 20-core matched radial pilot; no consistent three-seed LU-h gain |
| [GL638_X2_DIRECT_H.typ](GL638_X2_DIRECT_H.typ), [data](GL638_X2_DIRECT_H.json) | Actual selected/full936-catalogue conditional-H reference, native inverse rays and full12D determinant; physical comparisons remain separate |
| [GL638_X2_PHYSICAL_PILOT.typ](GL638_X2_PHYSICAL_PILOT.typ), [data](GL638_X2_PHYSICAL_PILOT.json) | All-936-orientation full-UV pilot, native physical rays and 48 exact maximum replays; direct H remains unbounded and global efficiency is inconclusive |
| [GL638_X4_ALPHA_REPLAY.typ](GL638_X4_ALPHA_REPLAY.typ), [data](GL638_X4_ALPHA_REPLAY.json) | Three full-orientation native hard/soft replays after the alpha projection change; tested physics agrees, with no precision or variance gain claimed |
| [GL638_HOSTED_JOINT_GATE.typ](GL638_HOSTED_JOINT_GATE.typ), [artifacts](gl638_hosted_joint_gate/) | Original hosted H/Z joint maps: 53 native physical checks pass; warmed 20-worker cost exceeds 10%, and one retained source fails Quad/Arb stability; no integration gain claimed |

The frozen inputs and equations preserve the numerical setup without requiring
the original machine. They are not executable acceptance tests. The original
`gl638_geometry_checks.py` is deliberately excluded because it imports local
scratch helpers and reads machine-specific files. Its output is retained as
historical evidence; implementation must add independently runnable tests.

Existing physical integration evidence is already tracked in
[the GL638 research records](../gl638/README.md), especially
[the runtime follow-up](../gl638/runtime-followup.typ) and its portable replay
inputs. These integrations precede the proposed maps and are not proof that
those maps improve GL638.

Markdown links to files in this repository or this archive have been made
relative. Historical source line numbers in prose, code snippets and labels
refer to the audited source version and may drift. Other absolute paths in
Markdown or JSON are provenance pointers to optional external local artifacts,
not portable dependencies or instructions to recreate those paths. In
particular, the user-supplied `ttH_defo.pdf` and thesis in `TMP_TO_IGNORE/`, and
the original prepared GammaLoop states, are not bundled. Public papers remain
available through their cited URLs. No PDFs or extracted full paper texts are
included.

## Provenance

[MANIFEST.json](MANIFEST.json) records the original paths, source hashes,
archived hashes and link-only transformations. Research source files were left
unchanged. The original audit identifies commit
`451a36212f7ec0271094eed21e8ff09ddf8b3255`; later discussion addenda refine the
proposal without claiming another implementation.

At archiving, the local baseline state at
`/common/dev/gl638_integration_validation/publication/smoke/execution/state`
existed with manifest version 7. Its `global_settings.toml`,
`model_parameters.json` and GL638 graph hashes matched
[the tracked pilot provenance](../gl638/runtime-followup-pilot-provenance.json).
This is a read-only provenance check, not a fresh load or integration.
