# Advanced sampling research archive

This is a portable snapshot of the research supporting
[ADVANCED_SAMPLING_PLAN.md](../../../ADVANCED_SAMPLING_PLAN.md). These are research
proposals and recorded arithmetic checks, not an implemented sampling system or
new GammaLoop integration results.

## Reading order and precedence

Agents implementing or auditing this work should consult the plan and relevant
research again at each implementation slice, before changing an invariant, and
when investigating a mismatch. Do not infer current requirements from an older
example in isolation. Where proposal details differ, use this precedence:

1. [ADVANCED_SAMPLING_PLAN.md](../../../ADVANCED_SAMPLING_PLAN.md): consolidated
   requirements, milestones and acceptance criteria.
2. [SOPER_AND_SAMPLING_API.md](SOPER_AND_SAMPLING_API.md): latest agreed string
   language, automatic presets, conditional existence, eager/dual derivatives,
   exact/proxy partitions and saved-state acceptance harness.
3. [API_AND_GENERICITY_ADDENDUM.md](API_AND_GENERICITY_ADDENDUM.md): genericity,
   amplitudes and physical-cut/left/right composition.
4. [REPORT.md](REPORT.md): original mathematical study, source audit and GL638
   derivations. Its original API syntax and staged scope are historical where
   the documents above supersede them.

The focused agent studies preserve independent reasoning and review:

- [api.md](api.md): existing code owners and estimator boundaries; API sketches
  are historical.
- [gl638.md](gl638.md): GL638 geometry, projected targets and numerical scope.
- [literature.md](literature.md): primary literature and mathematical derivations.
- [SOPER_ANGULAR_REVIEW.md](SOPER_ANGULAR_REVIEW.md): angular and correlated density
  construction.
- [TWO_NORMAL_PROPOSAL.md](TWO_NORMAL_PROPOSAL.md): exact H/Z chart and A-star
  pullback derivation, with explicit regularity and branch conditions.

## Numerical evidence and portability

The small records included here are:

| File | Meaning |
|---|---|
| [geometry_checks.json](geometry_checks.json) | Ellipsoid forward/inverse, independent Jacobian, volume and moment checks |
| [gl638_geometry_checks.json](gl638_geometry_checks.json) | Recorded GL638 rank, LU radial invariance and A-star affine arithmetic |
| [normal_plane_moments.csv](normal_plane_moments.csv) | Analytically normalized model comparisons for a codimension-two singularity |
| [hard_hz_points.json](hard_hz_points.json) | Frozen inputs and derived geometry for the H/Z approach |
| [a_star_pullback_points.json](a_star_pullback_points.json) | Frozen actual center, complements and A-star points |

The frozen inputs and equations preserve the numerical setup without requiring
the original machine. They are not executable acceptance tests. The original
`gl638_geometry_checks.py` is deliberately excluded because it imports local
scratch helpers and reads machine-specific files. Its output is retained as
historical evidence; implementation must add independently runnable tests.

Existing physical integration evidence is already tracked in
[the GL638 research records](../gl638/README.md), especially
[the runtime follow-up](../gl638/runtime-followup.md) and its portable replay
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
