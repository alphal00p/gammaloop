# GL638 all-orientation validation

This file records the reproducible all-orientation gate for the advanced
sampling work. It is deliberately a manual acceptance gate rather than a
normal CI test: unrestricted GL638 generation needs about 8.3 GiB peak RSS and
roughly six minutes on the recorded 16-thread host.

## What is covered

The gate must use the GL638 process card and DOT together, with an unset
orientation selector (the card's `run generate` block). This produces the
complete explicit sum over all 936 source orientations, all six physical cuts,
the 19 active threshold variants, direct three-dimensional local UV with
orientation localization, integrated UV, and the physical SM numerator. The
four `generate_o0` through `generate_o3` blocks in the card are representative
orientation diagnostics only; they are not a replacement for this gate.

The last successful source-stage run is archived in
`/common/dev/gl638_integration_validation/full-explicit/` (outside this
repository). Its `generation-summary.json` records:

- `returncode = 0`, `orientation_count = 936`, and
  `unique_orientation_count = 936`;
- binary SHA256 `ce79cac84b5dca80413441fa07fbbb4ab3d3ec60bdaba9dcb73624d6c48a4c5a`;
- elapsed time 335.414 s, state size 2.196 GB, and peak process-tree RSS
  8.316 GB;
- all 19 active variants and the complete component registry matching the
  reference state.

The archived hash values identify that run and must not be silently applied to
new source files. Check the current inputs before reusing the result:

```bash
sha256sum \
  examples/cli/epem_a_ttxh/NNLO/graphs/GL638.dot \
  examples/cli/epem_a_ttxh/NNLO/epem_a_tth_NNLO_test_GL638.toml
```

The tracked verifier performs this check and decodes the actual orientation
keys from the standalone archive; summary counts alone are rejected:

```bash
source /tmp/gammaloop-rebase-dev-env.sh
python3 docs/research/gl638/validate_all_orientations.py \
  --scope archive \
  --artifact /common/dev/gl638_integration_validation/full-explicit/generation-summary.json \
  --state /common/dev/gl638_integration_validation/full-explicit/eager/state \
  --card /common/dev/gl638_integration_validation/full-explicit/card.toml \
  --graph /common/dev/gl638_integration_validation/full/GL638.dot \
  --binary /common/dev/gl638_integration_validation/twenty-core/smooth-sliver/gammaloop-ce79cac84b5d \
  --standalone /common/dev/gl638_integration_validation/full-explicit/eager/standalone/processes/cross_sections/epem_a_tth/NNLO/standalone_cross_section.json \
  --metadata /common/dev/gl638_integration_validation/full-explicit/metadata-comparison.json
```

For a newly generated artifact, use `--scope current`. That mode requires a
recorded `source_commit` matching the current `HEAD` and a clean tracked
worktree, in addition to matching the card, DOT and binary SHA256 values. The
archived `generation-summary.json` predates this source-commit field, so it is
accepted only with the explicit `--scope archive` label. The verifier's
`--self-test` exercises both a successful synthetic archive and stale-DOT
rejection.

## Source-stage command

Build the release binary with the repository's normal Symbolica environment,
then run generation in a fresh state directory. Keep this run outside CI and
under the host's memory watchdog. The card imports the process itself; do not
add a `force_cuts` option or an orientation selector.

```bash
export GL638_BIN="$PWD/target/release/gammaloop"
export GL638_CARD="$PWD/examples/cli/epem_a_ttxh/NNLO/epem_a_tth_NNLO_test_GL638.toml"
export GL638_STATE="$(mktemp -d "${TMPDIR:-/tmp}/gl638-all-orientations.XXXXXX")"
"$GL638_BIN" --state-folder "$GL638_STATE" "$GL638_CARD" run generate
```

Record the binary hash, card hash, DOT hash, elapsed time, peak RSS, and the
orientation counts from the generated state. The generation result is not an
all-orientation *numerical* validation by itself: it only proves that every
source orientation and active metadata variant was constructed successfully.

## Numerical all-orientation checks

Load the generated state read-only and evaluate at least one fixed x-space
point with every physical cut retained. Compare the complete complex total and
each of the six cut contributions against an independently generated reference
state (or against a metadata-off control when testing a threshold change). Use
the same parent-frame coordinates and discrete channel ID in both evaluations.
Run normal and forced-Arb precision, inspect `is_nan`, the complete stability
history, and componentwise Re/Im errors. A successful process exit or a
`Stable` label alone is insufficient.

The portable point/replay protocol is in `replay-points.json` and `REPLAY.md`.
It includes the 936-orientation source count and all six cuts. The archived
full-explicit comparison currently covers one x-space point; it does **not**
constitute a multi-point convergence or bounded-weight proof. Exact coincident
A/P roots remain a known nonfinite case and must stay reported as a failure
until the common divided-difference/Hermite evaluation is implemented.

For the advanced sampler, the final gate additionally requires equal-budget
runs (independent seeds) with the ordinary baseline and each selected surface
catalogue. Report normalization, absolute-integral error, maximum weight,
finite/unstable/NaN counts, and per-cut cancellation. Keep the all-orientation
state and the representative-orientation diagnostics as separate artifacts.
