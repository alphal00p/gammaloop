# Tangled spindle CI

This repository's Tangled spindle port lives in `.tangled/workflows/`.

The mapping is intentionally conservative:

- `nix-ci.yml` runs the Linux Nix checks from `.github/workflows/crane-workflow.yml`.
- `nix-package.yml` keeps the package build separate, matching the GitHub workflow's non-blocking package job.
- The Linux CI follows Crane's workspace-oriented pattern: one shared dependency cache, then workspace `gammaloop-fmt`, workspace `gammaloop-clippy`, workspace `gammaloop-nextest`, and workspace `gammaloop-doctest`.
- Package builds stay separate from CI gating, matching the usual Crane examples where consumer-facing packages are distinct from lint/test checks.

The spindle workflow runs the Linux CI in a single step instead of splitting it into many steps. Tangled spindle executes each step in a fresh container and only preserves `/tangled/workspace`, so keeping the Nix builds in one step avoids discarding the warmed `/nix/store` between steps. The Nix commands are inlined directly in the workflow rather than delegated through a helper script.

Because Tangled's Nix builds run without the usual sandbox isolation, the workflows also remove `/homeless-shelter` before invoking `nix build`; otherwise stdenv aborts immediately on the fallback `HOME` purity check.

Not ported:

- `build-deps-macos`
- `test-macos`
- GitHub-specific cache plumbing (`cache-nix-action`, Cachix setup)
- GitHub-specific reporting and uploads (`GITHUB_STEP_SUMMARY`, Codecov action)

Required spindle repository secret:

- `SYMBOLICA_LICENSE`
