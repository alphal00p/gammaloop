# Tangled spindle CI

This repository's Tangled spindle port lives in `.tangled/workflows/`.

The mapping is intentionally conservative:

- `nix-ci.yml` runs the Linux Nix checks from `.github/workflows/crane-workflow.yml`.
- `nix-package.yml` keeps the package build separate, matching the GitHub workflow's non-blocking package job.
- The Linux lint phase keeps a single workspace `gammaloop-fmt` check, then runs package-local `clippy-*` checks in dependency order.
- The Linux test phase runs package-local `nextest-*` checks in dependency order and finishes with the integration-only `nextest-integration` check, avoiding a repeated full-workspace test pass.

The spindle workflow runs the Linux CI in a single step instead of splitting it into many steps. Tangled spindle executes each step in a fresh container and only preserves `/tangled/workspace`, so keeping the Nix builds in one step avoids discarding the warmed `/nix/store` between steps. The Nix commands are inlined directly in the workflow rather than delegated through a helper script.

Not ported:

- `build-deps-macos`
- `test-macos`
- GitHub-specific cache plumbing (`cache-nix-action`, Cachix setup)
- GitHub-specific reporting and uploads (`GITHUB_STEP_SUMMARY`, Codecov action)

Required spindle repository secret:

- `SYMBOLICA_LICENSE`
