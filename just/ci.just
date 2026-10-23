# Print tracked source, or validate CI_FLAKE_SOURCE supplied by an enclosing recipe.
# Exporting requires Python 3 and Git/jj; reusing a store snapshot needs only Python.
nix-source:
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--print-source)
    if [[ -n "${CI_FLAKE_SOURCE:-}" ]]; then
        args+=(--flake-source "$CI_FLAKE_SOURCE")
    fi
    python3 scripts/tracked_source.py "${args[@]}"

# Build an attribute from tracked source, keeping credentials out of Nix.
[positional-arguments]
nix-build ATTRIBUTE *ARGS:
    #!/usr/bin/env bash
    set -euo pipefail
    attribute=$1; shift
    flake=$(just nix-source)
    env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix build "$flake#$attribute" "$@"

# Execute a named runtime package after its credential-free preparation.
# Runtime launcher packages have a matching bin/PACKAGE executable.
[positional-arguments]
nix-run PACKAGE *ARGS:
    #!/usr/bin/env bash
    set -euo pipefail
    flake=$(just nix-source)
    tmp=$(mktemp -d)
    trap 'rm -rf "$tmp"' EXIT
    package=$1; shift
    launcher=$(env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix build --out-link "$tmp/launcher" --print-out-paths "$flake#$package")
    "$launcher/bin/$package" "$@"

# Prepare the pinned development environment without sending runtime keys to Nix.
# With no command, open Bash. Otherwise run e.g. `just nix-develop just fmt --check`.
[positional-arguments]
nix-develop *COMMAND:
    #!/usr/bin/env bash
    set -euo pipefail
    flake=$(python3 scripts/tracked_source.py --print-source)
    tmp=$(mktemp -d)
    trap 'rm -rf "$tmp"' EXIT
    env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix print-dev-env "$flake" > "$tmp/environment"
    TMPDIR="$tmp" bash -e -c '
        # print-dev-env restores caller paths and executes shellHook itself.
        source "$1"
        shift
        if [ "$#" -eq 0 ]; then bash; else "$@"; fi
    ' bash "$tmp/environment" "$@"

# Check the flake from tracked source (e.g. `just nix-flake-check --no-build`).
[positional-arguments]
nix-flake-check *ARGS:
    flake=$(python3 scripts/tracked_source.py --print-source) && env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix flake check "$flake" "$@"

# Run workspace clippy via Nix (same as CI)
clippy-nix:
    just nix-build "checks.$(env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-clippy"

# Check formatting via Nix (same as CI)
fmt-check-nix:
    just nix-build "checks.$(env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-fmt"

# Build documentation via Nix (same as CI)
doc-nix:
    just nix-build "checks.$(env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-doc"

# Run doctests via Nix (same as CI)
doctest-nix:
    just nix-run nix-ci-check-gammaloop-doctest

# Run workspace nextest via Nix (same as CI)
test-nix-all:
    just nix-run nix-ci-check-gammaloop-nextest

# Run workspace nextest via Nix
test-nix:
    just test-nix-all

# Generate coverage report via Nix
coverage-nix:
    just nix-build "packages.$(env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-llvm-coverage"

# Run all CI checks locally (same as CI)
ci-checks:
    #!/usr/bin/env bash
    set -euo pipefail
    CI_FLAKE_SOURCE=$(python3 scripts/tracked_source.py --print-source)
    export CI_FLAKE_SOURCE
    just ci-runtime-check "$CI_FLAKE_SOURCE"
    just clippy-nix fmt-check-nix ci-graph-check doctest-nix test-nix docs-check-nix

# Validate runtime guards with public fixtures and dummy credentials only.
# Reuse only an explicitly supplied snapshot; standalone checks export afresh.
[positional-arguments]
ci-runtime-check SOURCE="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--print-source)
    if [[ -n "$1" ]]; then
        args+=(--flake-source "$1")
    fi
    flake=$(python3 scripts/tracked_source.py "${args[@]}")
    system=$(env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix eval --impure --raw --expr 'builtins.currentSystem')
    env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix build "$flake#checks.$system.gammaloop-ci-runtime"
    env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix run "$flake#ci-runtime-license-boundary" -- --flake-source "$flake" --report-dir target/ci-results/license-boundary

# Execute the complete documentation validation/render/publication suite.
docs-check-nix:
    just nix-run nix-ci-check-alphal00p-docs

# Check static inputs and upload compiled producers; leave runtime tests to NixCI.
# See https://nix-ci.com/documentation/nix-ci-cache
ci-checks-and-upload:
    #!/usr/bin/env bash
    set -euo pipefail
    netrc="${NIXCI_NETRC:-$HOME/.netrc}"
    test -r "$netrc" || { echo 'Set up ~/.netrc or NIXCI_NETRC before uploading.' >&2; exit 1; }

    tmp=$(mktemp -d)
    trap 'rm -rf "$tmp"' EXIT

    # Scheduled builds include static checks, runtime launchers and their producers.
    # These outputs prepare runtime jobs; every NixCI launcher executes tests.
    started=$SECONDS
    # Keep all selected closures rooted until the upload finishes, including
    # across the gap between realization and nix copy on a shared daemon.
    # Resolve the unchanged onlyBuild list through one flake evaluation. This
    # avoids evaluating a separate flake installable for every selected job.
    # Resolve the tracked source before Nix can import ignored workspace files.
    CI_FLAKE_SOURCE=$(python3 scripts/tracked_source.py --print-source)
    export CI_FLAKE_SOURCE
    # A failing credential probe must stop this pre-push path before upload.
    env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix run "$CI_FLAKE_SOURCE#ci-runtime-license-boundary" -- \
        --flake-source "$CI_FLAKE_SOURCE" --report-dir target/ci-results/license-boundary
    # --impure permits currentSystem and this public source reference only.
    env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix build --impure --out-link "$tmp/result" --print-out-paths --expr '
      let
        flake = builtins.getFlake (builtins.getEnv "CI_FLAKE_SOURCE");
        system = builtins.currentSystem;
        ci = import (flake.outPath + "/nix/ci.nix") {
          inherit system;
          workspaceGraph = builtins.fromJSON (builtins.readFile (flake.outPath + "/nix/ci-workspace-graph.json"));
        };
        lookup = path: flake.inputs.nixpkgs.lib.getAttrFromPath
          (builtins.filter builtins.isString (builtins.split "\\." path)) flake;
      in map lookup ci.configuration.onlyBuild
    ' > "$tmp/outputs"
    printf 'CI static checks and artifact preparation completed (%ss).\n' "$((SECONDS - started))"

    # Keep a throwaway XDG_CACHE_HOME so locally remembered unsigned narinfo cannot
    # prevent later substitution of the signed paths from the cache.
    # Also override NIX_CACHE_HOME, which takes precedence if inherited.
    export XDG_CACHE_HOME="$tmp/cache" NIX_CACHE_HOME="$tmp/cache/nix"

    # Run as the calling user with their netrc; a daemon post-build-hook would
    # require Nix trusted-user privileges to select the script and credentials.
    # Upload failures are visible failures of this explicitly requested phase;
    # prepared outputs remain cached and do not need to be rebuilt.
    started=$SECONDS
    env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix copy --netrc-file "$netrc" \
        --to 'https://cache.nix-ci.com?compression=xz&parallel-compression=true' \
        --stdin < "$tmp/outputs"
    printf 'CI upload completed (%ss).\n' "$((SECONDS - started))"

# Check the workspace graph and generated NixCI scheduling configuration.
ci-graph-check:
    just nix-build "checks.$(env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-guppy-workspace-graph"

# Regenerate both committed CI definitions after changing workspace manifests or groups.
ci-update:
    #!/usr/bin/env bash
    set -euo pipefail
    flake=$(python3 scripts/tracked_source.py --print-source)
    graph=$(env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix build --no-link --print-out-paths "$flake#ci-workspace-graph-json")
    cp "$graph" nix/ci-workspace-graph.json
    flake=$(python3 scripts/tracked_source.py --print-source)
    configuration=$(env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix build --no-link --print-out-paths "$flake#nix-ci-config")
    cp "$configuration" nix-ci.nix

# Pin a compatible application revision for per-crate compiler-state reuse.
# The current CI implementation builds the seed once; source edits reuse it.
ci-cache-base REVISION:
    #!/usr/bin/env bash
    set -euo pipefail
    revision={{quote(REVISION)}}
    [[ "$revision" =~ ^[0-9a-f]{40}$ ]] || { echo 'Use a full commit hash' >&2; exit 1; }
    sed -E "s@github:alphal00p/gammaloop/[0-9a-f]{40}@github:alphal00p/gammaloop/$revision@" flake.nix > flake.nix.tmp
    mv flake.nix.tmp flake.nix
    # Update only the tracked snapshot's lock, writing the result back explicitly.
    flake=$(python3 scripts/tracked_source.py --print-source)
    env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix flake update --flake "$flake" --output-lock-file "$PWD/flake.lock" ci-cache-base

# Collect paired NixCI and GitHub Actions measurements.
ci-report MANIFEST OUTPUT_DIR:
    flake=$(python3 scripts/tracked_source.py --print-source) && env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix run "$flake#ci-report" -- {{quote(MANIFEST)}} {{quote(OUTPUT_DIR)}}

# Prepare the native CI dependency artifacts in the Nix store.
build-deps-nix:
    just nix-build "packages.$(env -u SYMBOLICA_LICENSE -u SYMBOLICA_LICENSE_SIGNED nix eval --impure --raw --expr 'builtins.currentSystem').cargoArtifacts"
