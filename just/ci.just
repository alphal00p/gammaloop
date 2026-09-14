# Run workspace clippy via Nix (same as CI)
clippy-nix:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-clippy

# Check formatting via Nix (same as CI)
fmt-check-nix:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-fmt

# Run security audit via Nix (same as CI)
audit-nix:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-audit

# Run license check via Nix (same as CI)
deny-nix:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-deny

# Build documentation via Nix (same as CI)
doc-nix:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-doc

# Run doctests via Nix (same as CI)
doctest-nix:
    nix build --impure .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-doctest

# Run workspace nextest via Nix (same as CI)
test-nix-all:
    nix build --impure .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-nextest

# Run workspace nextest via Nix
test-nix:
    just test-nix-all

# Generate coverage report via Nix
coverage-nix:
    nix build .#packages.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-llvm-coverage

# Run all CI checks locally (same as CI)
ci-checks: clippy-nix fmt-check-nix ci-graph-check doctest-nix test-nix

# Build selected CI outputs explicitly so NixCI can reuse matching local results.
# See https://nix-ci.com/documentation/nix-ci-cache
ci-checks-and-upload:
    #!/usr/bin/env bash
    set -euo pipefail
    netrc="${NIXCI_NETRC:-$HOME/.netrc}"
    test -r "$netrc" || { echo 'Set up ~/.netrc or NIXCI_NETRC before uploading.' >&2; exit 1; }

    started=$SECONDS
    just ci-checks
    printf 'CI checks passed (%ss).\n' "$((SECONDS - started))"

    selection=$(nix eval --impure --raw --expr '
      let
        system = builtins.currentSystem;
        ci = import ./nix/ci.nix {
          inherit system;
          workspaceGraph = builtins.fromJSON (builtins.readFile ./nix/ci-workspace-graph.json);
        };
        testResults = map (name: "checks.${system}.${name}") (builtins.attrNames ci.configuration.test);
      in builtins.concatStringsSep "\n" (ci.configuration.onlyBuild ++ testResults)
    ')
    targets=()
    while IFS= read -r target; do targets+=(".#$target"); done <<< "$selection"
    tmp=$(mktemp -d)
    trap 'rm -rf "$tmp"' EXIT

    # Test-result outputs alone do not retain the binaries and compiler artifacts.
    # Realize the scheduled producers too, including their stable publication links.
    started=$SECONDS
    nix build --impure --no-link --print-out-paths "${targets[@]}" > "$tmp/outputs"
    printf 'CI publication preparation completed (%ss).\n' "$((SECONDS - started))"

    # Keep a throwaway XDG_CACHE_HOME so locally remembered unsigned narinfo cannot
    # prevent later substitution of the signed paths from the cache.
    # Also override NIX_CACHE_HOME, which takes precedence if inherited.
    export XDG_CACHE_HOME="$tmp/cache" NIX_CACHE_HOME="$tmp/cache/nix"

    # Run as the calling user with their netrc; a daemon post-build-hook would
    # require Nix trusted-user privileges to select the script and credentials.
    # Upload failures are visible failures of this explicitly requested phase;
    # successful checks remain cached and do not need to be rebuilt.
    started=$SECONDS
    nix copy --netrc-file "$netrc" \
        --to 'https://cache.nix-ci.com?compression=xz&parallel-compression=true' \
        --stdin < "$tmp/outputs"
    printf 'CI upload completed (%ss).\n' "$((SECONDS - started))"

# Check the workspace graph and generated NixCI scheduling configuration.
ci-graph-check:
    nix build .#checks.$(nix eval --impure --raw --expr 'builtins.currentSystem').gammaloop-guppy-workspace-graph

# Regenerate both committed CI definitions after changing workspace manifests or groups.
ci-update:
    #!/usr/bin/env bash
    set -euo pipefail
    graph=$(nix build --no-link --print-out-paths .#ci-workspace-graph-json)
    cp "$graph" nix/ci-workspace-graph.json
    configuration=$(nix build --no-link --print-out-paths .#nix-ci-config)
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
    nix flake update ci-cache-base

# Collect paired NixCI and GitHub Actions measurements.
ci-report MANIFEST OUTPUT_DIR:
    nix run .#ci-report -- {{quote(MANIFEST)}} {{quote(OUTPUT_DIR)}}

# Build cargo dependencies via Nix (useful for caching) - uses release mode
build-deps-nix:
    nix build .#packages.$(nix eval --impure --raw --expr 'builtins.currentSystem').cargoArtifacts
