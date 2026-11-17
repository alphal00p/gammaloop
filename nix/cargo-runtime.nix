{
  pkgs,
  craneLib,
}:
# Keep Crane's toolchain, vendoring and artifact hooks, but run the licensed
# build phase as the invoking user. The derivation only prepares its inputs.
args:
assert pkgs.lib.assertMsg (
  !(args ? SYMBOLICA_LICENSE || args ? SYMBOLICA_LICENSE_SIGNED)
  && !((args.env or {}) ? SYMBOLICA_LICENSE || (args.env or {}) ? SYMBOLICA_LICENSE_SIGNED)
) "Cargo runtime licenses must come from the invoking user's environment"; let
  inherit (pkgs) lib;
  name = args.pname;
  environment = craneLib.mkCargoDerivation (args
    // {
      phases = ["prepareRuntimePhase"];
      prepareRuntimePhase = ''mkdir -p "$out"'';
      doInstallCargoArtifacts = false;
      doNotLinkInheritedArtifacts = true;
      # stdenv's runPhase normally writes exported variables to env-vars.
      noDumpEnvVars = true;
      runtimePhase = ''
        set +x
        set -eo pipefail
        export out="$GAMMALOOP_RUNTIME_TEMP/result"
        export HOME="$GAMMALOOP_RUNTIME_TEMP/home"
        export XDG_CACHE_HOME="$GAMMALOOP_RUNTIME_TEMP/cache"
        export XDG_CONFIG_HOME="$GAMMALOOP_RUNTIME_TEMP/config"
        export XDG_DATA_HOME="$GAMMALOOP_RUNTIME_TEMP/data"
        export CARGO_HOME="$GAMMALOOP_RUNTIME_TEMP/cargo-home"
        export CARGO_NET_OFFLINE=true
        export TMPDIR="$GAMMALOOP_RUNTIME_TEMP/tmp"
        export TMP="$TMPDIR" TEMP="$TMPDIR" TEMPDIR="$TMPDIR"
        export NIX_BUILD_TOP="$GAMMALOOP_RUNTIME_TEMP/build"
        mkdir -p "$out" "$HOME" "$XDG_CACHE_HOME" "$XDG_CONFIG_HOME" \
          "$XDG_DATA_HOME" "$CARGO_HOME" "$TMPDIR" "$NIX_BUILD_TOP"
        cd "$NIX_BUILD_TOP"
        runPhase unpackPhase
        # runPhase changes into sourceRoot after unpacking (unpackPhase alone
        # does not). Export the staged root before any patch or Cargo hooks.
        export GAMMALOOP_WORKSPACE_ROOT="$PWD"
        export NEXTEST_WORKSPACE_ROOT="$PWD"
        export INSTA_WORKSPACE_ROOT="$PWD"
        export CARGO_TARGET_DIR=${lib.escapeShellArg (args.CARGO_TARGET_DIR or "target")}
        runPhase patchPhase
        runPhase configurePhase
        runPhase buildPhase
      '';
    });
in
  pkgs.writeShellApplication {
    inherit name;
    # nix develop resolves bash through PATH when no registry shell is available.
    runtimeInputs = [pkgs.bash pkgs.nix pkgs.coreutils pkgs.findutils];
    text = ''
      set +x
      if [[ $# -ne 1 || -z "$1" ]]; then
        echo "Usage: ${name} OUTPUT_DIRECTORY" >&2
        exit 2
      fi
      destination=$(realpath -m -- "$1")
      if [[ -e "$destination" ]] && {
        [[ ! -d "$destination" ]] ||
          [[ -n "$(find "$destination" -mindepth 1 -maxdepth 1 -print -quit)" ]];
      }; then
        echo "Output directory must be absent or empty: $destination" >&2
        exit 2
      fi
      if [[ -z "''${SYMBOLICA_LICENSE:-}" ]]; then
        export SYMBOLICA_LICENSE="''${SYMBOLICA_LICENSE_SIGNED:-}"
      fi
      if [[ -z "''${SYMBOLICA_LICENSE:-}" ]]; then
        license_file="''${XDG_CONFIG_HOME:-''${HOME:-}/.config}/symbolica/license"
        if [[ -r "$license_file" ]]; then
          SYMBOLICA_LICENSE=$(< "$license_file")
          export SYMBOLICA_LICENSE
        fi
      fi
      if [[ -z "''${SYMBOLICA_LICENSE:-}" ]]; then
        echo "Missing SYMBOLICA_LICENSE; export it or configure symbolica/license" >&2
        exit 1
      fi
      umask 077
      # The derivation restores any Python configuration it explicitly needs.
      # Ambient virtual environments and Cargo paths belong to the caller.
      # Leave Cargo's build-dir default tied to each command's staged target dir.
      unset PYTHONHOME PYTHONPATH VIRTUAL_ENV CONDA_PREFIX
      unset CARGO_TARGET_DIR CARGO_BUILD_BUILD_DIR
      GAMMALOOP_RUNTIME_TEMP=$(mktemp -d -t ${name}.XXXXXXXX)
      export GAMMALOOP_RUNTIME_TEMP
      trap 'rm -rf -- "$GAMMALOOP_RUNTIME_TEMP"' EXIT
      trap 'exit 130' INT
      trap 'exit 143' TERM
      # Retain the exact build inputs without realizing every output of every
      # transitive compiler derivation. This is a secret-free Nix environment;
      # the user's license only enters the shell launched by nix develop.
      # Runtime dependency closure: ${environment.inputDerivation}
      nix --extra-experimental-features 'nix-command flakes' \
        develop ${builtins.unsafeDiscardOutputDependency environment.drvPath} --phase runtime
      mkdir -p -- "$destination"
      cp -R "$GAMMALOOP_RUNTIME_TEMP/result"/. "$destination"/
    '';
  }
