{
  description = "Gammaloop";

  # Substitute what NixCI has already built instead of building it again.
  # Reading from this cache needs a token in your netrc as well, and Nix asks
  # before it trusts these settings; see CONTRIBUTING.md.
  nixConfig = {
    extra-substituters = ["https://cache.nix-ci.com"];
    extra-trusted-public-keys = ["nix-ci:g3xV5BDTLtIBZr/A00IU1x0EtKKlb7YLgBN2SgYgM6A="];
  };

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

    crane = {
      url = "github:ipetkov/crane";
    };

    # Refresh deliberately with `just ci-cache-base REVISION` after a green run.
    ci-cache-base = {
      url = "github:alphal00p/gammaloop/5181661ec340ebfb181a0045dac79fcec4f35525";
      flake = false;
    };

    fenix = {
      url = "github:nix-community/fenix";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.rust-analyzer-src.follows = "";
    };

    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = {
    self,
    nixpkgs,
    crane,
    fenix,
    flake-utils,
    ...
  }: let
    # Nixpkgs unstable no longer supports x86_64-darwin, so exposing that
    # system makes `nix flake show` fail before NixCI can list jobs.
    supportedSystems = [
      "x86_64-linux"
      "aarch64-linux"
      "aarch64-darwin"
    ];
  in
    flake-utils.lib.eachSystem supportedSystems (system: let
      pkgs = nixpkgs.legacyPackages.${system};
      inherit (pkgs) lib;

      # NixCI memoizes successful top-level derivations across commits without
      # re-realizing their closures. Keep publication stable so an unchanged
      # commit does not restore and upload the same compilation artifacts again.
      nixCiArtifactBarrier = name: artifact:
        pkgs.runCommand "nix-ci-artifact-barrier-${name}" {
          passthru = {inherit artifact;};
        } ''
          ln -s ${artifact} "$out"
        '';

      baseCraneLib = crane.mkLib pkgs;
      stableToolchain = fenix.packages.${system}.stable;

      ciToolchain = stableToolchain.withComponents [
        "cargo"
        "clippy"
        "llvm-tools"
        "rust-std"
        "rustc"
        "rustfmt"
      ];

      craneLib =
        baseCraneLib.overrideToolchain
        ciToolchain;

      wasmTarget = "wasm32-unknown-unknown";

      wasmToolchain = fenix.packages.${system}.combine [
        (stableToolchain.withComponents [
          "cargo"
          "rust-std"
          "rustc"
        ])
        fenix.packages.${system}.targets.${wasmTarget}.stable.rust-std
      ];

      wasmCraneLib =
        baseCraneLib.overrideToolchain
        wasmToolchain;

      workspace = import ./nix/rust-workspace.nix {
        inherit pkgs craneLib wasmCraneLib ciToolchain wasmTarget system nixCiArtifactBarrier;
        workspaceRoot = ./.;
        incrementalBaselineRoot = /. + builtins.unsafeDiscardStringContext self.inputs.ci-cache-base.outPath;
      };
      inherit
        (workspace)
        allChecks
        hestiaChecks
        gammaloop-cli
        clinnet-cli
        gammaloop-python-module
        nixCiConfiguration
        guppyWorkspaceGraphJson
        linnest-wasm
        linnestWasmCargoArtifacts
        cargoArtifacts
        cargoCheckArtifacts
        ciCompilerState
        gammaloopApiPackageArtifacts
        workspaceBuildArtifacts
        nixCiPassed
        cranePackageDependencyOutputs
        cranePackageOutputs
        craneTestDependencyOutputs
        craneTestBinaryPackageOutputs
        nextestContextualTestOutputs
        impureCheckRunnerPackages
        commonArgs
        workspaceTestSrc
        nixCc
        nixCxx
        cargoLinkerVar
        runtimeLibPath
        ;

      rscls = pkgs.rustPlatform.buildRustPackage rec {
        pname = "rscls";
        version = "0.2.3";
        src = pkgs.fetchCrate {
          inherit pname version;
          sha256 = "sha256-tahAhWCjhIVjbJ1NzrtiHBwGb/FBmUdK4XP9VlSPqh0=";
        };
        cargoHash = "sha256-JikjBTFeDh4XHBm57yiorsCwZhKikz0aiWNOTaMn0Vo=";
      };

      devShellPackages = with pkgs;
        [
          tdf
          cargo-flamegraph
          yaml-language-server
          just
          dot-language-server
          cargo-insta
          cargo-udeps
          cargo-machete
          openssl
          pyright
          gmp
          mpfr
          libmpc
          form
          gnum4
          nickel
          nls
          typst
          cargo-nextest
          pkg-config
          cargo-deny
          cargo-edit
          cargo-guppy
          cargo-hakari
          cargo-watch
          bacon
          jq
          gfortran
          gcc
          rust-script
          uv
          graphviz
          mupdf
          tinymist
          typstyle
          poppler-utils
          rust-analyzer
          maturin
          virtualenv
        ]
        ++ lib.optionals (!pkgs.stdenv.isDarwin) [
          valgrind
        ];

      # Nix runs this after every build in the dev shell, so NixCI can reuse
      # matching outputs from checks already built locally.
      # See https://nix-ci.com/documentation/nix-ci-cache
      pushToNixCiCache = pkgs.writeShellScript "push-to-nix-ci-cache" ''
        set -eu
        set -f
        export IFS=' '

        # Keep the throwaway XDG_CACHE_HOME. Without it this machine remembers
        # the unsigned narinfo it built locally, and then refuses to
        # substitute back the paths it pushed itself.
        # Also override NIX_CACHE_HOME, which takes precedence if inherited.
        XDG_CACHE_HOME="$(mktemp -d)"
        NIX_CACHE_HOME="$XDG_CACHE_HOME/nix"
        export XDG_CACHE_HOME NIX_CACHE_HOME
        trap 'rm -rf "$XDG_CACHE_HOME"' EXIT

        # Nix fails the build whose post-build-hook fails, so being offline or
        # without a token must not come out of here non-zero.
        nix copy --to 'https://cache.nix-ci.com?compression=xz&parallel-compression=true' ''${OUT_PATHS-} \
          || echo "push-to-nix-ci-cache: could not push to the NixCI cache." >&2
      '';

      mkDevShell = extraPackages:
        craneLib.devShell {
          # checks = self.checks.${system};

          RUST_SRC_PATH = "${pkgs.rustPlatform.rustLibSrc}";
          GLIBC_TUNABLES = "glibc.rtld.optional_static_tls=10000";

          CC = nixCc;
          CXX = nixCxx;
          "${cargoLinkerVar}" = nixCc;
          RUSTFLAGS = "-C linker=${nixCc}";

          LD_LIBRARY_PATH = runtimeLibPath;
          DYLD_LIBRARY_PATH = runtimeLibPath;

          # shellHook = ''
          #   export CC="${nixCc}"
          #   export CXX="${nixCxx}"
          #   export ${cargoLinkerVar}="${nixCc}"
          # '';

          # The hook runs as the Nix daemon user, which finds the cache token
          # through the netrc-file setting it inherits from here. Nix honours
          # both settings only for a trusted user; see CONTRIBUTING.md.
          shellHook = ''
            if [ -r "$HOME/.netrc" ]; then
              NIX_CONFIG="$(printf '%s\nnetrc-file = %s\npost-build-hook = %s' "''${NIX_CONFIG-}" "$HOME/.netrc" "${pushToNixCiCache}")"
              export NIX_CONFIG
            fi
          '';

          packages = devShellPackages ++ extraPackages;
        };
    in {
      checks = allChecks;

      hydraJobs = hestiaChecks;

      packages =
        {
          default = gammaloop-cli;
          clinnet = clinnet-cli;
          gammaloop = gammaloop-cli;
          inherit clinnet-cli;
          "gammaloop-python-module" = nixCiArtifactBarrier "gammaloop-python-module" gammaloop-python-module;
          "ci-workspace-graph-json" = guppyWorkspaceGraphJson;
          "nix-ci-config" = nixCiConfiguration;
          inherit linnest-wasm;
          linnestWasmCargoArtifacts =
            nixCiArtifactBarrier "linnest-wasm-cargo-artifacts" linnestWasmCargoArtifacts;
          "crane-ci-prebuild" = cargoArtifacts;
          cargoArtifacts = nixCiArtifactBarrier "cargo-artifacts" cargoArtifacts;
          cargoCheckArtifacts = nixCiArtifactBarrier "cargo-check-artifacts" cargoCheckArtifacts;
          ci-compiler-state = ciCompilerState;
          gammaloopApiPackageArtifacts =
            nixCiArtifactBarrier "gammaloop-api-package-artifacts" gammaloopApiPackageArtifacts;
          inherit workspaceBuildArtifacts;
          "nix-ci-passed" = nixCiPassed;
        }
        // cranePackageDependencyOutputs
        // cranePackageOutputs
        // craneTestDependencyOutputs
        // craneTestBinaryPackageOutputs
        // nextestContextualTestOutputs
        // impureCheckRunnerPackages
        // lib.optionalAttrs (!pkgs.stdenv.isDarwin) {
          gammaloop-llvm-coverage = craneLib.cargoLlvmCov (commonArgs
            // {
              src = workspaceTestSrc;
              inherit cargoArtifacts;
              nativeBuildInputs = (commonArgs.nativeBuildInputs or []) ++ [pkgs.form];
            });
        };

      apps = {
        ci-report = flake-utils.lib.mkApp {
          drv = pkgs.writeShellApplication {
            name = "ci-report";
            runtimeInputs = [pkgs.nodejs];
            text = ''exec node ${./.github/scripts/ci-report.mjs} "$@"'';
          };
        };
        default = flake-utils.lib.mkApp {
          drv = gammaloop-cli;
          exePath = "/bin/gammaloop";
        };
        gammaloop = flake-utils.lib.mkApp {
          drv = gammaloop-cli;
          exePath = "/bin/gammaloop";
        };
        clinnet = flake-utils.lib.mkApp {
          drv = clinnet-cli;
          exePath = "/bin/linnet";
        };
        linnet = flake-utils.lib.mkApp {
          drv = clinnet-cli;
          exePath = "/bin/linnet";
        };
      };

      devShells = {
        default = mkDevShell [clinnet-cli];
        full = mkDevShell [clinnet-cli rscls];
        clinnet = mkDevShell [clinnet-cli];
      };
    });
}
