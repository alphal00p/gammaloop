{
  description = "Gammaloop";

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

      nixCiBarrierRevision =
        if self ? dirtyRev
        then self.dirtyRev
        else if self ? rev
        then self.rev
        else if self ? narHash
        then self.narHash
        else "local";
      # NixCI memoizes successful top-level derivations across commits without
      # re-realizing their closures. Salt only this zero-copy scheduling
      # wrapper so each commit primes the stable artifact in the shared cache.
      nixCiArtifactBarrier = name: artifact: let
        compilerState =
          if artifact ? incrementalBase
          then artifact.incrementalBase.incremental
          else artifact.incremental or null;
      in
        pkgs.runCommand "nix-ci-artifact-barrier-${name}" {
          NIX_CI_BARRIER_REVISION = nixCiBarrierRevision;
          passthru = { inherit artifact; };
        } (if compilerState != null then ''
          # Publish compiler state once per producer, keeping it out of the
          # dependency archives restored by other crates and runtime checks.
          mkdir -p "$out"
          for entry in ${artifact}/*; do
            ln -s "$entry" "$out/"
          done
          ln -s ${compilerState} "$out/incremental"
        '' else ''
          ln -s ${artifact} "$out"
        '');

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
