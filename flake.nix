{
  description = "Gammaloop";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

    crane = {
      url = "github:ipetkov/crane";
    };

    fenix = {
      url = "github:nix-community/fenix";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.rust-analyzer-src.follows = "";
    };

    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = {
    nixpkgs,
    crane,
    fenix,
    flake-utils,
    ...
  }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = nixpkgs.legacyPackages.${system};
      inherit (pkgs) lib;

      craneLib =
        (crane.mkLib pkgs).overrideToolchain
        fenix.packages.${system}.stable.toolchain;

      craneLibLLvmTools =
        craneLib.overrideToolchain
        (fenix.packages.${system}.stable.withComponents [
          "cargo"
          "llvm-tools"
          "rustc"
        ]);

      workspaceRoot = ./.;

      cargoSources = craneLib.fileset.commonCargoSources workspaceRoot;

      nonCargoBuildSources = lib.fileset.unions [
        ./.config
        ./assets
        ./crates/clinnet/templates
        ./crates/vakint/form_src
        ./crates/vakint/templates
      ];

      snapshotSources = lib.fileset.unions [
        ./crates/gammalooprs
        ./crates/linnet
        ./crates/spenso
      ];

      workspaceBuildSrc = lib.fileset.toSource {
        root = workspaceRoot;
        fileset = lib.fileset.unions [
          cargoSources
          nonCargoBuildSources
        ];
      };

      workspaceTestSrc = lib.fileset.toSource {
        root = workspaceRoot;
        fileset = lib.fileset.unions [
          cargoSources
          nonCargoBuildSources
          snapshotSources
          ./tests
          ./examples/cli
        ];
      };

      src = workspaceBuildSrc;

      apiMeta = craneLib.crateNameFromCargoToml {
        cargoToml = ./crates/gammaloop-api/Cargo.toml;
      };

      # Host Rust target triple, e.g. x86_64-unknown-linux-gnu
      rustTarget = pkgs.stdenv.hostPlatform.rust.rustcTarget;

      # Env var name Cargo uses to pick the linker for this target
      cargoLinkerVar = "CARGO_TARGET_${lib.toUpper (lib.replaceStrings ["-"] ["_"] rustTarget)}_LINKER";

      # Force GCC as both C/C++ compiler and Rust linker.
      nixCc = "${pkgs.gcc}/bin/gcc";
      nixCxx = "${pkgs.gcc}/bin/g++";

      # Runtime library search path for locally-built binaries and for maturin/auditwheel
      runtimeLibPath = lib.makeLibraryPath [
        pkgs.python313
        pkgs.gmp
        pkgs.mpfr
        pkgs.libmpc
        pkgs.openssl
        pkgs.stdenv.cc.cc.lib
      ];

      # Common arguments can be set here to avoid repeating them later
      commonArgs = {
        inherit src;
        pname = "gammaloop-workspace";
        inherit (apiMeta) version;
        strictDeps = true;

        nativeBuildInputs =
          [
            pkgs.pkg-config
            pkgs.gcc
            pkgs.git
            pkgs.python313
            pkgs.gnum4
          ]
          ++ lib.optionals pkgs.stdenv.isDarwin [
            pkgs.darwin.cctools
          ];

        buildInputs =
          [
            pkgs.openssl
            pkgs.gmp
            pkgs.gmp.dev
            pkgs.mpfr
            pkgs.mpfr.dev
            pkgs.libmpc
            pkgs.python313
          ]
          ++ lib.optionals pkgs.stdenv.isDarwin [
            pkgs.libiconv
          ];

        CC = nixCc;
        CXX = nixCxx;
        "${cargoLinkerVar}" = nixCc;
        RUSTFLAGS = "-C linker=${nixCc}";

        LD_LIBRARY_PATH = runtimeLibPath;
        DYLD_LIBRARY_PATH = runtimeLibPath;
      };

      ciArgs =
        commonArgs
        // {
          buildType = "dev-optim";
          cargoExtraArgs = "--locked";

          PYO3_PYTHON = "${pkgs.python313}/bin/python3";
          PYTHONPATH = "${pkgs.python313}/lib/python3.13/site-packages";
        };

      ciTestCargoProfile = "dev-optim";

      workspaceCrates = [
        {
          attr = "clinnet";
          package = "clinnet";
          path = ./crates/clinnet;
          usesSymbolica = false;
        }
        {
          attr = "linnet";
          package = "linnet";
          path = ./crates/linnet;
          usesSymbolica = false;
        }
        {
          attr = "spenso-macros";
          package = "spenso-macros";
          path = ./crates/spenso-macros;
          usesSymbolica = true;
        }
        {
          attr = "vakint";
          package = "vakint";
          path = ./crates/vakint;
          usesSymbolica = true;
        }
        {
          attr = "linnest";
          package = "linnest";
          path = ./crates/linnest;
          usesSymbolica = false;
        }
        {
          attr = "linnet-py";
          package = "linnet-py";
          path = ./crates/linnet-py;
          usesSymbolica = false;
        }
        {
          attr = "spenso";
          package = "spenso";
          path = ./crates/spenso;
          usesSymbolica = true;
        }
        {
          attr = "idenso";
          package = "idenso";
          path = ./crates/idenso;
          usesSymbolica = true;
        }
        {
          attr = "spenso-hep-lib";
          package = "spenso-hep-lib";
          path = ./crates/spenso-hep-lib;
          usesSymbolica = true;
        }
        {
          attr = "spynso3";
          package = "spynso3";
          path = ./crates/spynso3;
          usesSymbolica = true;
        }
        {
          attr = "gammalooprs";
          package = "gammalooprs";
          path = ./crates/gammalooprs;
          usesSymbolica = true;
        }
        {
          attr = "gammaloop-api";
          package = "gammaloop-api";
          path = ./crates/gammaloop-api;
          usesSymbolica = true;
        }
      ];

      licensePreCheck = ''
        if [ -z "''${SYMBOLICA_LICENSE:-}" ]; then
          echo "Missing SYMBOLICA_LICENSE environment variable" >&2
          exit 1
        fi
      '';

      # Package builds still need the full workspace crate tree because several crates embed
      # non-Rust assets (templates, FORM sources, wasm payloads) at compile time, but they do
      # not need integration tests in their source closure.
      fileSetForCrate = _: workspaceBuildSrc;

      cargoExtraArgsForPackages = packages:
        lib.concatStringsSep " " (["--locked"] ++ map (package: "-p ${package}") packages);

      mkDummySrc = source:
        craneLib.mkDummySrc {
          src = source;
          cleanCargoTomlFilter = craneLib.filters.cargoTomlAggressive;
        };

      mkSubsetArtifacts = {
        pname,
        packages,
        cargoArtifacts ? null,
        src ? workspaceBuildSrc,
        extraArgs ? {},
        cargoBuildCommand ? null,
        buildPhaseCargoCommand ? null,
      }:
        craneLib.buildDepsOnly ((ciArgs // extraArgs)
          // {
            inherit pname cargoArtifacts;
            inherit (apiMeta) version;
            dummySrc = mkDummySrc src;
            cargoExtraArgs = cargoExtraArgsForPackages packages;
            doCheck = false;
          }
          // lib.optionalAttrs (cargoBuildCommand != null) {
            inherit cargoBuildCommand;
          }
          // lib.optionalAttrs (buildPhaseCargoCommand != null) {
            inherit buildPhaseCargoCommand;
          });

      nextestDerivationArgs = {
        CARGO_PROFILE = ciTestCargoProfile;
        INSTA_WORKSPACE_ROOT = ".";
        nativeBuildInputs = (ciArgs.nativeBuildInputs or []) ++ [pkgs.form];
      };

      nextestExtraArgs = "--profile ci_gammaloop --no-fail-fast --final-status-level fail --no-tests=pass --run-ignored all";

      mkNextestArtifacts = {
        pname,
        packages,
        cargoArtifacts ? null,
        src ? workspaceTestSrc,
      }:
        mkSubsetArtifacts {
          inherit pname packages cargoArtifacts src;
          extraArgs = nextestDerivationArgs;
          cargoBuildCommand = "cargo test --profile ${ciTestCargoProfile} --no-run ${cargoExtraArgsForPackages packages}";
          buildPhaseCargoCommand = "cargo test --profile ${ciTestCargoProfile} --no-run ${cargoExtraArgsForPackages packages}";
        };

      stagedCrateArtifactsState = lib.lists.foldl' (
        state: crate: let
          artifacts = mkSubsetArtifacts {
            pname = "gammaloop-${crate.attr}-artifacts";
            packages = [crate.package];
            cargoArtifacts = state.previousArtifacts;
          };
        in {
          previousArtifacts = artifacts;
          artifactsByAttr = state.artifactsByAttr // {"${crate.attr}" = artifacts;};
        }
      ) {
        previousArtifacts = null;
        artifactsByAttr = {};
      } workspaceCrates;

      crateCargoArtifacts = stagedCrateArtifactsState.artifactsByAttr;

      # Top-of-tree artifacts cover the full workspace crate graph and are reused for the
      # higher-level clippy/doc/package/coverage builds.
      cargoArtifacts = stagedCrateArtifactsState.previousArtifacts;

      symbolicaCrateArgs = usesSymbolica:
        lib.optionalAttrs usesSymbolica {
          preBuild = licensePreCheck;
          SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
        };

      gammaloop-cli = craneLib.buildPackage (commonArgs
        // {
          cargoArtifacts = builtins.getAttr "gammaloop-api" crateCargoArtifacts;
          buildType = "dev-optim";
          doCheck = false;
          pname = "gammaloop";
          inherit (apiMeta) version;
          src = fileSetForCrate ./crates/gammaloop-api;
          cargoExtraArgs = "--locked -p gammaloop-api --bin gammaloop";
        });

      stagedNextestArtifactsState = lib.lists.foldl' (
        state: crate: let
          artifacts = mkNextestArtifacts {
            pname = "gammaloop-nextest-${crate.attr}-artifacts";
            packages = [crate.package];
            cargoArtifacts = state.previousArtifacts;
          };
        in {
          previousArtifacts = artifacts;
          artifactsByAttr = state.artifactsByAttr // {"${crate.attr}" = artifacts;};
        }
      ) {
        previousArtifacts = cargoArtifacts;
        artifactsByAttr = {};
      } workspaceCrates;

      nextestCargoArtifacts = stagedNextestArtifactsState.artifactsByAttr;

      integrationNextestCargoArtifacts = mkNextestArtifacts {
        pname = "gammaloop-nextest-integration-artifacts";
        packages = ["gammaloop-integration-tests"];
        cargoArtifacts = stagedNextestArtifactsState.previousArtifacts;
      };

      symbolicaNextestArgs = usesSymbolica:
        lib.optionalAttrs usesSymbolica {
          preCheck = licensePreCheck;
          SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
        };

      packageClippyChecks = lib.listToAttrs (map (crate: {
          name = "clippy-${crate.attr}";
          value = craneLib.cargoClippy ((ciArgs // symbolicaCrateArgs crate.usesSymbolica)
            // {
              cargoArtifacts = builtins.getAttr crate.attr crateCargoArtifacts;
              src = workspaceBuildSrc;
              cargoExtraArgs = "--locked -p ${crate.package}";
              cargoClippyExtraArgs = "--all-targets -- --deny warnings";
            });
        })
        workspaceCrates);

      nextestCheckForPackages = {
        checkName,
        packages,
        cargoArtifacts,
        usesSymbolica ? false,
      }:
        {
          name = checkName;
          value = craneLib.cargoNextest ((ciArgs // nextestDerivationArgs)
            // {
              inherit cargoArtifacts;
              src = workspaceTestSrc;
              cargoExtraArgs = cargoExtraArgsForPackages packages;
              cargoNextestExtraArgs = nextestExtraArgs;
            }
            // symbolicaNextestArgs usesSymbolica);
        };

      packageNextestChecks = lib.listToAttrs (map (crate:
          nextestCheckForPackages {
            checkName = "nextest-${crate.attr}";
            packages = [crate.package];
            cargoArtifacts = builtins.getAttr crate.attr nextestCargoArtifacts;
            usesSymbolica = crate.usesSymbolica;
          })
        workspaceCrates);

      integrationNextestCheck = nextestCheckForPackages {
        checkName = "nextest-integration";
        packages = ["gammaloop-integration-tests"];
        cargoArtifacts = integrationNextestCargoArtifacts;
        usesSymbolica = true;
      };

      impureCheckRunnerTargets =
        (map (crate: {
            runnerAttr = "nix-ci-check-clippy-${crate.attr}";
            checkAttr = "clippy-${crate.attr}";
          })
          (builtins.filter (crate: crate.usesSymbolica) workspaceCrates))
        ++
        (map (crate: {
            runnerAttr = "nix-ci-check-nextest-${crate.attr}";
            checkAttr = "nextest-${crate.attr}";
          })
          (builtins.filter (crate: crate.usesSymbolica) workspaceCrates))
        ++ [
          {
            runnerAttr = "nix-ci-check-nextest-integration";
            checkAttr = "nextest-integration";
          }
        ];

      impureCheckRunnerPackages = lib.listToAttrs (map (target: {
          name = target.runnerAttr;
          value = pkgs.writeShellApplication {
            name = target.runnerAttr;
            runtimeInputs = [pkgs.nix];
            text = ''
              set -euo pipefail
              exec nix \
                --extra-experimental-features nix-command \
                --extra-experimental-features flakes \
                build \
                --no-link \
                --print-build-logs \
                --impure \
                .#checks.${system}.${target.checkAttr}
            '';
          };
        })
        impureCheckRunnerTargets);
    in {
      checks =
        {
          # Keep existing check names for CI compatibility.
          gammaloop = gammaloop-cli;

          gammaloop-clippy = craneLib.cargoClippy (ciArgs
            // {
              inherit cargoArtifacts;
              cargoClippyExtraArgs = "--all-targets -- --deny warnings";
            });

          gammaloop-doc = craneLib.cargoDoc (ciArgs
            // {
              inherit cargoArtifacts;
            });

          gammaloop-fmt = craneLib.cargoFmt {
            inherit src;
            pname = "gammaloop-workspace";
            inherit (apiMeta) version;
          };

          ${integrationNextestCheck.name} = integrationNextestCheck.value;
        }
        // packageClippyChecks
        // packageNextestChecks;

      packages =
        {
          default = gammaloop-cli;
          gammaloop = gammaloop-cli;
          inherit cargoArtifacts;
        }
        // impureCheckRunnerPackages
        // lib.optionalAttrs (!pkgs.stdenv.isDarwin) {
          gammaloop-llvm-coverage = craneLibLLvmTools.cargoLlvmCov (commonArgs
            // {
              src = workspaceTestSrc;
              inherit cargoArtifacts;
              nativeBuildInputs = (commonArgs.nativeBuildInputs or []) ++ [pkgs.form];
            });
        };

      apps = {
        default = flake-utils.lib.mkApp {
          drv = gammaloop-cli;
        };
        gammaloop = flake-utils.lib.mkApp {
          drv = gammaloop-cli;
        };
      };

      devShells.default = craneLib.devShell {
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

        packages = with pkgs; [
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
          cargo-watch
          bacon
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
          (pkgs.rustPlatform.buildRustPackage rec {
            pname = "clinnet";
            version = "0.1.8";
            src = pkgs.fetchCrate {
              inherit pname version;
              sha256 = "sha256-CbZBHbf+8bIkdiSI5LMFO2Qc3zDr9UEBEry+fZOuep8=";
            };
            cargoHash = "sha256-GTixU2ZJZVMrEWLOfWjEnXMVLG2+cpkPbJuNnkTuFfo=";
          })
          (pkgs.rustPlatform.buildRustPackage rec {
            pname = "rscls";
            version = "0.2.3";
            src = pkgs.fetchCrate {
              inherit pname version;
              sha256 = "sha256-tahAhWCjhIVjbJ1NzrtiHBwGb/FBmUdK4XP9VlSPqh0=";
            };
            cargoHash = "sha256-JikjBTFeDh4XHBm57yiorsCwZhKikz0aiWNOTaMn0Vo=";
          })
        ];
      };
    });
}
