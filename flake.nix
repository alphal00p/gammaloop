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

      workspaceSrc = lib.fileset.toSource {
        root = ./.;
        fileset = lib.fileset.unions [
          ./.config
          ./Cargo.toml
          ./Cargo.lock
          ./assets
          ./crates
          ./tests
        ];
      };

      workspaceBuildSrc = lib.fileset.toSource {
        root = ./.;
        fileset = lib.fileset.unions [
          ./.config
          ./Cargo.toml
          ./Cargo.lock
          ./assets
          ./crates
        ];
      };

      src = workspaceSrc;

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

      gammaloopNextestPackages = [
        "gammaloop-api"
        "gammalooprs"
        "idenso"
        "linnest"
        "linnet"
        "spenso"
        "spenso-hep-lib"
        "spenso-macros"
        "vakint"
        "gammaloop-integration-tests"
      ];

      ciTestCargoProfile = "dev-optim";

      gammaloopNextestArgs =
        ciArgs
        // {
          CARGO_PROFILE = ciTestCargoProfile;
          cargoExtraArgs = lib.concatStringsSep " " (["--locked"] ++ map (package: "-p ${package}") gammaloopNextestPackages);
        };

      ciPartitionCount = 6;

      workspaceCrates = [
        {
          attr = "clinnet";
          package = "clinnet";
          path = ./crates/clinnet;
          usesSymbolica = false;
        }
        {
          attr = "gammaloop-api";
          package = "gammaloop-api";
          path = ./crates/gammaloop-api;
          usesSymbolica = true;
        }
        {
          attr = "gammalooprs";
          package = "gammalooprs";
          path = ./crates/gammalooprs;
          usesSymbolica = true;
        }
        {
          attr = "idenso";
          package = "idenso";
          path = ./crates/idenso;
          usesSymbolica = true;
        }
        {
          attr = "linnest";
          package = "linnest";
          path = ./crates/linnest;
          usesSymbolica = false;
        }
        {
          attr = "linnet";
          package = "linnet";
          path = ./crates/linnet;
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
          attr = "spenso-hep-lib";
          package = "spenso-hep-lib";
          path = ./crates/spenso-hep-lib;
          usesSymbolica = true;
        }
        {
          attr = "spenso-macros";
          package = "spenso-macros";
          path = ./crates/spenso-macros;
          usesSymbolica = true;
        }
        {
          attr = "spynso3";
          package = "spynso3";
          path = ./crates/spynso3;
          usesSymbolica = true;
        }
        {
          attr = "vakint";
          package = "vakint";
          path = ./crates/vakint;
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

      # Build workspace dependency artifacts once and reuse for downstream checks/packages.
      cargoArtifacts = craneLib.buildDepsOnly (ciArgs
        // {
          pname = "gammaloop-workspace-artifacts";
          inherit (apiMeta) version;
          cargoBuildCommand = "cargo test --workspace --all-targets --profile ${ciTestCargoProfile} --no-run --locked";
          buildPhaseCargoCommand = "cargo test --workspace --all-targets --profile ${ciTestCargoProfile} --no-run --locked";
          doCheck = false;
        });

      individualCrateArgs =
        commonArgs
        // {
          inherit cargoArtifacts;
          buildType = "dev-optim";
          doCheck = false;
        };

      symbolicaCrateArgs = usesSymbolica:
        lib.optionalAttrs usesSymbolica {
          preBuild = licensePreCheck;
          SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
        };

      crateChecks = lib.listToAttrs (map (crate: {
          name = "crate-${crate.attr}";
          value = craneLib.buildPackage (individualCrateArgs
            // {
              pname = crate.package;
              inherit (apiMeta) version;
              src = fileSetForCrate crate.path;
              cargoExtraArgs = "--locked -p ${crate.package}";
            }
            // symbolicaCrateArgs crate.usesSymbolica);
        })
        workspaceCrates);

      gammaloop-cli = craneLib.buildPackage (individualCrateArgs
        // {
          pname = "gammaloop";
          inherit (apiMeta) version;
          src = fileSetForCrate ./crates/gammaloop-api;
          cargoExtraArgs = "--locked -p gammaloop-api --bin gammaloop";
        });

      gammaloopNextestArchive = craneLib.mkCargoDerivation (gammaloopNextestArgs
        // {
          inherit cargoArtifacts;
          pnameSuffix = "-nextest-archive";
          preCheck = licensePreCheck;
          SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
          doCheck = true;
          nativeBuildInputs = (gammaloopNextestArgs.nativeBuildInputs or []) ++ [pkgs.cargo-nextest];
          buildPhaseCargoCommand = ''
            mkdir -p "$out"
            cargo nextest --version
          '';
          checkPhaseCargoCommand = ''
            cargo nextest archive \
              --cargo-profile ${ciTestCargoProfile} \
              ${gammaloopNextestArgs.cargoExtraArgs} \
              --profile ci_gammaloop \
              --archive-format tar-zst \
              --archive-file "$out/archive.tar.zst"
          '';
          installPhaseCommand = ''
            test -f "$out/archive.tar.zst"
          '';
        });

      nextestPartitionArgs = partition: commonArgs
        // {
          pname = "gammaloop-nextest-partition-${toString partition}";
          inherit (apiMeta) version;
          cargoArtifacts = null;
          cargoVendorDir = null;
          doInstallCargoArtifacts = false;
          doCheck = true;
          preCheck = licensePreCheck;
          SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
          nativeBuildInputs = (commonArgs.nativeBuildInputs or []) ++ [pkgs.cargo-nextest];
          buildPhaseCargoCommand = ''
            mkdir -p "$out"
            cargo nextest --version
          '';
          checkPhaseCargoCommand = ''
            cargo nextest run \
              --archive-format tar-zst \
              --archive-file ${gammaloopNextestArchive}/archive.tar.zst \
              --workspace-remap . \
              --profile ci_gammaloop \
              --partition hash:${toString partition}/${toString ciPartitionCount} \
              --no-fail-fast \
              --final-status-level fail \
              --no-tests=pass \
              --run-ignored all
          '';
          installPhaseCommand = ''
            mkdir -p "$out"
          '';
        };

      partitionedNextestChecks = lib.listToAttrs (map (partition: {
          name = "gammaloop-nextest-partition-${toString partition}";
          value = craneLib.mkCargoDerivation (nextestPartitionArgs partition);
        })
        (lib.range 1 ciPartitionCount));

      impureCheckRunnerTargets =
        (map (crate: {
            runnerAttr = "nix-ci-check-${crate.attr}";
            checkAttr = "crate-${crate.attr}";
          })
          (builtins.filter (crate: crate.usesSymbolica) workspaceCrates))
        ++ [
          {
            runnerAttr = "nix-ci-check-gammaloop-nextest";
            checkAttr = "gammaloop-nextest";
          }
        ]
        ++ map (partition: {
          runnerAttr = "nix-ci-check-gammaloop-nextest-partition-${toString partition}";
          checkAttr = "gammaloop-nextest-partition-${toString partition}";
        }) (lib.range 1 ciPartitionCount);

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

          gammaloop-nextest = craneLib.cargoNextest (gammaloopNextestArgs
            // {
              inherit cargoArtifacts;
              preCheck = licensePreCheck;
              SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
              partitions = ciPartitionCount;
              partitionType = "hash";
              cargoNextestPartitionsExtraArgs = "--profile ci_gammaloop --no-fail-fast --final-status-level fail --no-tests=pass --run-ignored all";
            });
        }
        // crateChecks
        // partitionedNextestChecks;

      packages =
        {
          default = gammaloop-cli;
          gammaloop = gammaloop-cli;
          gammaloop-nextest-archive = gammaloopNextestArchive;
          inherit cargoArtifacts;
        }
        // impureCheckRunnerPackages
        // lib.optionalAttrs (!pkgs.stdenv.isDarwin) {
          gammaloop-llvm-coverage = craneLibLLvmTools.cargoLlvmCov (commonArgs
            // {
              inherit cargoArtifacts;
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

      ci = {
        partitionCount = toString ciPartitionCount;
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
