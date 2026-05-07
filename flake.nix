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

      workspaceRoot = ./.;

      cargoSources = craneLib.fileset.commonCargoSources workspaceRoot;

      cargoVendorDir = craneLib.vendorCargoDeps {
        cargoLock = ./Cargo.lock;
        overrideVendorGitCheckout = packages: drv:
          if lib.any (package: package.name == "symbolica") packages
          then
            drv.overrideAttrs (old: {
              postInstall =
                (old.postInstall or "")
                + ''
                  for crate in ${lib.concatMapStringsSep " " (package: lib.escapeShellArg "${package.name}-${package.version}") packages}; do
                    if [ -d "$out/$crate" ]; then
                      mkdir -p "$out/$crate/.git"
                      printf 'ref: refs/heads/nix-vendor\n' > "$out/$crate/.git/HEAD"
                    fi
                  done
                '';
            })
          else drv;
      };

      nonCargoBuildSources = lib.fileset.unions [
        ./.config
        ./assets
        ./crates/clinnet/templates
        ./crates/kurvst/typst/kurvst.wasm
        ./crates/kurvst/typst/src
        ./crates/kurvst/typst/typst.toml
        ./crates/linnest/typst/linnest.wasm
        ./crates/linnest/typst/src
        ./crates/linnest/typst/typst.toml
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

      workspaceFmtSrc = lib.fileset.toSource {
        root = workspaceRoot;
        fileset = cargoSources;
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

      linnestWasmSrc = lib.fileset.toSource {
        root = workspaceRoot;
        fileset = lib.fileset.unions [
          cargoSources
          ./crates/clinnet/templates/figure.typ
          ./crates/clinnet/templates/grid.typ
          ./crates/clinnet/templates/layout.typ
          ./crates/kurvst/typst/src
          ./crates/kurvst/typst/typst.toml
          ./crates/linnest/typst/src
          ./crates/linnest/typst/typst.toml
        ];
      };

      workspaceMemberDirs = let
        crateEntries = builtins.readDir ./crates;
        crateMemberDirs = map (name: "crates/${name}") (
          lib.filter (
            name:
              crateEntries.${name}
              == "directory"
              && builtins.pathExists (workspaceRoot + "/crates/${name}/Cargo.toml")
          ) (builtins.attrNames crateEntries)
        );
      in
        crateMemberDirs ++ ["tests"];

      autoCargoTargetDirs =
        lib.concatMap (
          member:
            lib.filter (
              dir: builtins.pathExists (workspaceRoot + "/${dir}")
            ) [
              "${member}/benches"
              "${member}/examples"
              "${member}/tests"
            ]
        )
        workspaceMemberDirs;

      autoCargoTargetPaths = lib.sort (left: right: left < right) (
        lib.concatMap (
          dir: let
            entries = builtins.readDir (workspaceRoot + "/${dir}");
          in
            map (name: "${dir}/${name}") (
              lib.filter (
                name:
                  entries.${name}
                  == "regular"
                  && lib.hasSuffix ".rs" name
                  && name != "mod.rs"
              ) (builtins.attrNames entries)
            )
        )
        autoCargoTargetDirs
      );

      dummyCargoTarget = pkgs.writeText "crane-dummy-cargo-target.rs" ''
        #![allow(clippy::all)]
        #![allow(dead_code)]

        pub fn main() {}
      '';

      src = workspaceBuildSrc;

      apiMeta = craneLib.crateNameFromCargoToml {
        cargoToml = ./crates/gammaloop-api/Cargo.toml;
      };

      clinnetMeta = craneLib.crateNameFromCargoToml {
        cargoToml = ./crates/clinnet/Cargo.toml;
      };

      linnestMeta = wasmCraneLib.crateNameFromCargoToml {
        cargoToml = ./crates/linnest/Cargo.toml;
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
        inherit cargoVendorDir;

        nativeBuildInputs =
          [
            pkgs.pkg-config
            pkgs.gcc
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

      ciCargoProfile = "dev-optim";

      ciArgs =
        commonArgs
        // {
          buildType = ciCargoProfile;
          CARGO_PROFILE = ciCargoProfile;
          # The workspace sets default-members to gammaloop-api, so CI checks must
          # opt into the full workspace explicitly.
          cargoExtraArgs = "--locked --workspace";

          PYO3_PYTHON = "${pkgs.python313}/bin/python3";
          PYTHONPATH = "${pkgs.python313}/lib/python3.13/site-packages";
        };

      licensePreCheck = ''
        if [ -z "''${SYMBOLICA_LICENSE:-}" ]; then
          echo "Missing SYMBOLICA_LICENSE environment variable" >&2
          exit 1
        fi
      '';

      # Crane's documented workspace pattern is to build one shared dependency cache and
      # reuse it across workspace lint/test/doc/package checks.
      cargoArtifacts = craneLib.buildDepsOnly (ciArgs
        // {
          pname = "gammaloop-workspace-deps";
          src = workspaceTestSrc;
          extraDummyScript =
            lib.concatMapStringsSep "\n" (path: ''
              install -D -m 0644 ${dummyCargoTarget} "$out/${path}"
            '')
            autoCargoTargetPaths;
        });

      symbolicaCrateArgs = usesSymbolica:
        lib.optionalAttrs usesSymbolica {
          preBuild = licensePreCheck;
          SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
        };

      linnestWasmArgs = {
        inherit cargoVendorDir;
        src = linnestWasmSrc;
        pname = "linnest-wasm";
        inherit (linnestMeta) version;
        strictDeps = true;
        doCheck = false;
        buildType = "release";
        CARGO_BUILD_TARGET = wasmTarget;
        cargoExtraArgs = "--locked -p linnest -p kurvst --features linnest/custom --target ${wasmTarget}";
      };

      linnestWasmCargoArtifacts = wasmCraneLib.buildDepsOnly (linnestWasmArgs
        // {
          pname = "linnest-wasm-deps";
        });

      linnest-wasm = wasmCraneLib.buildPackage (linnestWasmArgs
        // {
          cargoArtifacts = linnestWasmCargoArtifacts;
          cargoBuildCommand = "cargo build --release";
          installPhaseCommand = ''
            mkdir -p \
              "$out/templates" \
              "$out/templates/crates/linnest/typst" \
              "$out/templates/crates/kurvst/typst"
            cp "target/${wasmTarget}/release/linnest.wasm" "$out/linnest.wasm"
            cp "target/${wasmTarget}/release/kurvst.wasm" "$out/kurvst.wasm"
            cp crates/clinnet/templates/*.typ "$out/templates/"
            cp -R crates/linnest/typst/src "$out/templates/crates/linnest/typst/"
            cp crates/linnest/typst/typst.toml "$out/templates/crates/linnest/typst/typst.toml"
            cp -R crates/kurvst/typst/src "$out/templates/crates/kurvst/typst/"
            cp crates/kurvst/typst/typst.toml "$out/templates/crates/kurvst/typst/typst.toml"
            cp "$out/linnest.wasm" "$out/templates/crates/linnest/typst/linnest.wasm"
            cp "$out/kurvst.wasm" "$out/templates/crates/kurvst/typst/kurvst.wasm"
          '';
        });

      drawingTypstBundleAssets = ''
        mkdir -p crates/linnest/typst/src crates/kurvst/typst/src
        cp -R ${linnest-wasm}/templates/crates/linnest/typst/src/. crates/linnest/typst/src/
        cp ${linnest-wasm}/templates/crates/linnest/typst/typst.toml crates/linnest/typst/typst.toml
        cp ${linnest-wasm}/templates/crates/linnest/typst/linnest.wasm crates/linnest/typst/linnest.wasm
        cp -R ${linnest-wasm}/templates/crates/kurvst/typst/src/. crates/kurvst/typst/src/
        cp ${linnest-wasm}/templates/crates/kurvst/typst/typst.toml crates/kurvst/typst/typst.toml
        cp ${linnest-wasm}/templates/crates/kurvst/typst/kurvst.wasm crates/kurvst/typst/kurvst.wasm
      '';

      gammaloop-cli = craneLib.buildPackage (ciArgs
        // {
          inherit cargoArtifacts;
          buildType = ciCargoProfile;
          doCheck = false;
          pname = "gammaloop";
          inherit (apiMeta) version;
          cargoBuildCommand = "cargo build --profile ${ciCargoProfile}";
          cargoExtraArgs = "--locked -p gammaloop-api --bin gammaloop";
          preBuild = drawingTypstBundleAssets;
        });

      clinnetArgs = commonArgs
        // {
          buildType = ciCargoProfile;
          CARGO_PROFILE = ciCargoProfile;
          doCheck = false;
          pname = "clinnet";
          inherit (clinnetMeta) version;
          cargoBuildCommand = "cargo build --profile ${ciCargoProfile}";
          cargoExtraArgs = "--locked -p clinnet --bin linnet";
        };

      clinnetCargoArtifacts = craneLib.buildDepsOnly (clinnetArgs
        // {
          pname = "clinnet-deps";
          preBuild = drawingTypstBundleAssets;
        });

      clinnet-cli = craneLib.buildPackage (clinnetArgs
        // {
          cargoArtifacts = clinnetCargoArtifacts;
          preBuild = drawingTypstBundleAssets;
        });

      rscls = pkgs.rustPlatform.buildRustPackage rec {
        pname = "rscls";
        version = "0.2.3";
        src = pkgs.fetchCrate {
          inherit pname version;
          sha256 = "sha256-tahAhWCjhIVjbJ1NzrtiHBwGb/FBmUdK4XP9VlSPqh0=";
        };
        cargoHash = "sha256-JikjBTFeDh4XHBm57yiorsCwZhKikz0aiWNOTaMn0Vo=";
      };

      devShellPackages = with pkgs; [
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

          packages = devShellPackages ++ extraPackages;
        };

      impureCheckRunnerTargets = [
        {
          runnerAttr = "nix-ci-check-gammaloop-doctest";
          checkAttr = "gammaloop-doctest";
        }
        {
          runnerAttr = "nix-ci-check-gammaloop-nextest";
          checkAttr = "gammaloop-nextest";
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
              src = workspaceTestSrc;
              cargoClippyExtraArgs = "--all-targets -- --deny warnings";
            });

          gammaloop-doc = craneLib.cargoDoc (ciArgs
            // {
              inherit cargoArtifacts;
              src = workspaceTestSrc;
            });

          gammaloop-doctest = craneLib.cargoDocTest (ciArgs
            // {
              inherit cargoArtifacts;
              src = workspaceTestSrc;
            }
            // symbolicaCrateArgs true);

          gammaloop-fmt = craneLib.cargoFmt {
            src = workspaceFmtSrc;
            pname = "gammaloop-workspace";
            inherit (apiMeta) version;
          };

          linnest-wasm = pkgs.runCommand "linnest-wasm-check" {
            nativeBuildInputs = [pkgs.wasm-tools];
          } ''
            test -s ${linnest-wasm}/linnest.wasm
            test -s ${linnest-wasm}/kurvst.wasm
            test -s ${linnest-wasm}/templates/crates/linnest/typst/linnest.wasm
            test -s ${linnest-wasm}/templates/crates/kurvst/typst/kurvst.wasm
            cmp ${linnest-wasm}/linnest.wasm ${linnest-wasm}/templates/crates/linnest/typst/linnest.wasm
            cmp ${linnest-wasm}/kurvst.wasm ${linnest-wasm}/templates/crates/kurvst/typst/kurvst.wasm
            wasm-tools validate ${linnest-wasm}/linnest.wasm
            wasm-tools validate ${linnest-wasm}/kurvst.wasm
            test -s ${linnest-wasm}/templates/layout.typ
            test -s ${linnest-wasm}/templates/crates/linnest/typst/src/lib.typ
            test -s ${linnest-wasm}/templates/crates/linnest/typst/src/curve.typ
            test -s ${linnest-wasm}/templates/crates/kurvst/typst/src/lib.typ
            mkdir -p "$out"
          '';
        }
        // {
          gammaloop-nextest = craneLib.cargoNextest (ciArgs
            // {
              inherit cargoArtifacts;
              src = workspaceTestSrc;
              nativeBuildInputs = (ciArgs.nativeBuildInputs or []) ++ [pkgs.form];
              cargoNextestExtraArgs = "--profile ci_gammaloop --no-fail-fast --final-status-level fail --no-tests=pass --run-ignored all";
            }
            // {
              preCheck = ''
                ${licensePreCheck}
                export INSTA_WORKSPACE_ROOT="$PWD"
              '';
              SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
            });
        };

      packages =
        {
          default = gammaloop-cli;
          clinnet = clinnet-cli;
          gammaloop = gammaloop-cli;
          inherit clinnetCargoArtifacts linnest-wasm linnestWasmCargoArtifacts;
          inherit cargoArtifacts;
        }
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
        default = flake-utils.lib.mkApp {
          drv = gammaloop-cli;
        };
        gammaloop = flake-utils.lib.mkApp {
          drv = gammaloop-cli;
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
        default = mkDevShell [];
        full = mkDevShell [clinnet-cli rscls];
        clinnet = mkDevShell [clinnet-cli];
      };
    });
}
