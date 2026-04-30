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

      workspaceMemberPackages = map (
        member:
          (builtins.fromTOML (builtins.readFile (workspaceRoot + "/${member}/Cargo.toml"))).package.name
      )
      workspaceMemberDirs;

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

      nextestProfile = "ci_gammaloop";
      nextestJunitPath = "target/nextest/${nextestProfile}/junit.xml";

      nextestFailureSummary = pkgs.writeShellApplication {
        name = "nextest-failure-summary";
        runtimeInputs = [pkgs.python313];
        text = ''
          exec python3 - "$@" <<'PY'
          import pathlib
          import sys
          import xml.etree.ElementTree as ET

          MAX_FAILURES = 20
          MAX_LINES_PER_FAILURE = 120


          def local_name(tag):
              return tag.rsplit("}", 1)[-1]


          def int_attr(node, name):
              try:
                  return int(node.attrib.get(name, "0"))
              except ValueError:
                  return 0


          def clean_text(text):
              if not text:
                  return ""
              lines = text.replace("\r\n", "\n").replace("\r", "\n").splitlines()
              return "\n".join(line.rstrip() for line in lines).strip()


          def unique_chunks(chunks):
              seen = set()
              result = []
              for chunk in chunks:
                  chunk = clean_text(chunk)
                  if chunk and chunk not in seen:
                      seen.add(chunk)
                      result.append(chunk)
              return result


          def failing_cases(root):
              failure_tags = {
                  "failure",
                  "error",
                  "rerunFailure",
                  "rerunError",
                  "flakyFailure",
                  "flakyError",
              }
              for case in root.iter():
                  if local_name(case.tag) != "testcase":
                      continue
                  failures = [
                      child
                      for child in list(case)
                      if local_name(child.tag) in failure_tags
                  ]
                  if failures:
                      yield case, failures


          def test_name(case):
              classname = case.attrib.get("classname", "").strip()
              name = case.attrib.get("name", "").strip()
              if classname and name:
                  return f"{classname}::{name}"
              return name or classname or "<unknown test>"


          def failure_chunks(case, failures):
              chunks = []
              for failure in failures:
                  label = local_name(failure.tag)
                  kind = failure.attrib.get("type", "").strip()
                  header = f"{label}: {kind}" if kind else label
                  text = clean_text(failure.text)
                  chunks.append(f"{header}\n{text}" if text else header)

              for child in list(case):
                  if local_name(child.tag) in {"system-err", "system-out"}:
                      text = clean_text(child.text)
                      if text:
                          chunks.append(f"{local_name(child.tag)}\n{text}")

              return unique_chunks(chunks)


          path = pathlib.Path(sys.argv[1])
          if not path.exists():
              print(f"nextest JUnit report not found: {path}", file=sys.stderr)
              sys.exit(0)

          try:
              root = ET.parse(path).getroot()
          except ET.ParseError as error:
              print(f"could not parse nextest JUnit report {path}: {error}", file=sys.stderr)
              sys.exit(0)

          failures = list(failing_cases(root))
          total = int_attr(root, "tests")
          failure_count = int_attr(root, "failures")
          error_count = int_attr(root, "errors")

          print()
          print("========== nextest test summary ==========")
          print(f"report: {path}")
          print(f"tests: {total}, failures: {failure_count}, errors: {error_count}")

          if not failures:
              print("No failing testcases found.")
              print("=========================================")
              sys.exit(0)

          for index, (case, case_failures) in enumerate(failures[:MAX_FAILURES], start=1):
              print()
              print(f"{index}. {test_name(case)}")
              for chunk in failure_chunks(case, case_failures):
                  lines = chunk.splitlines()
                  truncated = len(lines) > MAX_LINES_PER_FAILURE
                  if truncated:
                      omitted = len(lines) - MAX_LINES_PER_FAILURE
                      lines = lines[:MAX_LINES_PER_FAILURE]
                  for line in lines:
                      print(f"   {line}")
                  if truncated:
                      print(f"   ... omitted {omitted} more lines; see full nextest output above")

          if len(failures) > MAX_FAILURES:
              print()
              print(f"... omitted {len(failures) - MAX_FAILURES} more failing tests")

          print("=========================================")
          PY
        '';
      };

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
        cargoExtraArgs = "--locked -p linnest --features custom --target ${wasmTarget}";
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
            mkdir -p "$out/templates"
            cp "target/${wasmTarget}/release/linnest.wasm" "$out/linnest.wasm"
            cp crates/clinnet/templates/*.typ "$out/templates/"
            cp "$out/linnest.wasm" "$out/templates/linnest.wasm"
          '';
        });

      gammaloop-cli = craneLib.buildPackage (ciArgs
        // {
          inherit cargoArtifacts;
          buildType = ciCargoProfile;
          doCheck = false;
          pname = "gammaloop";
          inherit (apiMeta) version;
          cargoBuildCommand = "cargo build --profile ${ciCargoProfile}";
          cargoExtraArgs = "--locked -p gammaloop-api --bin gammaloop";
        });

      nextestPackageGroups = [
        {
          name = "core";
          packages = [
            "gammaloop-api"
            "gammaloop-tracing-filter"
            "gammalooprs"
          ];
        }
        {
          name = "integration";
          packages = ["gammaloop-integration-tests"];
        }
        {
          name = "linnet";
          packages = [
            "clinnet"
            "linnet"
            "linnet-py"
            "linnest"
          ];
        }
        {
          name = "spenso";
          packages = [
            "idenso"
            "spenso"
            "spenso-hep-lib"
            "spenso-macros"
            "spynso3"
          ];
        }
        {
          name = "vakint";
          packages = ["vakint"];
        }
      ];

      sortedUnique = list: lib.sort (left: right: left < right) (lib.unique list);

      nextestSplitPackages = sortedUnique (lib.concatMap (target: target.packages) nextestPackageGroups);
      workspacePackages = sortedUnique workspaceMemberPackages;
      missingNextestPackages = lib.subtractLists nextestSplitPackages workspacePackages;
      extraNextestPackages = lib.subtractLists workspacePackages nextestSplitPackages;

      checkedNextestPackageGroups =
        assert lib.asserts.assertMsg (
          missingNextestPackages == [] && extraNextestPackages == []
        ) "nextest split package coverage mismatch: missing [${lib.concatStringsSep ", " missingNextestPackages}], extra [${lib.concatStringsSep ", " extraNextestPackages}]";
          nextestPackageGroups;

      nextestBaseExtraArgs = "--profile ${nextestProfile} --no-fail-fast --final-status-level fail --no-tests=pass --run-ignored all";

      nextestPackageArgs = packages:
        lib.concatMapStringsSep " " (package: "-p ${package}") packages;

      nextestCheckFor = target:
        craneLib.cargoNextest (ciArgs
          // {
            inherit cargoArtifacts;
            src = workspaceTestSrc;
            pname = "gammaloop-nextest-${target.name}";
            nativeBuildInputs = (ciArgs.nativeBuildInputs or []) ++ [pkgs.form];
            cargoExtraArgs = "--locked ${nextestPackageArgs target.packages}";
            cargoNextestExtraArgs = nextestBaseExtraArgs;
            checkPhase = ''
              runHook preCheck

              set +e
              cargo nextest run \
                ''${CARGO_PROFILE:+--cargo-profile $CARGO_PROFILE} \
                --locked ${nextestPackageArgs target.packages} ${nextestBaseExtraArgs}
              nextest_status=$?
              set -e

              ${nextestFailureSummary}/bin/nextest-failure-summary ${lib.escapeShellArg nextestJunitPath} || true

              runHook postCheck
              exit "$nextest_status"
            '';
            doInstallCargoArtifacts = false;
            RUST_BACKTRACE = "1";
            RUST_LIB_BACKTRACE = "1";
          }
          // {
            preCheck = ''
              ${licensePreCheck}
              export INSTA_WORKSPACE_ROOT="$PWD"
            '';
            SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
          });

      nextestChecks = lib.listToAttrs (map (target: {
          name = "gammaloop-nextest-${target.name}";
          value = nextestCheckFor target;
        })
        checkedNextestPackageGroups);

      nextestAggregate = pkgs.runCommand "gammaloop-nextest" {} ''
        ${lib.concatMapStringsSep "\n" (check: "test -d ${check}") (builtins.attrValues nextestChecks)}
        mkdir -p "$out"
      '';

      impureCheckRunnerTargets = [
        {
          runnerAttr = "nix-ci-check-gammaloop-doctest";
          checkAttr = "gammaloop-doctest";
        }
        {
          runnerAttr = "nix-ci-check-gammaloop-nextest";
          checkAttr = "gammaloop-nextest";
        }
      ]
      ++ map (target: {
        runnerAttr = "nix-ci-check-gammaloop-nextest-${target.name}";
        checkAttr = "gammaloop-nextest-${target.name}";
      })
      checkedNextestPackageGroups;

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

      nixCiPassed = pkgs.writeShellApplication {
        name = "nix-ci-passed";
        text = ''
          echo "All NixCI build and test jobs passed."
        '';
      };

      gungraunRunner = pkgs.rustPlatform.buildRustPackage rec {
        pname = "gungraun-runner";
        version = "0.18.2";
        src = pkgs.fetchCrate {
          inherit pname version;
          hash = "sha256-DiJq9TZCZdWKSstIyMjkLuxaYXua0WKD2AVbEIxM590=";
        };
        cargoHash = "sha256-eb9U1MgCg7MpwzS2RnFXMWdPitweKMMty0n3SC0F6+I=";
        doCheck = false;
      };
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
            test -s ${linnest-wasm}/templates/linnest.wasm
            cmp ${linnest-wasm}/linnest.wasm ${linnest-wasm}/templates/linnest.wasm
            wasm-tools validate ${linnest-wasm}/linnest.wasm
            test -s ${linnest-wasm}/templates/layout.typ
            mkdir -p "$out"
          '';
        }
        // nextestChecks
        // {
          gammaloop-nextest = nextestAggregate;
        };

      packages =
        {
          default = gammaloop-cli;
          gammaloop = gammaloop-cli;
          inherit linnest-wasm linnestWasmCargoArtifacts;
          inherit cargoArtifacts;
          "nix-ci-passed" = nixCiPassed;
        }
        // impureCheckRunnerPackages
        // lib.optionalAttrs (!pkgs.stdenv.isDarwin) {
          inherit gungraunRunner;
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
        ]
        ++ lib.optionals (!pkgs.stdenv.isDarwin) [
          gungraunRunner
          valgrind
        ]
        ++ [
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
