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
        ./crates/idenso
        ./crates/linnet
        ./crates/spenso
      ];

      integrationTestTargetSources = lib.fileset.unions [
        ./tests/resources
        ./tests/tests
      ];

      nonIntegrationCargoSources = lib.fileset.difference cargoSources integrationTestTargetSources;

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

      workspaceNonIntegrationTestSrc = lib.fileset.toSource {
        root = workspaceRoot;
        fileset = lib.fileset.unions [
          nonIntegrationCargoSources
          nonCargoBuildSources
          snapshotSources
          ./tests/resources
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

      workspaceMemberPackageDirs =
        lib.listToAttrs (map (member: {
            name = (builtins.fromTOML (builtins.readFile (workspaceRoot + "/${member}/Cargo.toml"))).package.name;
            value = member;
          })
          workspaceMemberDirs);

      workspaceMemberPackages = builtins.attrNames workspaceMemberPackageDirs;

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

      workspaceDependencyManifestFiles =
        [
          ./Cargo.lock
          ./Cargo.toml
        ]
        ++ map (member: workspaceRoot + "/${member}/Cargo.toml") workspaceMemberDirs;

      workspaceDependencyBuildScripts =
        lib.filter builtins.pathExists (
          [./build.rs]
          ++ map (member: workspaceRoot + "/${member}/build.rs") workspaceMemberDirs
        );

      workspaceDependencySrc = lib.fileset.toSource {
        root = workspaceRoot;
        fileset = lib.fileset.unions (
          workspaceDependencyManifestFiles
          ++ workspaceDependencyBuildScripts
          ++ [
            ./tests/resources/fjcore
          ]
        );
      };

      crate2nixSourceRootFiles = [
        ./Cargo.lock
        ./Cargo.toml
      ];

      crate2nixSourceWith = fileset:
        lib.fileset.toSource {
          root = workspaceRoot;
          fileset = lib.fileset.unions (crate2nixSourceRootFiles ++ fileset);
        };

      crate2nixGammaloopApiSrc = crate2nixSourceWith [
        ./crates/gammaloop-api
        ./assets
      ];

      crate2nixGammalooprsSrc = crate2nixSourceWith [
        ./crates/gammalooprs
        ./assets
      ];

      crate2nixIntegrationTestsSrc = crate2nixSourceWith [
        ./tests
        ./assets
        ./examples/cli
      ];

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

      clinnetMeta = craneLib.crateNameFromCargoToml {
        cargoToml = ./crates/clinnet/Cargo.toml;
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
      nextestPython = pkgs.python313.withPackages (pythonPackages: [
        pythonPackages.numpy
      ]);

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

      ciCargoProfile = "ci-optim";

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

      crate2nixCommonOverride = old: {
        nativeBuildInputs = (old.nativeBuildInputs or []) ++ (commonArgs.nativeBuildInputs or []);
        buildInputs = (old.buildInputs or []) ++ (commonArgs.buildInputs or []);

        CC = nixCc;
        CXX = nixCxx;
        CARGO_CRATE_NAME = lib.replaceStrings [ "-" ] [ "_" ] old.crateName;
        "${cargoLinkerVar}" = nixCc;
        RUSTFLAGS = "-C linker=${nixCc}";

        LD_LIBRARY_PATH = runtimeLibPath;
        DYLD_LIBRARY_PATH = runtimeLibPath;
        PYO3_PYTHON = "${pkgs.python313}/bin/python3";
        PYTHONPATH = "${pkgs.python313}/lib/python3.13/site-packages";
        SYMBOLICA_OEM_LICENSE = "SYMBOLICA_OEM_GAMMALOOP";
      };

      crate2nixSourceOverride = old: src: workspace_member:
        (crate2nixCommonOverride old)
        // {
          inherit src workspace_member;
        };

      crate2nixSymbolicaGitHeadOverride = old:
        (crate2nixCommonOverride old)
        // {
          # GammaLoop only depends on Symbolica as a library. Declaring an empty
          # binary set keeps buildRustCrate from auto-detecting Symbolica's CLI
          # and putting it in the default out output for every feature set.
          crateBin = [];
          postPatch =
            (old.postPatch or "")
            + ''
              mkdir -p .git
              printf 'ref: refs/heads/nix-vendor\n' > .git/HEAD
            '';
        };

      crate2nixWorkspaceCrateOverrides = lib.genAttrs workspaceMemberPackages (_: crate2nixCommonOverride);

      crate2nixDefaultCrateOverrides =
        pkgs.defaultCrateOverrides
        // crate2nixWorkspaceCrateOverrides
        // {
          "gammaloop-api" = old: crate2nixSourceOverride old crate2nixGammaloopApiSrc "crates/gammaloop-api";
          "gammalooprs" = old: crate2nixSourceOverride old crate2nixGammalooprsSrc "crates/gammalooprs";
          "gammaloop-integration-tests" = old: crate2nixSourceOverride old crate2nixIntegrationTestsSrc "tests";
          "gmp-mpfr-sys" = crate2nixCommonOverride;
          "pyo3-build-config" = crate2nixCommonOverride;
          "rug" = old:
            (crate2nixCommonOverride old)
            // {
              DEP_GMP_LIMB_BITS = "64";
            };
          "symbolica" = crate2nixSymbolicaGitHeadOverride;
          "graphica" = crate2nixSymbolicaGitHeadOverride;
          "numerica" = crate2nixSymbolicaGitHeadOverride;
        };

      crate2nixBuildRustCrateForPkgs = pkgs':
        pkgs'.buildRustCrate.override {
          rustc = ciToolchain;
          cargo = ciToolchain;
        };
      crate2nixBuildRustCrateWithDefaultOverridesForPkgs = pkgs':
        (crate2nixBuildRustCrateForPkgs pkgs').override {
          defaultCrateOverrides = crate2nixDefaultCrateOverrides;
        };

      crate2nixPackageSet = import ./Cargo.nix {
        inherit pkgs;
        release = true;
        buildRustCrateForPkgs = crate2nixBuildRustCrateForPkgs;
        defaultCrateOverrides = crate2nixDefaultCrateOverrides;
      };

      crate2nixWorkspaceMembers = crate2nixPackageSet.workspaceMembers;
      crate2nixBuildDefault = package: crate2nixWorkspaceMembers.${package}.build;
      crate2nixBuildWithFeatures = package: features:
        (crate2nixBuildDefault package).override {
          inherit features;
        };
      crate2nixRootFeatures = package: crate2nixPackageSet.internal.crates.${package}.features or {};
      crate2nixRootFeatureIf = package: feature:
        lib.optional (builtins.hasAttr feature (crate2nixRootFeatures package)) feature;
      crate2nixCiCommonFeaturesFor = package:
        ["default"] ++ lib.concatMap (crate2nixRootFeatureIf package) ["bincode" "serde"];
      crate2nixCiExtraFeatureSets = {
        "gammaloop-tracing-filter" = ["clap" "symbolica"];
        idenso = ["reference-cases"];
        linnet = ["symbolica"];
        spenso = ["shadowing"];
        "spenso-macros" = ["shadowing"];
      };
      crate2nixCiFeaturesFor = package:
        sortedUnique (crate2nixCiCommonFeaturesFor package ++ (crate2nixCiExtraFeatureSets.${package} or []));
      crate2nixBuild = package: crate2nixBuildWithFeatures package (crate2nixCiFeaturesFor package);
      crate2nixCiSanitizePackageId = packageId: lib.replaceStrings [" "] ["-"] packageId;
      crate2nixCiFeatureSetName = features:
        if features == []
        then "plain"
        else lib.concatStringsSep "-" features;
      crate2nixCiPrebuildCrateFor = packageId: features:
        (crate2nixPackageSet.internal.builtRustCratesWithFeatures {
          inherit packageId features;
          buildRustCrateForPkgsFunc = crate2nixBuildRustCrateWithDefaultOverridesForPkgs;
          runTests = false;
        })
        .crates.${packageId};
      crate2nixCiDependencyEntriesFor = packageId: featureSets:
        map (features: {
          name = "${crate2nixCiSanitizePackageId packageId}-${crate2nixCiFeatureSetName features}";
          path = crate2nixCiPrebuildCrateFor packageId features;
        })
        featureSets;

      gammaloop-cli = crate2nixBuildWithFeatures "gammaloop-api" ["default"];
      gammaloop-python-lib = crate2nixBuildWithFeatures "gammaloop-api" [
        "python_abi"
        "pyo3-extension-module"
      ];
      gammaloop-python-lib-output = lib.getLib gammaloop-python-lib;
      pythonSitePackages = "${pkgs.python313.sitePackages}";
      gammaloop-python-module = pkgs.runCommand "gammaloop-python-module" {} ''
        mkdir -p "$out/${pythonSitePackages}/gammaloop"
        cp ${crate2nixGammaloopApiSrc}/crates/gammaloop-api/python/gammaloop/__init__.py \
          "$out/${pythonSitePackages}/gammaloop/__init__.py"

        extension="$(
          find ${gammaloop-python-lib-output} -type f \
            \( -name 'libgammaloop_api*.so' -o -name 'gammaloop_api*.so' -o -name 'libgammaloop_api*.dylib' -o -name 'gammaloop_api*.dylib' \) \
            | sort \
            | head -n 1
        )"
        if [ -z "$extension" ]; then
          echo "Could not find crate2nix-built gammaloop-api Python extension in ${gammaloop-python-lib-output}" >&2
          exit 1
        fi
        cp "$extension" "$out/${pythonSitePackages}/gammaloop/_gammaloop.so"
      '';
      clinnet-cli = crate2nixBuildWithFeatures "clinnet" ["default"];

      crate2nixPackageOutputs = lib.listToAttrs (map (package: {
          name = "crate-${package}";
          value = crate2nixBuild package;
        })
        workspaceMemberPackages);

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
      # Keep this input to manifests and build-script inputs so source-only commits
      # can reuse the dependency artifact from the NixCI cache.
      cargoArtifacts = craneLib.buildDepsOnly (ciArgs
        // {
          pname = "gammaloop-workspace-deps";
          src = workspaceDependencySrc;
          extraDummyScript =
            lib.concatMapStringsSep "\n" (path: ''
              install -D -m 0644 ${dummyCargoTarget} "$out/${path}"
            '')
            autoCargoTargetPaths;
        });

      workspaceBuildArtifacts = craneLib.cargoBuild (ciArgs
        // {
          inherit cargoArtifacts;
          pname = "gammaloop-workspace-build-artifacts";
          src = workspaceNonIntegrationTestSrc;
          cargoExtraArgs = "--locked --workspace --exclude gammaloop-integration-tests --tests";
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

      nextestPackageGroups = [
        {
          name = "core";
          packages = [
            "gammaloop-api"
            "gammaloop-tracing-filter"
            "gammaloop-tracing-filter-macros"
            "gammalooprs"
          ];
        }
        {
          name = "integration";
          packages = ["gammaloop-integration-tests"];
        }
        {
          name = "python-api";
          packages = ["gammaloop-integration-tests"];
          filter = "package(gammaloop-integration-tests) & binary(test_python_api)";
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

      checkedNextestPackageGroups = assert lib.asserts.assertMsg (
        missingNextestPackages == [] && extraNextestPackages == []
      ) "nextest split package coverage mismatch: missing [${lib.concatStringsSep ", " missingNextestPackages}], extra [${lib.concatStringsSep ", " extraNextestPackages}]"; nextestPackageGroups;

      nextestTargetTriple = pkgs.stdenv.hostPlatform.rust.rustcTargetSpec or pkgs.stdenv.hostPlatform.config;
      nextestRustLibDir = "${ciToolchain}/lib/rustlib/${nextestTargetTriple}/lib";

      nextestBaseExtraArgs = "--profile ${nextestProfile} --no-fail-fast --final-status-level fail --no-tests=pass";

      nextestPackageFilter = packages: "-E ${lib.escapeShellArg (lib.concatMapStringsSep " | " (package: "package(${package})") packages)}";
      nextestFilterFor = target:
        if target ? filter
        then "-E ${lib.escapeShellArg target.filter}"
        else nextestPackageFilter target.packages;
      nextestUsesIntegrationTests = target: builtins.elem "gammaloop-integration-tests" target.packages;
      nextestSrcFor = target:
        if nextestUsesIntegrationTests target
        then workspaceTestSrc
        else workspaceNonIntegrationTestSrc;

      crate2nixReplacePackageDependency = dependencyName: newPackageId: dependencies:
        map (
          dependency:
            if dependency.name == dependencyName
            then dependency // {packageId = newPackageId;}
            else dependency
        )
        dependencies;

      crate2nixTestSharedSources = [
        ./.config
        ./assets
        ./crates/clinnet/templates
        ./crates/vakint/form_src
        ./crates/vakint/templates
        ./examples/cli
        ./tests/resources
      ];

      crate2nixTestSrcForPackage = package:
        crate2nixSourceWith (
          [
            (workspaceRoot + "/${workspaceMemberPackageDirs.${package}}")
          ]
          ++ crate2nixTestSharedSources
        );

      crate2nixTestCrateConfigs = let
        crates = crate2nixPackageSet.internal.crates;
        workspaceTestCrates = lib.genAttrs workspacePackages (package:
          crates.${package}
          // {
            src = crate2nixTestSrcForPackage package;
            workspace_member = workspaceMemberPackageDirs.${package};
          });
      in
        crates
        // workspaceTestCrates
        // {
          "linnet-nontest" = workspaceTestCrates."linnet";
          "linnest-for-linnet-test" =
            workspaceTestCrates."linnest"
            // {
              dependencies = crate2nixReplacePackageDependency "linnet" "linnet-nontest" workspaceTestCrates."linnest".dependencies;
            };
          "linnet" =
            workspaceTestCrates."linnet"
            // {
              devDependencies = crate2nixReplacePackageDependency "linnest" "linnest-for-linnet-test" workspaceTestCrates."linnet".devDependencies;
            };
          "spenso-macros-nontest" = workspaceTestCrates."spenso-macros";
          "spenso-for-spenso-macros-test" =
            workspaceTestCrates."spenso"
            // {
              dependencies = crate2nixReplacePackageDependency "spenso-macros" "spenso-macros-nontest" workspaceTestCrates."spenso".dependencies;
            };
          "spenso-macros" =
            workspaceTestCrates."spenso-macros"
            // {
              devDependencies = crate2nixReplacePackageDependency "spenso" "spenso-for-spenso-macros-test" workspaceTestCrates."spenso-macros".devDependencies;
            };
        };

      crate2nixTestBuildRustCrateForPkgs = pkgs':
        crate2nixBuildRustCrateWithDefaultOverridesForPkgs pkgs';

      crate2nixBuiltTestCratesFor = package:
        crate2nixPackageSet.internal.builtRustCratesWithFeatures {
          packageId = package;
          features = crate2nixCiFeaturesFor package;
          crateConfigs = crate2nixTestCrateConfigs;
          buildRustCrateForPkgsFunc = crate2nixTestBuildRustCrateForPkgs;
          runTests = true;
        };

      crate2nixTestBinaryCrateFor = package:
        (crate2nixBuiltTestCratesFor package).crates.${package}.override (old: {
          buildTests = true;
          preBuild =
            (old.preBuild or "")
            + lib.optionalString (crate2nixTestCrateConfigs.${package}.procMacro or false) ''
              # buildRustCrate wires the current crate as an rlib for integration
              # tests. Proc-macro crates emit a build-platform shared library
              # instead, so correct the self --extern path for those tests.
              build_lib() {
                lib_src=$1
                echo_build_heading $lib_src "$LIB_NAME"

                noisily rustc \
                  --crate-name $CRATE_NAME \
                  $lib_src \
                  --out-dir target/lib \
                  -L dependency=target/deps \
                  --cap-lints allow \
                  $LINK \
                  $EXTRA_LINK_ARGS \
                  $EXTRA_LINK_ARGS_LIB \
                  $LIB_RUSTC_OPTS \
                  $BUILD_OUT_DIR \
                  $EXTRA_BUILD \
                  $EXTRA_FEATURES \
                  $EXTRA_RUSTC_FLAGS \
                  --color $colors

                if [ -e target/lib/lib$CRATE_NAME-$metadata$LIB_EXT ]; then
                  EXTRA_LIB=" --extern $CRATE_NAME=target/lib/lib$CRATE_NAME-$metadata$LIB_EXT"
                else
                  EXTRA_LIB=" --extern $CRATE_NAME=target/lib/lib$CRATE_NAME-$metadata.rlib"
                fi
              }
            '';
        });

      crate2nixTestBinaryCrates = lib.genAttrs workspacePackages crate2nixTestBinaryCrateFor;

      crate2nixTestBinaryPackageOutputs = lib.listToAttrs (map (package: {
          name = "crate-test-binaries-${package}";
          value = crate2nixTestBinaryCrates.${package};
        })
        workspacePackages);

      crate2nixCiSymbolicaFeatureSets = [
        ["bincode" "gmp" "serde" "tracing_max_level_info"]
        ["bincode" "gmp" "python_export" "serde"]
      ];
      crate2nixCiPrebuild = pkgs.linkFarm "gammaloop-crate2nix-ci-prebuild" (
        crate2nixCiDependencyEntriesFor "symbolica" crate2nixCiSymbolicaFeatureSets
      );

      nextestCargoMetadata = pkgs.runCommand "gammaloop-nextest-cargo-metadata.json" {
        nativeBuildInputs = [ciToolchain];
      } ''
        cp -R ${workspaceTestSrc} source
        chmod -R u+w source
        cd source
        cargo metadata --format-version 1 --no-deps > "$out"
      '';

      nextestBinarySetFor = target: let
        packageNamesFile = pkgs.writeText "gammaloop-nextest-packages-${target.name}.json" (builtins.toJSON target.packages);
        testCratePathsFile = pkgs.writeText "gammaloop-nextest-test-crates-${target.name}.json" (builtins.toJSON (lib.genAttrs target.packages (package: "${crate2nixTestBinaryCrates.${package}}")));
      in
        pkgs.runCommand "gammaloop-nextest-binaries-${target.name}" {
          nativeBuildInputs = [pkgs.python3];
          packageNames = packageNamesFile;
          testCratePaths = testCratePathsFile;
          cargoMetadata = nextestCargoMetadata;
          rustLibDir = nextestRustLibDir;
          targetTriple = nextestTargetTriple;
        } ''
          mkdir -p "$out/target/debug/deps"
          export OUT_DIR="$out"
          python3 <<'PY'
          import json
          import os
          import shutil
          import stat
          import sys
          from pathlib import Path

          out = Path(os.environ["OUT_DIR"])
          target_dir = out / "target"
          deps_dir = target_dir / "debug" / "deps"

          with open(os.environ["cargoMetadata"], "r", encoding="utf-8") as handle:
              metadata = json.load(handle)
          with open(os.environ["packageNames"], "r", encoding="utf-8") as handle:
              package_names = json.load(handle)
          with open(os.environ["testCratePaths"], "r", encoding="utf-8") as handle:
              test_crate_paths = json.load(handle)

          metadata["target_directory"] = str(target_dir)
          metadata["build_directory"] = str(target_dir)
          with open(out / "cargo-metadata.json", "w", encoding="utf-8") as handle:
              json.dump(metadata, handle, separators=(",", ":"))

          packages_by_name = {package["name"]: package for package in metadata["packages"]}

          def nextest_kind(kinds):
              if any(kind in kinds for kind in ("lib", "rlib", "cdylib", "proc-macro")):
                  return "lib"
              for kind in ("bin", "test", "bench", "example"):
                  if kind in kinds:
                      return kind
              return None

          def executable_files(test_dir):
              if not test_dir.exists():
                  return []
              return [path for path in sorted(test_dir.iterdir()) if path.is_file() and os.access(path, os.X_OK)]

          def find_binary(files, package_name, target_name, kind):
              names = [target_name, target_name.replace("-", "_")]
              if kind == "lib":
                  names.extend([package_name, package_name.replace("-", "_")])
              candidates = []
              for name in dict.fromkeys(names):
                  candidates.extend(path for path in files if path.name == name or path.name.startswith(f"{name}-"))
              unique = []
              seen = set()
              for candidate in candidates:
                  if candidate not in seen:
                      unique.append(candidate)
                      seen.add(candidate)
              exact = [path for path in unique if path.name in names]
              if len(exact) == 1:
                  return exact[0]
              if len(unique) == 1:
                  return unique[0]
              if unique:
                  return sorted(unique, key=lambda path: (len(path.name), path.name))[0]
              return None

          rust_binaries = {}
          missing = []
          for package_name in package_names:
              package = packages_by_name[package_name]
              test_dir = Path(test_crate_paths[package_name]) / "tests"
              files = executable_files(test_dir)
              for target in package["targets"]:
                  if not target.get("test", False):
                      continue
                  kind = nextest_kind(target["kind"])
                  if kind is None:
                      continue
                  source = find_binary(files, package_name, target["name"], kind)
                  if source is None:
                      missing.append(f"{package_name}:{target['name']} ({'/'.join(target['kind'])}) in {test_dir}")
                      continue
                  binary_id = f"{package_name}::{target['name']}"
                  safe_binary_id = binary_id.replace("/", "_").replace(":", "_")
                  destination = deps_dir / f"{safe_binary_id}-{source.name}"
                  shutil.copy2(source, destination)
                  destination.chmod(destination.stat().st_mode | stat.S_IXUSR)
                  rust_binaries[binary_id] = {
                      "binary-id": binary_id,
                      "binary-name": target["name"],
                      "package-id": package["id"],
                      "kind": kind,
                      "binary-path": str(destination),
                      "build-platform": "target",
                  }

          if missing:
              print("missing crate2nix test binaries for nextest metadata:", file=sys.stderr)
              for item in missing:
                  print(f"  - {item}", file=sys.stderr)
              sys.exit(1)

          build_meta = {
              "target-directory": str(target_dir),
              "base-output-directories": ["debug"],
              "non-test-binaries": {},
              "build-script-out-dirs": {},
              "linked-paths": [],
              "platforms": {
                  "host": {
                      "platform": {
                          "triple": os.environ["targetTriple"],
                          "target-features": "unknown",
                      },
                      "libdir": {
                          "status": "available",
                          "path": os.environ["rustLibDir"],
                      },
                  },
                  "targets": [],
              },
              "target-platforms": [
                  {
                      "triple": os.environ["targetTriple"],
                      "target-features": "unknown",
                  }
              ],
              "target-platform": None,
          }
          binaries_metadata = {
              "rust-build-meta": build_meta,
              "rust-binaries": rust_binaries,
          }
          with open(out / "binaries-metadata.json", "w", encoding="utf-8") as handle:
              json.dump(binaries_metadata, handle, separators=(",", ":"))
          PY
        '';

      nextestBinarySets = lib.listToAttrs (map (target: {
          name = "gammaloop-nextest-binaries-${target.name}";
          value = nextestBinarySetFor target;
        })
        checkedNextestPackageGroups);

      nextestBinarySetForTarget = target: nextestBinarySets."gammaloop-nextest-binaries-${target.name}";

      nextestBinarySetAggregate = pkgs.linkFarm "gammaloop-nextest-binaries" (map (target: {
          name = target.name;
          path = nextestBinarySetForTarget target;
        })
        checkedNextestPackageGroups);

      nextestCheckFor = target:
        pkgs.stdenv.mkDerivation ({
          pname = "gammaloop-nextest-${target.name}";
          version = "0.1.0";
          src = nextestSrcFor target;
          nativeBuildInputs =
            [ciToolchain pkgs.cargo-nextest pkgs.form]
            ++ lib.optionals (nextestUsesIntegrationTests target) [nextestPython];
          dontConfigure = true;
          dontBuild = true;
          doCheck = true;
          dontFixup = true;
          LD_LIBRARY_PATH = runtimeLibPath;
          DYLD_LIBRARY_PATH = runtimeLibPath;
          RUST_BACKTRACE = "1";
          RUST_LIB_BACKTRACE = "1";
          SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
          preCheck = ''
            ${licensePreCheck}
            # crate2nix-built test binaries run under nextest with the workspace
            # root as cwd, while insta snapshots are stored under each crate src.
            # Mirror those snapshot directories into the workspace-level src tree
            # in the disposable Nix build directory so insta can find them.
            if [ -d crates ]; then
              while IFS= read -r snapshots; do
                rel="''${snapshots#crates/*/src/}"
                mkdir -p "src/$rel"
                cp -R "$snapshots/." "src/$rel/"
              done < <(find crates -type d -name snapshots | sort)
            fi
          '';
          checkPhase = ''
            runHook preCheck

            cp -R ${nextestBinarySetForTarget target}/target target
            chmod -R u+w target

            cat > nextest-nix.toml <<EOF
            [store]
            dir = "$PWD/target/nextest"

            EOF
            cat "$PWD/.config/nextest.toml" >> nextest-nix.toml

            set +e
            cargo nextest run \
              --cargo-metadata ${nextestBinarySetForTarget target}/cargo-metadata.json \
              --binaries-metadata ${nextestBinarySetForTarget target}/binaries-metadata.json \
              --target-dir-remap "$PWD/target" \
              --workspace-remap "$PWD" \
              --config-file "$PWD/nextest-nix.toml" \
              ${nextestFilterFor target} ${nextestBaseExtraArgs}
            nextest_status=$?
            set -e

            ${nextestFailureSummary}/bin/nextest-failure-summary ${lib.escapeShellArg nextestJunitPath} || true

            runHook postCheck
            if [ "$nextest_status" -ne 0 ]; then
              exit "$nextest_status"
            fi
          '';
          installPhase = ''
            mkdir -p "$out"
          '';
        } // lib.optionalAttrs (nextestUsesIntegrationTests target) {
          PYO3_PYTHON = "${nextestPython}/bin/python3";
          PYTHON = "${nextestPython}/bin/python3";
          PYTHONPATH = "${gammaloop-python-module}/${pythonSitePackages}:${nextestPython}/${pythonSitePackages}";
        });

      nextestRunChecks = lib.listToAttrs (map (target: {
          name = "gammaloop-nextest-${target.name}";
          value = nextestCheckFor target;
        })
        checkedNextestPackageGroups);

      nextestChecks =
        {
          gammaloop-nextest-binaries = nextestBinarySetAggregate;
        }
        // nextestBinarySets
        // nextestRunChecks;

      nextestAggregate = pkgs.runCommand "gammaloop-nextest" {} ''
        ${lib.concatMapStringsSep "\n" (check: "test -d ${check}") (builtins.attrValues nextestRunChecks)}
        mkdir -p "$out"
      '';

      impureCheckRunnerTargets =
        [
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

          linnest-wasm =
            pkgs.runCommand "linnest-wasm-check" {
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
          inherit clinnet-cli;
          "gammaloop-python-module" = gammaloop-python-module;
          inherit linnest-wasm linnestWasmCargoArtifacts;
          "crate2nix-ci-prebuild" = crate2nixCiPrebuild;
          inherit cargoArtifacts workspaceBuildArtifacts;
          "nix-ci-passed" = nixCiPassed;
        }
        // crate2nixPackageOutputs
        // crate2nixTestBinaryPackageOutputs
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
          exePath = "/bin/gammaloop";
        };
        gammaloop = flake-utils.lib.mkApp {
          drv = gammaloop-cli;
          exePath = "/bin/gammaloop";
        };
      };

      devShells.default = craneLib.devShell ({
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

          packages = with pkgs;
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
              clinnet-cli
            ]
            ++ lib.optionals (!pkgs.stdenv.isDarwin) [
              valgrind
            ];
        });
    });
}
