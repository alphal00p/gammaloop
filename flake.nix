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
      workspaceGraph = builtins.fromJSON (builtins.readFile ./nix/ci-workspace-graph.json);
      workspaceHackPackage = "gammaloop-workspace-hack";

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

      workspaceManifestFor = member:
        builtins.fromTOML (builtins.readFile (workspaceRoot + "/${member}/Cargo.toml"));

      workspaceMemberPackageDirs =
        lib.listToAttrs (map (member: {
            name = (workspaceManifestFor member).package.name;
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

      workspaceExplicitCargoTargetPaths =
        lib.concatMap (
          member: let
            manifest = workspaceManifestFor member;
            targetPaths = targets:
              lib.concatMap (
                target:
                  lib.optional (target ? path) "${member}/${target.path}"
              )
              targets;
          in
            lib.optional ((manifest ? lib) && (manifest.lib ? path)) "${member}/${manifest.lib.path}"
            ++ targetPaths (manifest.bin or [])
            ++ targetPaths (manifest.example or [])
            ++ targetPaths (manifest.test or [])
            ++ targetPaths (manifest.bench or [])
        )
        workspaceMemberDirs;

      workspaceDefaultCargoTargetPaths =
        lib.concatMap (
          member: [
            "${member}/src/lib.rs"
            "${member}/src/main.rs"
          ]
        )
        workspaceMemberDirs;

      workspaceCargoTargetRelPaths =
        lib.filter (path: builtins.pathExists (workspaceRoot + "/${path}")) (
          lib.sort (left: right: left < right) (lib.unique (
            workspaceDefaultCargoTargetPaths
            ++ workspaceExplicitCargoTargetPaths
            ++ autoCargoTargetPaths
          ))
        );

      workspaceCargoTargetEntrypoints =
        map (path: workspaceRoot + "/${path}") workspaceCargoTargetRelPaths;

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

      cargoGraphGenerationSrc = lib.fileset.toSource {
        root = workspaceRoot;
        fileset = lib.fileset.unions (
          workspaceDependencyManifestFiles
          ++ workspaceDependencyBuildScripts
          ++ workspaceCargoTargetEntrypoints
          ++ [
            ./crate-hashes.json
            ./.config/hakari.toml
          ]
        );
      };

      dummyCargoTarget = pkgs.writeText "crane-dummy-cargo-target.rs" ''
        #![allow(clippy::all)]
        #![allow(dead_code)]

        pub fn main() {}
      '';

      workspaceDependencyNamesFor = package:
        workspaceGraph.normal_dependencies.${package} or [];

      workspacePackageExtraFilesets = {
        "gammaloop-api" = [
          ./assets
        ];
        gammalooprs = [
          ./assets
        ];
        "gammaloop-integration-tests" = [
          ./assets
          ./examples/cli
        ];
      };

      workspacePackageSrcFor = package:
        lib.fileset.toSource {
          root = workspaceRoot;
          fileset = lib.fileset.unions (
            workspaceDependencyManifestFiles
            ++ workspaceDependencyBuildScripts
            ++ [
              (workspaceRoot + "/${workspaceMemberPackageDirs.${package}}")
              ./tests/resources/fjcore
            ]
            ++ (workspacePackageExtraFilesets.${package} or [])
          );
        };

      workspaceDummyCargoTargetPathsFor = package:
        lib.filter (
          path: !(lib.hasPrefix "${workspaceMemberPackageDirs.${package}}/" path)
        )
        workspaceCargoTargetRelPaths;

      workspaceDummyCargoTargetsScriptFor = package:
        lib.concatMapStringsSep "\n" (path: ''
          if [ ! -e ${lib.escapeShellArg path} ]; then
            install -D -m 0644 ${dummyCargoTarget} ${lib.escapeShellArg path}
          fi
        '')
        (workspaceDummyCargoTargetPathsFor package);

      workspaceAllDummyCargoTargetsScript =
        lib.concatMapStringsSep "\n" (path: ''
          install -D -m 0644 ${dummyCargoTarget} ${lib.escapeShellArg path}
        '')
        workspaceCargoTargetRelPaths;

      packageCargoFeatures = package:
        (workspaceManifestFor workspaceMemberPackageDirs.${package}).features or {};

      packageFeatureIf = package: feature:
        lib.optional (builtins.hasAttr feature (packageCargoFeatures package)) feature;

      craneCiCommonFeaturesFor = package:
        lib.concatMap (packageFeatureIf package) ["bincode" "serde"];

      craneCiExtraFeatureSets = {
        "gammaloop-tracing-filter" = ["clap" "symbolica"];
        idenso = ["reference-cases"];
        linnet = ["symbolica"];
        spenso = ["shadowing"];
        "spenso-macros" = ["shadowing"];
      };

      craneCiFeaturesFor = package:
        sortedUnique (craneCiCommonFeaturesFor package ++ (craneCiExtraFeatureSets.${package} or []));

      craneTestExtraFeatureSets = {
        "gammaloop-integration-tests" = ["python-api-tests"];
      };

      craneTestFeaturesFor = package:
        sortedUnique (craneCiFeaturesFor package ++ (craneTestExtraFeatureSets.${package} or []));

      cargoFeatureArgs = features:
        lib.optionalString (features != []) "--features ${lib.escapeShellArg (lib.concatStringsSep "," features)}";

      cargoPackageArgsFor = package: features:
        lib.concatStringsSep " " (
          [
            "--locked"
            "-p ${lib.escapeShellArg package}"
          ]
          ++ lib.optional (features != []) (cargoFeatureArgs features)
        );

      cargoQualifiedFeatureArgsFor = packages: featuresFor:
        let
          features =
            lib.concatMap (
              package:
                map (feature: "${package}/${feature}") (featuresFor package)
            )
            packages;
        in
          cargoFeatureArgs (sortedUnique features);

      craneWorkspacePrebuildFeatureArgs =
        cargoQualifiedFeatureArgsFor workspaceMemberPackages craneTestFeaturesFor;

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

      gammaloop-cli = craneLib.buildPackage (ciArgs
        // {
          cargoArtifacts = cranePackageBuildArtifacts."gammaloop-api";
          doNotLinkInheritedArtifacts = true;
          pname = "gammaloop-api";
          src = workspacePackageSrcFor "gammaloop-api";
          cargoExtraArgs = cargoPackageArgsFor "gammaloop-api" [];
          doCheck = false;
          postPatch = workspaceDummyCargoTargetsScriptFor "gammaloop-api";
        });
      gammaloop-python-lib = craneLib.buildPackage (ciArgs
        // {
          cargoArtifacts = cranePythonBuildArtifacts;
          doNotLinkInheritedArtifacts = true;
          pname = "gammaloop-api-python";
          src = workspacePackageSrcFor "gammaloop-api";
          cargoExtraArgs = "--locked -p gammaloop-api --no-default-features --features python_abi,pyo3-extension-module";
          doCheck = false;
          postPatch = workspaceDummyCargoTargetsScriptFor "gammaloop-api";
        });
      gammaloop-python-lib-output = lib.getLib gammaloop-python-lib;
      pythonSitePackages = "${pkgs.python313.sitePackages}";
      gammaloop-python-module = pkgs.runCommand "gammaloop-python-module" {} ''
        mkdir -p "$out/${pythonSitePackages}/gammaloop"
        cp ${./crates/gammaloop-api/python/gammaloop/__init__.py} \
          "$out/${pythonSitePackages}/gammaloop/__init__.py"

        extension="$(
          find ${gammaloop-python-lib-output} -type f \
            \( -name 'libgammaloop_api*.so' -o -name 'gammaloop_api*.so' -o -name 'libgammaloop_api*.dylib' -o -name 'gammaloop_api*.dylib' \) \
            | sort \
            | head -n 1
        )"
        if [ -z "$extension" ]; then
          echo "Could not find Crane-built gammaloop-api Python extension in ${gammaloop-python-lib-output}" >&2
          exit 1
        fi
        cp "$extension" "$out/${pythonSitePackages}/gammaloop/_gammaloop.so"
      '';
      clinnet-cli = craneLib.buildPackage (ciArgs
        // {
          cargoArtifacts = cranePackageBuildArtifacts.clinnet;
          doNotLinkInheritedArtifacts = true;
          pname = "clinnet";
          inherit (clinnetMeta) version;
          src = workspacePackageSrcFor "clinnet";
          cargoExtraArgs = cargoPackageArgsFor "clinnet" [];
          doCheck = false;
          postPatch = workspaceDummyCargoTargetsScriptFor "clinnet";
        });

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
          cargoExtraArgs = "--locked --workspace ${craneWorkspacePrebuildFeatureArgs}";
          extraDummyScript =
            lib.concatMapStringsSep "\n" (path: ''
              install -D -m 0644 ${dummyCargoTarget} "$out/${path}"
            '')
            autoCargoTargetPaths;
        });

      guppyWorkspaceGraphNativeBuildInputs = [
        ciToolchain
        pkgs.cargo-guppy
        pkgs.coreutils
        pkgs.gawk
        pkgs.gnugrep
        pkgs.jq
      ];

      guppyWorkspaceGraphJson = pkgs.runCommand "gammaloop-guppy-workspace-graph.json" {
        nativeBuildInputs = guppyWorkspaceGraphNativeBuildInputs;
      } ''
        set -euo pipefail

        cp -R ${cargoGraphGenerationSrc} source
        chmod -R u+w source
        cd source
        mkdir -p .cargo .cargo-home
        cp ${cargoVendorDir}/config.toml .cargo/config.toml
        export CARGO_HOME="$PWD/.cargo-home"

        manifest_path="Cargo.toml"
        tmp="$(mktemp -d)"
        trap 'rm -rf "$tmp"' EXIT

        cargo metadata \
          --format-version 1 \
          --no-deps \
          --manifest-path "$manifest_path" \
          > "$tmp/metadata.json"

        jq '
          def workspace_deps($path_map; $kinds):
            [
              .dependencies[]
              | select((.kind // "normal") as $kind | $kinds | index($kind))
              | select(.path != null)
              | $path_map[.path] // empty
            ]
            | unique
            | sort;

          . as $metadata
          | [
              $metadata.packages[]
              | select(.id as $id | $metadata.workspace_members | index($id))
            ]
            | sort_by(.name) as $packages
          | ($packages | map({key: (.manifest_path | rtrimstr("/Cargo.toml")), value: .name}) | from_entries) as $path_map
          | {
              workspace_root: ".",
              packages: ($packages | map(.name)),
              package_dirs: (
                $packages
                | map({
                    key: .name,
                    value: (.manifest_path | ltrimstr($metadata.workspace_root + "/") | rtrimstr("/Cargo.toml"))
                  })
                | from_entries
              ),
              normal_dependencies: (
                $packages
                | map({key: .name, value: workspace_deps($path_map; ["normal", "build"])})
                | from_entries
              ),
              test_dependencies: (
                $packages
                | map({key: .name, value: workspace_deps($path_map; ["normal", "build", "dev"])})
                | from_entries
              )
            }
        ' "$tmp/metadata.json" > "$tmp/base.json"

        jq -r '.packages[]' "$tmp/base.json" > "$tmp/workspace-packages"

        guppy_names() {
          local package="$1"
          local kind="$2"
          local include_dev="$3"
          if [[ "$include_dev" == "true" ]]; then
            cargo-guppy resolve-cargo \
              --manifest-path "$manifest_path" \
              --resolver-version v2 \
              --target-platform any \
              --host-platform any \
              -k "$kind" \
              -p "$package" \
              --include-dev
          else
            cargo-guppy resolve-cargo \
              --manifest-path "$manifest_path" \
              --resolver-version v2 \
              --target-platform any \
              --host-platform any \
              -k "$kind" \
              -p "$package"
          fi | awk '{ print $1 }' | sort -u
        }

        workspace_guppy_list() {
          local package="$1"
          local include_dev="$2"
          guppy_names "$package" workspace "$include_dev" \
            | comm -12 - "$tmp/workspace-packages" \
            | awk -v self="$package" '$0 != self' \
            | jq -R . \
            | jq -s 'sort'
        }

        write_dependency_object() {
          local include_dev="$1"
          local output="$2"
          local fragments="$tmp/$(basename "$output").fragments"
          : > "$fragments"
          while IFS= read -r package; do
            deps="$(workspace_guppy_list "$package" "$include_dev")"
            jq -n --arg package "$package" --argjson deps "$deps" '{($package): $deps}' >> "$fragments"
          done < "$tmp/workspace-packages"
          jq -s 'add' "$fragments" > "$output"
        }

        write_symbolica_packages() {
          local include_dev="$1"
          local output="$2"
          local names="$tmp/$(basename "$output").names"
          : > "$names"
          while IFS= read -r package; do
            if guppy_names "$package" all "$include_dev" | grep -Eq '^(symbolica|numerica|graphica)$'; then
              printf '%s\n' "$package" >> "$names"
            fi
          done < "$tmp/workspace-packages"
          jq -R . "$names" | jq -s 'sort' > "$output"
        }

        write_dependency_object false "$tmp/resolved-normal.json"
        write_dependency_object true "$tmp/resolved-test.json"
        write_symbolica_packages false "$tmp/symbolica-normal.json"
        write_symbolica_packages true "$tmp/symbolica-test.json"

        jq -S \
          --slurpfile resolvedNormal "$tmp/resolved-normal.json" \
          --slurpfile resolvedTest "$tmp/resolved-test.json" \
          --slurpfile symbolicaNormal "$tmp/symbolica-normal.json" \
          --slurpfile symbolicaTest "$tmp/symbolica-test.json" \
          '. + {
            resolved_normal_dependencies: $resolvedNormal[0],
            resolved_test_dependencies: $resolvedTest[0],
            symbolica_normal_packages: $symbolicaNormal[0],
            symbolica_test_packages: $symbolicaTest[0]
          }' \
          "$tmp/base.json" \
          > "$out"
      '';

      guppyWorkspaceGraphCheck = pkgs.runCommand "gammaloop-guppy-workspace-graph-check" {
        nativeBuildInputs = [
          pkgs.diffutils
        ];
      } ''
        diff -u ${./nix/ci-workspace-graph.json} ${guppyWorkspaceGraphJson}
        mkdir -p "$out"
      '';

      mergeCargoArtifacts = name: artifacts: let
        artifactList = lib.filter (artifact: artifact != null) artifacts;
      in
        if builtins.length artifactList == 1
        then builtins.head artifactList
        else
          pkgs.runCommand name {
            nativeBuildInputs = [pkgs.rsync pkgs.zstd pkgs.gnutar];
          } ''
            mkdir -p "$out/target"

            unpack_artifact() {
              local artifact="$1"

              if [ -d "$artifact" ] && [ -f "$artifact/target.tar.zst" ]; then
                artifact="$artifact/target.tar.zst"
              elif [ -d "$artifact" ] && [ -d "$artifact/target" ]; then
                artifact="$artifact/target"
              fi

              if [ -f "$artifact" ]; then
                if [ -f "$artifact.prev" ]; then
                  unpack_artifact "$(realpath "$artifact.prev")"
                fi
                zstd -d "$artifact" --stdout | tar -x -C "$out/target"
              elif [ -d "$artifact" ]; then
                rsync -a "$artifact/" "$out/target/"
              else
                echo "unsupported cargo artifact path: $artifact" >&2
                exit 1
              fi
            }

            ${lib.concatMapStringsSep "\n" (artifact: "unpack_artifact ${artifact}") artifactList}
          '';

      mergeCargoArtifactsOrNull = name: artifacts: let
        artifactList = lib.filter (artifact: artifact != null) artifacts;
      in
        if artifactList == []
        then null
        else mergeCargoArtifacts name artifactList;

      # Crane's buildDepsOnly deliberately clears cargoArtifacts. This local
      # variant keeps the same dummy-source behavior while allowing per-crate
      # dependency archives to build incrementally from earlier archives.
      buildDepsOnlyWithArtifacts = {
        cargoBuildCommand ? "cargoWithProfile build",
        cargoCheckCommand ? "cargoWithProfile check",
        cargoExtraArgs ? "--locked",
        cargoTestCommand ? "cargoWithProfile test",
        cargoTestExtraArgs ? "--no-run",
        ...
      } @ args: let
        cleanedArgs = builtins.removeAttrs args [
          "cargoBuildCommand"
          "cargoCheckCommand"
          "cargoCheckExtraArgs"
          "cargoExtraArgs"
          "cargoTestCommand"
          "cargoTestExtraArgs"
          "dummySrc"
          "outputHashes"
          "outputs"
        ];
        doCheck = args.doCheck or true;
        cargoCheckExtraArgs = args.cargoCheckExtraArgs or (if doCheck then "--all-targets" else "");
        dummySrc =
          if args ? dummySrc
          then
            lib.warnIf (
              args ? src && args.src != null
            ) "buildDepsOnlyWithArtifacts will ignore `src` when `dummySrc` is specified" args.dummySrc
          else
            craneLib.mkDummySrc args;
        argsMaybeDummySrcOverride =
          if args ? dummySrc
          then args // {src = args.dummySrc;}
          else args;
        crateName = craneLib.crateNameFromCargoToml argsMaybeDummySrcOverride;
      in
        craneLib.mkCargoDerivation (
          cleanedArgs
          // {
            inherit doCheck;

            src = dummySrc;
            pnameSuffix = "-deps";
            pname = args.pname or crateName.pname;
            version = args.version or crateName.version;

            cargoArtifacts = args.cargoArtifacts or null;
            doNotLinkInheritedArtifacts = args.doNotLinkInheritedArtifacts or true;
            cargoVendorDir = args.cargoVendorDir or (craneLib.vendorCargoDeps argsMaybeDummySrcOverride);

            env = (args.env or {}) // {
              CRANE_BUILD_DEPS_ONLY = ((args.env or {}).CRANE_BUILD_DEPS_ONLY or 1);
            };

            buildPhaseCargoCommand =
              args.buildPhaseCargoCommand or ''
                ${cargoCheckCommand} ${cargoExtraArgs} ${cargoCheckExtraArgs}
                ${cargoBuildCommand} ${cargoExtraArgs}
              '';

            checkPhaseCargoCommand =
              args.checkPhaseCargoCommand or ''
                ${cargoTestCommand} ${cargoExtraArgs} ${cargoTestExtraArgs}
              '';

            doInstallCargoArtifacts = true;
          }
        );

      rootPackageDependencyArtifactsFor = package:
        lib.optional (
          package != workspaceHackPackage
          && (workspaceDependencyNamesFor package) == []
          && builtins.elem package workspaceGraph.symbolica_normal_packages
        )
        cargoArtifacts;

      cranePackageDependencyArtifacts = lib.fix (self:
        lib.genAttrs workspacePackages (package:
          if package == workspaceHackPackage
          then cargoArtifacts
          else
            buildDepsOnlyWithArtifacts (ciArgs
              // {
                cargoArtifacts = mergeCargoArtifactsOrNull "gammaloop-crate-${package}-deps-inputs" (
                  rootPackageDependencyArtifactsFor package
                  ++ map (dependency: self.${dependency}) (workspaceDependencyNamesFor package)
                );
                pname = "gammaloop-crate-${package}";
                src = workspacePackageSrcFor package;
                cargoExtraArgs = cargoPackageArgsFor package (craneCiFeaturesFor package);
                doCheck = false;
                extraDummyScript = workspaceAllDummyCargoTargetsScript;
              })));

      cranePackageBuildArtifacts = lib.fix (self:
        lib.genAttrs workspacePackages (package:
          if package == workspaceHackPackage
          then cargoArtifacts
          else
            craneLib.cargoBuild (ciArgs
              // {
                cargoArtifacts = mergeCargoArtifacts "gammaloop-crate-${package}-deps" (
                  [cranePackageDependencyArtifacts.${package}]
                  ++ map (dependency: self.${dependency}) (workspaceDependencyNamesFor package)
                );
                doNotLinkInheritedArtifacts = true;
                pname = "gammaloop-crate-${package}";
                src = workspacePackageSrcFor package;
                cargoExtraArgs = cargoPackageArgsFor package (craneCiFeaturesFor package);
                postPatch = workspaceDummyCargoTargetsScriptFor package;
              })));

      craneTestDependencyArtifacts = lib.fix (self:
        lib.genAttrs workspacePackages (package:
          if package == workspaceHackPackage
          then cargoArtifacts
          else
            buildDepsOnlyWithArtifacts (ciArgs
              // {
                cargoArtifacts = mergeCargoArtifactsOrNull "gammaloop-crate-test-binaries-${package}-deps-inputs" (
                  [cranePackageDependencyArtifacts.${package}]
                  ++ map (dependency: self.${dependency}) (workspaceDependencyNamesFor package)
                );
                pname = "gammaloop-crate-test-binaries-${package}";
                src = workspacePackageSrcFor package;
                cargoExtraArgs = cargoPackageArgsFor package (craneTestFeaturesFor package);
                cargoTestExtraArgs = "--no-run";
                extraDummyScript = workspaceAllDummyCargoTargetsScript;
              })));

      craneTestBinaryArtifacts = lib.genAttrs workspacePackages (package:
        if package == workspaceHackPackage
        then cargoArtifacts
        else
          craneLib.mkCargoDerivation (ciArgs
            // {
              cargoArtifacts = mergeCargoArtifacts "gammaloop-crate-test-binaries-${package}-deps" (
                [craneTestDependencyArtifacts.${package} cranePackageBuildArtifacts.${package}]
                ++ map (dependency: cranePackageBuildArtifacts.${dependency}) (workspaceDependencyNamesFor package)
              );
              doNotLinkInheritedArtifacts = true;
              pname = "gammaloop-crate-test-binaries-${package}";
              src = workspacePackageSrcFor package;
              buildPhaseCargoCommand = "cargoWithProfile test --no-run ${cargoPackageArgsFor package (craneTestFeaturesFor package)}";
              checkPhaseCargoCommand = "";
              doCheck = false;
              postPatch = workspaceDummyCargoTargetsScriptFor package;
            }));

      cranePythonBuildArtifacts = craneLib.cargoBuild (ciArgs
        // {
          cargoArtifacts = mergeCargoArtifacts "gammaloop-python-deps" (
            [cargoArtifacts]
            ++ map (dependency: cranePackageBuildArtifacts.${dependency}) (workspaceDependencyNamesFor "gammaloop-api")
          );
          doNotLinkInheritedArtifacts = true;
          pname = "gammaloop-api-python-build";
          src = workspacePackageSrcFor "gammaloop-api";
          cargoExtraArgs = "--locked -p gammaloop-api --no-default-features --features python_abi,pyo3-extension-module";
          postPatch = workspaceDummyCargoTargetsScriptFor "gammaloop-api";
        });

      cranePackageOutputs = lib.listToAttrs (map (package: {
          name = "crate-${package}";
          value = cranePackageBuildArtifacts.${package};
        })
        workspacePackages);

      cranePackageDependencyOutputs = lib.listToAttrs (map (package: {
          name = "crate-deps-${package}";
          value = cranePackageDependencyArtifacts.${package};
        })
        workspacePackages);

      craneTestBinaryPackageOutputs = lib.listToAttrs (map (package: {
          name = "crate-test-binaries-${package}";
          value = craneTestBinaryArtifacts.${package};
        })
        workspacePackages);

      craneTestBinaryDependencyOutputs = lib.listToAttrs (map (package: {
          name = "crate-test-binaries-deps-${package}";
          value = craneTestDependencyArtifacts.${package};
        })
        workspacePackages);

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
          extraFeatures."gammaloop-integration-tests" = ["python-api-tests"];
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
      nextestCoverageIgnoredPackages = [
        "gammaloop-workspace-hack"
      ];
      nextestCoveredWorkspacePackages = sortedUnique (lib.subtractLists nextestCoverageIgnoredPackages workspaceMemberPackages);
      missingNextestPackages = lib.subtractLists nextestSplitPackages nextestCoveredWorkspacePackages;
      extraNextestPackages = lib.subtractLists nextestCoveredWorkspacePackages nextestSplitPackages;

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

      nextestFeatureArgsFor = target:
        cargoQualifiedFeatureArgsFor target.packages (
          package:
            sortedUnique (
              craneCiFeaturesFor package
              ++ lib.optionals ((target ? extraFeatures) && builtins.hasAttr package target.extraFeatures) target.extraFeatures.${package}
            )
        );

      nextestCargoArgsFor = target:
        lib.concatStringsSep " " (
          [
            "--locked"
          ]
          ++ map (package: "-p ${lib.escapeShellArg package}") target.packages
          ++ lib.optional (nextestFeatureArgsFor target != "") (nextestFeatureArgsFor target)
        );

      nextestBinarySetFor = target:
        mergeCargoArtifacts "gammaloop-nextest-binaries-${target.name}"
        (map (package: craneTestBinaryArtifacts.${package}) target.packages);

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
        craneLib.cargoNextest (ciArgs
          // {
          pname = "gammaloop-nextest-${target.name}";
          src = nextestSrcFor target;
          cargoArtifacts = nextestBinarySetForTarget target;
          doNotLinkInheritedArtifacts = true;
          cargoExtraArgs = nextestCargoArgsFor target;
          cargoNextestExtraArgs = "${nextestFilterFor target} ${nextestBaseExtraArgs}";
          nativeBuildInputs =
            (ciArgs.nativeBuildInputs or [])
            ++ [pkgs.form]
            ++ lib.optionals (nextestUsesIntegrationTests target) [nextestPython];
          RUST_BACKTRACE = "1";
          RUST_LIB_BACKTRACE = "1";
          SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
          preCheck = ''
            ${licensePreCheck}
            # Nextest runs with the workspace root as cwd, while some insta
            # snapshots are stored under each crate src. Mirror those snapshot
            # directories into the workspace-level src tree in the disposable
            # Nix build directory so insta can find them.
            if [ -d crates ]; then
              while IFS= read -r snapshots; do
                rel="''${snapshots#crates/*/src/}"
                mkdir -p "src/$rel"
                cp -R "$snapshots/." "src/$rel/"
              done < <(find crates -type d -name snapshots | sort)
            fi
          '';
          postCheck = ''
            ${nextestFailureSummary}/bin/nextest-failure-summary ${lib.escapeShellArg nextestJunitPath} || true
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

          gammaloop-guppy-workspace-graph = guppyWorkspaceGraphCheck;

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
          "ci-workspace-graph-json" = guppyWorkspaceGraphJson;
          inherit linnest-wasm linnestWasmCargoArtifacts;
          "crane-ci-prebuild" = cargoArtifacts;
          inherit cargoArtifacts workspaceBuildArtifacts;
          "nix-ci-passed" = nixCiPassed;
        }
        // cranePackageDependencyOutputs
        // cranePackageOutputs
        // craneTestBinaryDependencyOutputs
        // craneTestBinaryPackageOutputs
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
              clinnet-cli
            ]
            ++ lib.optionals (!pkgs.stdenv.isDarwin) [
              valgrind
            ];
        });
    });
}
