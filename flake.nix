{
  description = "Gammaloop";

  # Substitute what NixCI has already built instead of building it again.
  # Reading from this cache needs a token in your netrc as well, and Nix asks
  # before it trusts these settings; see CONTRIBUTING.typ.
  nixConfig = {
    extra-substituters = ["https://cache.nix-ci.com"];
    extra-trusted-public-keys = ["nix-ci:g3xV5BDTLtIBZr/A00IU1x0EtKKlb7YLgBN2SgYgM6A="];
  };

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

    # Keep generated documentation assets stable across general Nixpkgs updates.
    nixpkgs-docs.url = "github:NixOS/nixpkgs/1559d3daa3ecc813a650b79375ea61b6741b8746";

    crane = {
      url = "github:ipetkov/crane";
    };

    # Refresh deliberately with `just ci-cache-base REVISION` after a green run.
    ci-cache-base = {
      url = "github:alphal00p/gammaloop/9935a2b52d047daf5bbe6616c798fa07e364da6a";
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
    nixpkgs-docs,
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
      docsPkgs = nixpkgs-docs.legacyPackages.${system};
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
        inherit self pkgs docsPkgs craneLib wasmCraneLib ciToolchain wasmTarget system nixCiArtifactBarrier;
        workspaceRoot = ./.;
        incrementalBaselineRoot = /. + builtins.unsafeDiscardStringContext self.inputs.ci-cache-base.outPath;
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
        fileset = lib.fileset.unions (
          [
            cargoSources
            nonCargoBuildSources
            snapshotSources
            ./tests
            ./examples/cli
          ]
          ++ (workspacePackageTestCompileTimeExtraFilesetsForSourcePackages workspaceMemberPackages)
        );
      };

      workspaceNonIntegrationTestSrc = lib.fileset.toSource {
        root = workspaceRoot;
        fileset = lib.fileset.unions (
          [
            nonIntegrationCargoSources
            nonCargoBuildSources
            snapshotSources
            ./tests/resources
            ./examples/cli
          ]
          ++ (workspacePackageTestCompileTimeExtraFilesetsForSourcePackages workspaceMemberPackages)
        );
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
              "${member}/src/bin"
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

      workspacePackageBuildScriptsForSourcePackages = sourcePackages:
        lib.filter builtins.pathExists (
          map (package: workspaceRoot + "/${workspaceMemberPackageDirs.${package}}/build.rs") sourcePackages
        );

      workspaceDependencySrc = lib.fileset.toSource {
        root = workspaceRoot;
        fileset = lib.fileset.unions workspaceDependencyManifestFiles;
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

      dummyProcMacroCargoTarget = pkgs.writeText "crane-dummy-proc-macro-cargo-target.rs" ''
        #![allow(clippy::all)]
        #![allow(dead_code)]
      '';

      normalizeWorkspaceHackBuildScriptTimestampScriptFor = prefix: ''
        if [ -e ${prefix}crates/${workspaceHackPackage}/Cargo.toml ]; then
          touch -d @0 ${prefix}crates/${workspaceHackPackage}/Cargo.toml
        fi
        if [ -e ${prefix}crates/${workspaceHackPackage}/src/lib.rs ]; then
          touch -d @0 ${prefix}crates/${workspaceHackPackage}/src/lib.rs
        fi
        if [ -e ${prefix}crates/${workspaceHackPackage}/build.rs ]; then
          touch -d @1 ${prefix}crates/${workspaceHackPackage}/build.rs
        fi
      '';

      normalizeWorkspaceHackBuildScriptTimestampScript =
        normalizeWorkspaceHackBuildScriptTimestampScriptFor "";

      normalizeWorkspaceHackBuildScriptTimestampInDummySrcScript =
        normalizeWorkspaceHackBuildScriptTimestampScriptFor "$out/";

      workspaceDependencyNamesFor = package:
        workspaceGraph.normal_dependencies.${package} or [];

      workspaceResolvedDependencyNamesFor = package:
        workspaceGraph.resolved_normal_dependencies.${package} or (workspaceDependencyNamesFor package);

      workspaceResolvedTestDependencyNamesFor = package:
        workspaceGraph.resolved_test_dependencies.${package} or (workspaceResolvedDependencyNamesFor package);

      workspaceTestDependencyNamesFor = package:
        workspaceGraph.test_dependencies.${package} or (workspaceDependencyNamesFor package);

      workspaceDependencyClosureFor = dependencyNamesFor: package:
        sortedUnique (map (entry: entry.key) (builtins.genericClosure {
          startSet = [{key = package;}];
          operator = entry: map (dependency: {key = dependency;}) (dependencyNamesFor entry.key);
        }));

      workspaceTestDependencyClosureFor =
        workspaceDependencyClosureFor workspaceTestDependencyNamesFor;

      workspaceTestComponentMembersFor = package:
        sortedUnique (lib.filter (
            other:
              builtins.elem other (workspaceTestDependencyClosureFor package)
              && builtins.elem package (workspaceTestDependencyClosureFor other)
          )
          workspaceMemberPackages);

      workspaceTestComponentRepresentativeFor = package:
        builtins.head (workspaceTestComponentMembersFor package);

      workspaceTestComponentRepresentatives =
        sortedUnique (map workspaceTestComponentRepresentativeFor workspaceMemberPackages);

      workspaceTestDependencyComponentRepresentatives =
        lib.filter (representative: representative != workspaceHackPackage) workspaceTestComponentRepresentatives;

      workspaceTestComponentMembers =
        lib.listToAttrs (map (representative: {
            name = representative;
            value = workspaceTestComponentMembersFor representative;
          })
          workspaceTestComponentRepresentatives);

      workspaceTestComponentDependencyRepresentativesFor = representative:
        sortedUnique (lib.filter (dependencyRepresentative: dependencyRepresentative != representative) (
          map workspaceTestComponentRepresentativeFor (
            lib.concatMap workspaceTestDependencyNamesFor workspaceTestComponentMembers.${representative}
          )
        ));

      workspaceSourcePackageNamesFor = package: dependencies:
        sortedUnique ([package] ++ dependencies);

      workspaceNormalSourcePackageNamesFor =
        workspaceDependencyClosureFor workspaceResolvedDependencyNamesFor;

      workspaceTestSourcePackageNamesFor = package:
        sortedUnique (
          lib.concatMap workspaceNormalSourcePackageNamesFor (
            workspaceSourcePackageNamesFor package (workspaceResolvedTestDependencyNamesFor package)
          )
        );

      workspacePackageProductionExtraSourceRoots = {
        "gammaloop-api" = [
          "assets/embedded"
          "assets/models"
          "crates/kurvst/typst/kurvst.wasm"
          "crates/kurvst/typst/src"
          "crates/kurvst/typst/typst.toml"
          "crates/linnest/typst/linnest.wasm"
          "crates/linnest/typst/src"
          "crates/linnest/typst/typst.toml"
        ];
        gammalooprs = [
          "assets/models/json"
        ];
        "gammaloop-integration-tests" = [
          "tests/resources/fjcore"
        ];
        "feynkit-py" = [
          "crates/kurvst/typst/kurvst.wasm"
          "crates/kurvst/typst/src"
          "crates/linnest/typst/linnest.wasm"
          "crates/linnest/typst/src"
        ];
        clinnet = [
          "crates/clinnet/templates"
          "crates/kurvst/typst/kurvst.wasm"
          "crates/kurvst/typst/src"
          "crates/kurvst/typst/typst.toml"
          "crates/linnest/typst/linnest.wasm"
          "crates/linnest/typst/src"
          "crates/linnest/typst/typst.toml"
        ];
        vakint = [
          "crates/vakint/form_src"
          "crates/vakint/templates"
        ];
      };

      workspacePackageTestCompileTimeExtraSourceRoots = {
        "feynkit-cff" = [
          "crates/feynkit-model/tests/fixtures"
        ];
        "feynkit-generator" = [
          "crates/feynkit-model/tests/fixtures"
        ];
        "feynkit-model" = [
          "crates/feynkit-model/tests/fixtures"
        ];
        "feynkit-tensor" = [
          "crates/feynkit-model/tests/fixtures"
        ];
        "feynkit-py" = [
          "crates/feynkit-model/tests/fixtures"
          "crates/feynkit-py/python/symbolica/community/feynkit/__init__.py"
          "crates/feynkit-py/tests/fixtures"
        ];
        gammalooprs = [
          "tests/resources/graphs/scalar/dod2_bubble.dot"
        ];
      };

      workspacePackageRuntimeTestExtraSourceRoots = {
        "gammaloop-api" = [
          "tests/resources/graphs/scalar_bubble.dot"
        ];
        "gammaloop-tracing-filter" = [
          "tests/resources/run_cards"
        ];
        gammalooprs = [
          "tests/resources/graphs/uv_tests/rqft_a_3l_no_ghost.dot"
          "tests/resources/graphs/uv_tests/rqft_ghG_3l.dot"
        ];
        "gammaloop-integration-tests" = [
          "assets/plot_approach_result.py"
          "examples/api"
          "examples/cli"
          "tests/resources"
        ];
      };

      workspacePackageOwnTestSourceRoots = {
        gammalooprs = [
          "crates/gammalooprs/src/graph/parse/tests.rs"
          "crates/gammalooprs/src/numerator/spensotests.rs"
          "crates/gammalooprs/src/numerator/tests.rs"
          "crates/gammalooprs/src/utils/test_utils.rs"
          "crates/gammalooprs/src/uv/tests.rs"
        ];
        idenso = [
          "crates/idenso/src/color/test"
          "crates/idenso/src/dirac/test"
          "crates/idenso/src/shorthands/schoonschip/test"
          "crates/idenso/src/tensor/tests"
          "crates/idenso/src/test_support.rs"
        ];
        linnet = [
          "crates/linnet/src/half_edge/involution/test.rs"
          "crates/linnet/src/half_edge/nodestore/test.rs"
          "crates/linnet/src/half_edge/test_graphs.rs"
          "crates/linnet/src/half_edge/tests.rs"
          "crates/linnet/src/union_find/test.rs"
        ];
        spenso = [
          "crates/spenso/src/data/tests.rs"
          "crates/spenso/src/iterators/tests"
          "crates/spenso/src/network/parsing/test.rs"
          "crates/spenso/src/network/shadowing_tests.rs"
          "crates/spenso/src/network/test.rs"
          "crates/spenso/src/network/tests.rs"
          "crates/spenso/src/shadowing/tests.rs"
          "crates/spenso/src/tests.rs"
        ];
      };

      workspacePackageExtraSourceRootsForSourcePackages = sourcePackages:
        sortedUnique (lib.concatMap (sourcePackage: workspacePackageProductionExtraSourceRoots.${sourcePackage} or []) sourcePackages);

      workspacePackageTestCompileTimeExtraSourceRootsForSourcePackages = sourcePackages:
        sortedUnique (lib.concatMap (sourcePackage: workspacePackageTestCompileTimeExtraSourceRoots.${sourcePackage} or []) sourcePackages);

      workspacePackageRuntimeTestExtraSourceRootsForSourcePackages = sourcePackages:
        sortedUnique (lib.concatMap (sourcePackage: workspacePackageRuntimeTestExtraSourceRoots.${sourcePackage} or []) sourcePackages);

      workspacePackageExtraFilesetsForSourcePackages = sourcePackages:
        map (sourceRoot: workspaceRoot + "/${sourceRoot}") (workspacePackageExtraSourceRootsForSourcePackages sourcePackages);

      workspacePackageTestCompileTimeExtraFilesetsForSourcePackages = sourcePackages:
        map (sourceRoot: workspaceRoot + "/${sourceRoot}") (workspacePackageTestCompileTimeExtraSourceRootsForSourcePackages sourcePackages);

      workspacePackageRuntimeTestExtraFilesetsForSourcePackages = sourcePackages:
        map (sourceRoot: workspaceRoot + "/${sourceRoot}") (workspacePackageRuntimeTestExtraSourceRootsForSourcePackages sourcePackages);

      workspacePackageOwnTestSourceRootsForSourcePackages = sourcePackages:
        sortedUnique (lib.concatMap (sourcePackage: workspacePackageOwnTestSourceRoots.${sourcePackage} or []) sourcePackages);

      workspacePackageOwnTestFilesetsForSourcePackages = sourcePackages:
        map (sourceRoot: workspaceRoot + "/${sourceRoot}") (workspacePackageOwnTestSourceRootsForSourcePackages sourcePackages);

      workspacePackageExtraSourceRestoreInDummySrcScriptFor = sourcePackages:
        ''
          ${lib.concatMapStringsSep "\n" (sourceRoot: let
            sourceParentDir = builtins.dirOf sourceRoot;
            source = workspaceRoot + "/${sourceRoot}";
          in ''
            rm -rf "$out"/${lib.escapeShellArg sourceRoot}
            mkdir -p "$out"/${lib.escapeShellArg sourceParentDir}
            cp -R --no-preserve=ownership ${source} "$out"/${lib.escapeShellArg sourceRoot}
            chmod -R u+w "$out"/${lib.escapeShellArg sourceRoot}
          '')
          (workspacePackageExtraSourceRootsForSourcePackages sourcePackages)}
        '';

      workspacePackageSrcForSourcePackages = {
        sourcePackages,
        packageSourcePackages ? [],
        testSourcePackages ? [],
        runtimeTestSourcePackages ? [],
        extraFilesets ? [],
      }: let
        productionSourceFilesets = map (sourcePackage: let
          packageRoot = workspaceRoot + "/${workspaceMemberPackageDirs.${sourcePackage}}";
          libraryTargetPath = workspacePackageLibTargetRelPath sourcePackage;
          testRoots = lib.filter builtins.pathExists (map (root: packageRoot + "/${root}") [
            "benches"
            "examples"
            "tests"
          ]);
          nonLibraryTargetEntrypoints = map (path: workspaceRoot + "/${path}") (lib.filter (
              path:
                lib.hasPrefix "${workspaceMemberPackageDirs.${sourcePackage}}/" path
                && path != libraryTargetPath
            )
            workspaceCargoTargetRelPaths);
          ownTestSources = workspacePackageOwnTestFilesetsForSourcePackages [sourcePackage];
          rustSources = lib.fileset.fileFilter (file: file.hasExt "rs") packageRoot;
        in
          lib.fileset.difference rustSources (lib.fileset.unions (testRoots ++ nonLibraryTargetEntrypoints ++ ownTestSources)))
        sourcePackages;
        packageSourceFilesets = map (sourcePackage: let
          libraryTargetPath = workspacePackageLibTargetRelPath sourcePackage;
        in
          map (path: workspaceRoot + "/${path}") (lib.filter (
              path:
                lib.hasPrefix "${workspaceMemberPackageDirs.${sourcePackage}}/" path
                && path != libraryTargetPath
            )
            workspaceCargoTargetRelPaths))
        packageSourcePackages;
        testSourceFilesets = map (sourcePackage: let
          packageRoot = workspaceRoot + "/${workspaceMemberPackageDirs.${sourcePackage}}";
          testRoots = lib.filter builtins.pathExists (map (root: packageRoot + "/${root}") [
            "benches"
            "examples"
            "tests"
          ]);
        in
          map (testRoot: lib.fileset.fileFilter (file: file.hasExt "rs") testRoot) testRoots)
        testSourcePackages;
        runtimeTestSourceFilesets = map (sourcePackage: let
          packageRoot = workspaceRoot + "/${workspaceMemberPackageDirs.${sourcePackage}}";
        in
          lib.fileset.fileFilter (file: file.hasExt "snap") packageRoot)
        runtimeTestSourcePackages;
      in
        lib.fileset.toSource {
          root = workspaceRoot;
          fileset = lib.fileset.unions (
            workspaceDependencyManifestFiles
            ++ (workspacePackageBuildScriptsForSourcePackages sourcePackages)
            ++ productionSourceFilesets
            ++ lib.concatLists packageSourceFilesets
            ++ lib.concatLists testSourceFilesets
            ++ runtimeTestSourceFilesets
            ++ (workspacePackageExtraFilesetsForSourcePackages sourcePackages)
            ++ (workspacePackageTestCompileTimeExtraFilesetsForSourcePackages testSourcePackages)
            ++ (workspacePackageRuntimeTestExtraFilesetsForSourcePackages runtimeTestSourcePackages)
            ++ (workspacePackageOwnTestFilesetsForSourcePackages testSourcePackages)
            ++ extraFilesets
          );
        };

      workspacePackageSrcFor = package:
        workspacePackageSrcForSourcePackages {
          sourcePackages = workspaceNormalSourcePackageNamesFor package;
          packageSourcePackages = [package];
        };

      workspaceTestPackageSrcFor = package:
        workspacePackageSrcForSourcePackages {
          sourcePackages = workspaceTestSourcePackageNamesFor package;
          packageSourcePackages = [package];
          testSourcePackages = [package];
        };

      workspacePackageIsProcMacro = package: let
        manifest = workspaceManifestFor workspaceMemberPackageDirs.${package};
        libManifest = manifest.lib or {};
        crateTypes = libManifest."crate-type" or [];
      in
        (libManifest."proc-macro" or false)
        || (libManifest.proc_macro or false)
        || builtins.elem "proc-macro" crateTypes
        || builtins.elem "proc_macro" crateTypes;

      workspacePackageHasLibTarget = package: let
        manifest = workspaceManifestFor workspaceMemberPackageDirs.${package};
        packageDir = workspaceMemberPackageDirs.${package};
      in
        manifest ? lib || builtins.pathExists (workspaceRoot + "/${packageDir}/src/lib.rs");

      workspacePackageLibTargetRelPath = package: let
        manifest = workspaceManifestFor workspaceMemberPackageDirs.${package};
      in
        "${workspaceMemberPackageDirs.${package}}/${manifest.lib.path or "src/lib.rs"}";

      workspacePackageForCargoTargetPath = path:
        lib.findFirst (
          package:
            lib.hasPrefix "${workspaceMemberPackageDirs.${package}}/" path
        )
        null
        workspaceMemberPackages;

      workspaceDummyCargoTargetForPath = path: let
        package = workspacePackageForCargoTargetPath path;
      in
        if package != null
        && workspacePackageIsProcMacro package
        && path == workspacePackageLibTargetRelPath package
        then dummyProcMacroCargoTarget
        else dummyCargoTarget;

      workspaceMissingCargoTargetsScript =
        ''
          ${lib.concatMapStringsSep "\n" (path: ''
            if [ ! -e ${lib.escapeShellArg path} ]; then
              install -D -m 0644 ${workspaceDummyCargoTargetForPath path} ${lib.escapeShellArg path}
            fi
          '')
          workspaceCargoTargetRelPaths}

          ${normalizeWorkspaceHackBuildScriptTimestampScript}
        '';

      workspaceAllDummyCargoTargetsScript =
        ''
          ${lib.concatMapStringsSep "\n" (path: ''
            install -D -m 0644 ${workspaceDummyCargoTargetForPath path} "$out"/${lib.escapeShellArg path}
          '')
          workspaceCargoTargetRelPaths}

          ${normalizeWorkspaceHackBuildScriptTimestampInDummySrcScript}
        '';

      workspaceAllDummyCargoTargetsPreservingProcMacrosScript =
        ''
          ${lib.concatMapStringsSep "\n" (path: let
            package = workspacePackageForCargoTargetPath path;
            keepRealProcMacro =
              package != null
              && workspacePackageIsProcMacro package
              && path == workspacePackageLibTargetRelPath package;
            source = workspaceRoot + "/${path}";
          in
            if keepRealProcMacro
            then ''
              install -D -m 0644 ${source} "$out"/${lib.escapeShellArg path}
            ''
            else ''
              install -D -m 0644 ${workspaceDummyCargoTargetForPath path} "$out"/${lib.escapeShellArg path}
            '')
          workspaceCargoTargetRelPaths}

          ${normalizeWorkspaceHackBuildScriptTimestampInDummySrcScript}
        '';

      workspacePackageSourceRestoreInDummySrcScriptFor = sourcePackages:
        ''
          ${lib.concatMapStringsSep "\n" (sourcePackage: let
            packageDir = workspaceMemberPackageDirs.${sourcePackage};
            packageParentDir = builtins.dirOf packageDir;
            sourceRoot = workspacePackageSrcForSourcePackages {
              sourcePackages = [sourcePackage];
            };
            source = sourceRoot + "/${packageDir}";
          in ''
            rm -rf "$out"/${lib.escapeShellArg packageDir}
            mkdir -p "$out"/${lib.escapeShellArg packageParentDir}
            cp -R --no-preserve=ownership ${source} "$out"/${lib.escapeShellArg packageDir}
            chmod -R u+w "$out"/${lib.escapeShellArg packageDir}
          '')
          sourcePackages}

          ${workspacePackageExtraSourceRestoreInDummySrcScriptFor sourcePackages}
        '';

      workspaceDependencyDummyCargoTargetsScriptFor = package: let
        dependencySourcePackages =
          lib.filter (sourcePackage: sourcePackage != package) (workspaceNormalSourcePackageNamesFor package);
        packageTargetPaths =
          lib.filter (
            path:
              lib.hasPrefix "${workspaceMemberPackageDirs.${package}}/" path
          )
          workspaceCargoTargetRelPaths;
      in
        ''
          ${workspacePackageSourceRestoreInDummySrcScriptFor dependencySourcePackages}

          ${lib.concatMapStringsSep "\n" (path: ''
            install -D -m 0644 ${workspaceDummyCargoTargetForPath path} "$out"/${lib.escapeShellArg path}
          '')
          packageTargetPaths}

          ${normalizeWorkspaceHackBuildScriptTimestampInDummySrcScript}
        '';

      packageCargoFeatures = package:
        (workspaceManifestFor workspaceMemberPackageDirs.${package}).features or {};

      packageFeatureIf = package: feature:
        lib.optional (builtins.hasAttr feature (packageCargoFeatures package)) feature;

      craneCiCommonFeaturesFor = package:
        lib.concatMap (packageFeatureIf package) ["bincode" "serde"];

      craneCiExtraFeatureSets = {
        "gammaloop-workspace-hack" = ["symbolica/tracing_max_level_info"];
        "gammaloop-tracing-filter" = ["clap" "symbolica"];
        idenso = ["reference-cases"];
        linnet = ["drawing" "symbolica"];
        spenso = ["shadowing"];
        "spenso-macros" = ["shadowing"];
      };

      workspaceFeatureUnificationExcludedPackages = [
        "feynkit-py"
        "linnet-py"
        "spynso3"
      ];

      workspaceIncomingNormalDependencyFeaturesFor = dependency:
        sortedUnique (lib.concatMap (
            package: workspaceGraph.normal_dependency_features.${package}.${dependency} or []
          ) (lib.subtractLists workspaceFeatureUnificationExcludedPackages workspaceMemberPackages));

      workspaceIncomingTestDependencyFeaturesFor = packages: dependency:
        sortedUnique (lib.concatMap (
            package: workspaceGraph.test_dependency_features.${package}.${dependency} or []
          )
          packages);

      craneCiFeaturesFor = package:
        sortedUnique (
          craneCiCommonFeaturesFor package
          ++ (craneCiExtraFeatureSets.${package} or [])
          ++ (workspaceIncomingNormalDependencyFeaturesFor package)
        );

      craneTestExtraFeatureSets = {
        "spenso-macros" = ["spenso/shadowing"];
      };

      craneTestFeaturesFor = package:
        sortedUnique (craneCiFeaturesFor package ++ (craneTestExtraFeatureSets.${package} or []));

      craneTestContextFeaturesFor = sourcePackages: package:
        sortedUnique (
          craneCiCommonFeaturesFor package
          ++ (craneCiExtraFeatureSets.${package} or [])
          ++ (craneTestExtraFeatureSets.${package} or [])
          ++ (workspaceIncomingTestDependencyFeaturesFor sourcePackages package)
        );

      workspaceTestContextFor = {
        packages,
        extraFeatures ? {},
      }: let
        componentRepresentatives =
          sortedUnique (map workspaceTestComponentRepresentativeFor packages);
        componentPackages =
          sortedUnique (lib.concatMap (representative: workspaceTestComponentMembers.${representative}) componentRepresentatives);
        sourcePackages =
          sortedUnique (lib.concatMap workspaceTestSourcePackageNamesFor componentPackages);
        anchorPackages =
          lib.optionals (
            builtins.any (
              package:
                package != workspaceHackPackage
                && builtins.elem package workspaceGraph.symbolica_test_packages
            )
            sourcePackages
          ) [workspaceHackPackage];
        featurePackages = sortedUnique (sourcePackages ++ anchorPackages);
        features = lib.listToAttrs (map (package: {
            name = package;
            value = sortedUnique (craneTestContextFeaturesFor sourcePackages package ++ (extraFeatures.${package} or []));
          })
          featurePackages);
        resolvedFeatureVector = map (package: {
            inherit package;
            features = features.${package};
          })
          featurePackages;
        usesPythonModule = builtins.any (
          package: builtins.elem "python-api-tests" features.${package}
        ) featurePackages;
        compileEnvironment = {
          inherit (ciArgs) NO_SYMBOLICA_OEM_LICENSE;
          inherit (commonArgs) CC CXX RUSTFLAGS;
          PYO3_PYTHON =
            if usesPythonModule
            then "${nextestPython}/bin/python3"
            else ciArgs.PYO3_PYTHON;
        };
        key = builtins.substring 0 16 (builtins.hashString "sha256" (builtins.toJSON {
          inherit compileEnvironment resolvedFeatureVector;
          profile = ciCargoProfile;
          target = nextestTargetTriple;
        }));
      in {
        inherit
          anchorPackages
          componentPackages
          componentRepresentatives
          extraFeatures
          featurePackages
          features
          key
          resolvedFeatureVector
          sourcePackages
          usesPythonModule
          ;
      };

      workspaceTestDependencyContextsFor = context: let
        dependencyRepresentatives = sortedUnique (lib.concatMap workspaceTestComponentDependencyRepresentativesFor context.componentRepresentatives);
        dependencyContexts = map (representative:
          workspaceTestContextFor {
            packages = workspaceTestComponentMembers.${representative};
            inherit (context) extraFeatures;
          })
        dependencyRepresentatives;
      in
        builtins.attrValues (builtins.removeAttrs (lib.listToAttrs (map (dependencyContext: {
              name = dependencyContext.key;
              value = dependencyContext;
            })
            dependencyContexts)) [context.key]);

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

      cargoQualifiedFeaturesFor = packages: featuresFor:
        let
          features =
            lib.concatMap (
              package:
                map (
                  feature:
                    if lib.hasInfix "/" feature
                    then feature
                    else "${package}/${feature}"
                ) (featuresFor package)
            )
            packages;
        in
          sortedUnique features;

      cargoQualifiedFeatureArgsFor = packages: featuresFor:
        cargoFeatureArgs (cargoQualifiedFeaturesFor packages featuresFor);

      cargoPackageArgsWithFeaturePackagesFor = {
        package,
        featurePackages,
        featuresFor,
        selectFeaturePackages ? true,
      }: let
        anchorPackages =
          lib.optionals (
            package != workspaceHackPackage
            && builtins.elem package workspaceGraph.symbolica_normal_packages
          ) [workspaceHackPackage];
        selectedFeaturePackages =
          if selectFeaturePackages
          then
            lib.filter (
              featurePackage:
                (featuresFor featurePackage) != []
                && !workspacePackageIsProcMacro featurePackage
            )
            featurePackages
          else [];
        selectedPackages = sortedUnique ([package] ++ anchorPackages ++ selectedFeaturePackages);
        featureArgs =
          cargoQualifiedFeatureArgsFor (sortedUnique (featurePackages ++ anchorPackages)) featuresFor;
      in
        lib.concatStringsSep " " (
          [
            "--locked"
          ]
          ++ map (selectedPackage: "-p ${lib.escapeShellArg selectedPackage}") selectedPackages
          ++ lib.optional (featureArgs != "") featureArgs
        );

      cargoPackagesArgsFor = packages: featuresFor: let
        featureArgs = cargoQualifiedFeatureArgsFor packages featuresFor;
      in
        lib.concatStringsSep " " (
          [
            "--locked"
          ]
          ++ map (package: "-p ${lib.escapeShellArg package}") packages
          ++ lib.optional (featureArgs != "") featureArgs
        );

      craneWorkspacePrebuildFeatureArgs =
        cargoQualifiedFeatureArgsFor workspaceMemberPackages craneTestFeaturesFor;
      workspacePrebuildExcludedPackages = workspaceFeatureUnificationExcludedPackages;
      workspacePrebuildDependencyPackages =
        lib.filter (
          package:
            package != workspaceHackPackage
            && !(builtins.elem package workspacePrebuildExcludedPackages)
            && workspacePackageHasLibTarget package
        )
        workspaceMemberPackages;
      workspacePrebuildDependencyPathFor = package: let
        packageDir = workspaceMemberPackageDirs.${package};
      in
        if lib.hasPrefix "crates/" packageDir
        then "../${lib.removePrefix "crates/" packageDir}"
        else "../../${packageDir}";
      workspacePrebuildCargoToml = pkgs.writeText "${workspacePrebuildPackage}-Cargo.toml" ''
        [package]
        name = "${workspacePrebuildPackage}"
        version = "0.1.0"
        edition = "2024"
        publish = false

        [lib]
        path = "src/lib.rs"

        [dependencies]
        ${lib.concatMapStringsSep "\n" (package: let
          features = lib.filter (feature: !lib.hasInfix "/" feature) (craneTestFeaturesFor package);
          featureEntry = lib.optionalString (features != []) ", features = ${builtins.toJSON features}";
        in ''
          ${package} = { path = "${workspacePrebuildDependencyPathFor package}"${featureEntry} }
        '') workspacePrebuildDependencyPackages}
      '';
      workspacePrebuildSourceScript = ''
        install -D -m 0644 ${workspacePrebuildCargoToml} "$out/${workspacePrebuildPackageDir}/Cargo.toml"
        install -D -m 0644 ${dummyCargoTarget} "$out/${workspacePrebuildPackageDir}/src/lib.rs"
      '';
      workspacePrebuildCargoArgs =
        lib.concatStringsSep " " (
          [
            "--offline"
            "-p ${lib.escapeShellArg workspacePrebuildPackage}"
            "-p ${lib.escapeShellArg workspaceHackPackage}"
          ]
          ++ lib.optional (craneWorkspacePrebuildFeatureArgs != "") craneWorkspacePrebuildFeatureArgs
        );
      workspaceConsumerPackageFor = package: "gammaloop-ci-consumer-${package}";
      workspaceConsumerPackageDirFor = package: "crates/${workspaceConsumerPackageFor package}";
      workspaceConsumerCargoTomlFor = package: dependencyPackages: pkgs.writeText "${workspaceConsumerPackageFor package}-Cargo.toml" ''
        [package]
        name = "${workspaceConsumerPackageFor package}"
        version = "0.1.0"
        edition = "2024"
        publish = false

        [lib]
        path = "src/lib.rs"

        [dependencies]
        ${lib.concatMapStringsSep "\n" (dependencyPackage: let
          features = lib.filter (feature: !lib.hasInfix "/" feature) (craneCiFeaturesFor dependencyPackage);
          featureEntry = lib.optionalString (features != []) ", features = ${builtins.toJSON features}";
        in ''
          ${dependencyPackage} = { path = "${workspacePrebuildDependencyPathFor dependencyPackage}"${featureEntry} }
        '') dependencyPackages}
      '';
      workspaceConsumerSourceScriptFor = package: dependencyPackages: ''
        install -D -m 0644 ${workspaceConsumerCargoTomlFor package dependencyPackages} "$out/${workspaceConsumerPackageDirFor package}/Cargo.toml"
        install -D -m 0644 ${dummyCargoTarget} "$out/${workspaceConsumerPackageDirFor package}/src/lib.rs"
      '';
      cargoPackageDependencyModeArgsFor = package: let
        sourcePackages = workspaceNormalSourcePackageNamesFor package;
        anchorPackages =
          lib.optionals (
            package != workspaceHackPackage
            && builtins.elem package workspaceGraph.symbolica_normal_packages
          ) [workspaceHackPackage];
        crossFeatures =
          sortedUnique (lib.concatMap (
              featurePackage:
                lib.filter (feature: lib.hasInfix "/" feature) (craneCiFeaturesFor featurePackage)
            )
            (sortedUnique (sourcePackages ++ anchorPackages)));
        featureArgs = cargoFeatureArgs crossFeatures;
      in
        lib.concatStringsSep " " (
          [
            "--offline"
            "-p ${lib.escapeShellArg (workspaceConsumerPackageFor package)}"
          ]
          ++ map (anchorPackage: "-p ${lib.escapeShellArg anchorPackage}") anchorPackages
          ++ lib.optional (featureArgs != "") featureArgs
        );
      cargoPackageCiArgsFor = package:
        cargoPackageArgsWithFeaturePackagesFor {
          inherit package;
          featurePackages = workspaceNormalSourcePackageNamesFor package;
          featuresFor = craneCiFeaturesFor;
        };
      testDependencyFeatureAnchorPackageFor = context: "gammaloop-ci-test-dependencies-${context.key}";
      testDependencyFeatureAnchorPackageDirFor = context: "crates/${testDependencyFeatureAnchorPackageFor context}";
      testDependencyFeatureAnchorDependencyPackagesFor = context:
        lib.filter workspacePackageHasLibTarget context.sourcePackages;
      testDependencyFeatureAnchorCargoTomlFor = context: pkgs.writeText "${testDependencyFeatureAnchorPackageFor context}-Cargo.toml" ''
        [package]
        name = "${testDependencyFeatureAnchorPackageFor context}"
        version = "0.1.0"
        edition = "2024"
        publish = false

        [lib]
        path = "src/lib.rs"

        [dependencies]
        ${lib.concatMapStringsSep "\n" (package: let
          features = lib.filter (feature: !lib.hasInfix "/" feature) context.features.${package};
          featureEntry = lib.optionalString (features != []) ", features = ${builtins.toJSON features}";
        in ''
          ${package} = { path = "${workspacePrebuildDependencyPathFor package}"${featureEntry} }
        '') (testDependencyFeatureAnchorDependencyPackagesFor context)}
      '';
      testDependencyFeatureAnchorSourceScriptFor = context: prefix: ''
        install -D -m 0644 ${testDependencyFeatureAnchorCargoTomlFor context} "${prefix}${testDependencyFeatureAnchorPackageDirFor context}/Cargo.toml"
        install -D -m 0644 ${dummyCargoTarget} "${prefix}${testDependencyFeatureAnchorPackageDirFor context}/src/lib.rs"
        touch -d @0 "${prefix}${testDependencyFeatureAnchorPackageDirFor context}/Cargo.toml" "${prefix}${testDependencyFeatureAnchorPackageDirFor context}/src/lib.rs"
      '';
      testBinaryFeatureAnchorPackageFor = context: "gammaloop-ci-test-binary-dependencies-${context.key}";
      testBinaryFeatureAnchorPackageDirFor = context: "crates/${testBinaryFeatureAnchorPackageFor context}";
      testBinaryFeatureAnchorDevDependenciesFor = context:
        lib.foldl'
        lib.recursiveUpdate
        {}
        (map (package:
            builtins.removeAttrs (
              (workspaceManifestFor workspaceMemberPackageDirs.${package})."dev-dependencies" or {}
            )
            context.componentPackages)
          context.componentPackages);
      testBinaryFeatureAnchorCrossFeaturesFor = context:
        lib.foldl' (
          features: feature: let
            parts = lib.splitString "/" feature;
            dependency = builtins.head parts;
          in
            features
            // {
              ${dependency} = sortedUnique ((features.${dependency} or []) ++ [(builtins.elemAt parts 1)]);
            }
        ) {} (lib.concatMap (
            package: lib.filter (feature: lib.hasInfix "/" feature) context.features.${package}
          )
          context.featurePackages);
      testBinaryFeatureAnchorDependenciesFor = context: let
        inheritedDependencies = lib.listToAttrs (map (package: {
            name = package;
            value.path = workspacePrebuildDependencyPathFor package;
          })
          (lib.subtractLists context.componentPackages (testDependencyFeatureAnchorDependencyPackagesFor context)));
        rawDependencies =
          inheritedDependencies
          // (testBinaryFeatureAnchorDevDependenciesFor context)
          // lib.optionalAttrs (context.anchorPackages != []) {
            ${workspaceHackPackage} = {path = "../${workspaceHackPackage}";};
          };
        crossFeatures = testBinaryFeatureAnchorCrossFeaturesFor context;
        dependencyNames = sortedUnique (builtins.attrNames rawDependencies ++ builtins.attrNames crossFeatures);
      in
        lib.listToAttrs (map (dependency: let
            rawDependency = rawDependencies.${dependency} or {workspace = true;};
            dependencyAttrs =
              if builtins.isAttrs rawDependency
              then rawDependency
              else {version = rawDependency;};
            features = sortedUnique (
              (dependencyAttrs.features or [])
              ++ lib.filter (feature: !lib.hasInfix "/" feature) (context.features.${dependency} or [])
              ++ (crossFeatures.${dependency} or [])
            );
          in {
            name = dependency;
            value =
              dependencyAttrs
              // lib.optionalAttrs (features != []) {inherit features;};
          })
          dependencyNames);
      testBinaryFeatureAnchorCargoTomlFor = context:
        (pkgs.formats.toml {}).generate "${testBinaryFeatureAnchorPackageFor context}-Cargo.toml" {
          package = {
            name = testBinaryFeatureAnchorPackageFor context;
            version = "0.1.0";
            edition = "2024";
            publish = false;
          };
          lib.path = "src/lib.rs";
          dependencies = testBinaryFeatureAnchorDependenciesFor context;
        };
      testBinaryFeatureAnchorSourceScriptFor = context: prefix: ''
        install -D -m 0644 ${testBinaryFeatureAnchorCargoTomlFor context} "${prefix}${testBinaryFeatureAnchorPackageDirFor context}/Cargo.toml"
        install -D -m 0644 ${dummyCargoTarget} "${prefix}${testBinaryFeatureAnchorPackageDirFor context}/src/lib.rs"
        touch -d @0 "${prefix}${testBinaryFeatureAnchorPackageDirFor context}/Cargo.toml" "${prefix}${testBinaryFeatureAnchorPackageDirFor context}/src/lib.rs"
      '';
      testBinaryFeatureAnchorDependencyPathFor = context: package:
        if lib.hasPrefix "crates/" workspaceMemberPackageDirs.${package}
        then "../${testBinaryFeatureAnchorPackageFor context}"
        else "../${testBinaryFeatureAnchorPackageDirFor context}";
      testBinaryFeatureAnchorDevDependencyFor = context: package: pkgs.writeText "${testBinaryFeatureAnchorPackageFor context}-${package}-dev-dependency.toml" ''

        [dev-dependencies.${testBinaryFeatureAnchorPackageFor context}]
        path = "${testBinaryFeatureAnchorDependencyPathFor context package}"
      '';
      testBinaryFeatureAnchorDevDependencyScriptFor = context: package: prefix: ''
        cat ${testBinaryFeatureAnchorDevDependencyFor context package} >> "${prefix}${workspaceMemberPackageDirs.${package}}/Cargo.toml"
        touch -d @0 "${prefix}${workspaceMemberPackageDirs.${package}}/Cargo.toml"
      '';
      cargoTestDependencyArgsFor = context: let
        selectedPackages = sortedUnique (context.anchorPackages ++ [
          (testDependencyFeatureAnchorPackageFor context)
        ]);
        featureArgs = cargoQualifiedFeatureArgsFor context.featurePackages (package: context.features.${package});
      in
        lib.concatStringsSep " " (
          [
            "--offline"
          ]
          ++ map (selectedPackage: "-p ${lib.escapeShellArg selectedPackage}") selectedPackages
          ++ lib.optional (featureArgs != "") featureArgs
        );
      cargoTestContextArgsFor = context: packages: let
        featureArgs = cargoQualifiedFeatureArgsFor context.featurePackages (featurePackage: context.features.${featurePackage});
      in
        lib.concatStringsSep " " (
          ["--offline"]
          ++ map (package: "-p ${lib.escapeShellArg package}") (sortedUnique packages)
          ++ lib.optional (featureArgs != "") featureArgs
        );
      cargoTestBinaryDependencyArgsFor = context:
        cargoTestContextArgsFor context ([(testBinaryFeatureAnchorPackageFor context)] ++ context.componentPackages);
      cargoTestBinaryArgsFor = context: package: let
        featureArgs = cargoQualifiedFeatureArgsFor [package] (featurePackage: context.features.${featurePackage});
      in
        lib.concatStringsSep " " (
          [
            "--offline"
            "-p ${lib.escapeShellArg package}"
          ]
          ++ lib.optional (featureArgs != "") featureArgs
        );
      cranePythonExtraFeatureSets = {
        "gammaloop-api" = ["python_abi" "pyo3-extension-module"];
      };
      cranePythonFeaturesFor = package:
        sortedUnique (craneCiFeaturesFor package ++ (cranePythonExtraFeatureSets.${package} or []));
      # The Python feature set enables optional workspace dependencies that are
      # absent from the default resolved closure used by ordinary package builds.
      cranePythonSourcePackageNames =
        workspaceDependencyClosureFor workspaceDependencyNamesFor "gammaloop-api";
      cranePythonDependencySourcePackageNames =
        lib.filter (sourcePackage: sourcePackage != "gammaloop-api") cranePythonSourcePackageNames;
      cranePythonExtraDependencySourcePackageNames =
        lib.subtractLists
        (lib.filter
          (sourcePackage: sourcePackage != "gammaloop-api")
          (workspaceNormalSourcePackageNamesFor "gammaloop-api"))
        cranePythonDependencySourcePackageNames;
      cranePythonSrc = workspacePackageSrcForSourcePackages {
        sourcePackages = cranePythonSourcePackageNames;
        packageSourcePackages = ["gammaloop-api"];
      };
      cranePythonCargoArgs = let
        featurePackages = cranePythonSourcePackageNames;
        selectedFeaturePackages =
          lib.filter (
            featurePackage:
              (cranePythonFeaturesFor featurePackage) != []
              && !workspacePackageIsProcMacro featurePackage
          )
          featurePackages;
        selectedPackages = sortedUnique (["gammaloop-api" workspaceHackPackage] ++ selectedFeaturePackages);
        featureArgs = cargoQualifiedFeatureArgsFor featurePackages cranePythonFeaturesFor;
      in
        lib.concatStringsSep " " (
          [
            "--locked"
            "--no-default-features"
          ]
          ++ map (package: "-p ${lib.escapeShellArg package}") selectedPackages
          ++ lib.optional (featureArgs != "") featureArgs
        );

      guppyFeatureMapFor = featuresFor:
        builtins.toJSON (lib.listToAttrs (map (package: {
            name = package;
            value = featuresFor package;
          })
          workspaceMemberPackages));

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
          cargoExtraArgs = "--locked --workspace ${craneWorkspacePrebuildFeatureArgs}";
          # NixCI provides the runtime Symbolica license, not the compile-time
          # OEM key consumed by gammalooprs' activate_oem_license! path.
          NO_SYMBOLICA_OEM_LICENSE = "1";

          PYO3_PYTHON = "${pkgs.python313}/bin/python3";
          PYTHONPATH = "${pkgs.python313}/lib/python3.13/site-packages";
        };

      licensePreCheck = ''
        if [ -z "''${SYMBOLICA_LICENSE:-}" ]; then
          echo "Missing SYMBOLICA_LICENSE environment variable" >&2
          exit 1
        fi
      '';

      gammaloop-cli = pkgs.runCommand "gammaloop-api-${apiMeta.version}" {
        nativeBuildInputs = [
          pkgs.coreutils
          pkgs.findutils
        ];
      } ''
        target="${gammaloopApiPackageArtifacts}/target/${ciCargoProfile}"
        binary="$target/gammaloop"

        if [ ! -x "$binary" ]; then
          echo "Could not find Crane-built gammaloop binary in $target" >&2
          exit 1
        fi

        library="$(
          find "$target" -maxdepth 2 -type f \
            \( -name 'libgammaloop_api*.so' -o -name 'libgammaloop_api*.dylib' \) \
            | sort \
            | head -n 1
        )"
        if [ -z "$library" ]; then
          echo "Could not find Crane-built gammaloop-api shared library in $target" >&2
          exit 1
        fi

        install -D -m 0755 "$binary" "$out/bin/gammaloop"
        install -D -m 0644 "$library" "$out/lib/$(basename "$library")"
      '';
      gammaloop-python-lib = craneLib.buildPackage (ciArgs
        // {
          cargoArtifacts = cranePythonBuildArtifacts;
          pname = "gammaloop-api-python";
          src = cranePythonSrc;
          cargoExtraArgs = cranePythonCargoArgs;
          doCheck = false;
          postPatch = workspaceMissingCargoTargetsScript;
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
      clinnetArgs = ciArgs
        // {
          pname = "clinnet";
          inherit (clinnetMeta) version;
          src = workspacePackageSrcFor "clinnet";
          cargoExtraArgs = cargoPackageCiArgsFor "clinnet";
          doCheck = false;
          postPatch = workspaceMissingCargoTargetsScript;
        };
      drawingTypstBundleAssets = ''
        mkdir -p crates/linnest/typst/src crates/kurvst/typst/src
        cp -R ${linnest-wasm}/templates/crates/linnest/typst/src/. crates/linnest/typst/src/
        cp ${linnest-wasm}/templates/crates/linnest/typst/typst.toml crates/linnest/typst/typst.toml
        cp ${linnest-wasm}/templates/crates/linnest/typst/linnest.wasm crates/linnest/typst/linnest.wasm
        cp -R ${linnest-wasm}/templates/crates/kurvst/typst/src/. crates/kurvst/typst/src/
        cp ${linnest-wasm}/templates/crates/kurvst/typst/typst.toml crates/kurvst/typst/typst.toml
        cp ${linnest-wasm}/templates/crates/kurvst/typst/kurvst.wasm crates/kurvst/typst/kurvst.wasm
      '';

      clinnetCargoArtifacts = craneLib.buildDepsOnly (clinnetArgs
        // {
          preBuild = drawingTypstBundleAssets;
        });

      clinnet-cli = craneLib.buildPackage (clinnetArgs
        // {
          cargoArtifacts = clinnetCargoArtifacts;
          doNotLinkInheritedArtifacts = true;
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

      devShellPackages = with pkgs;
        [
          tdf
          cargo-flamegraph
          yaml-language-server
          just
          gitMinimal
          jujutsu
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
          docsTypst
          linnetPython
          docsPkgs.roboto
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
        ++ lib.optionals (!pkgs.stdenv.hostPlatform.isDarwin) [
          valgrind
        ];

      mkDevShell = craneLibForShell: extraPackages: rustFlags:
        craneLibForShell.devShell ({
            # checks = self.checks.${system};

            RUST_SRC_PATH = "${pkgs.rustPlatform.rustLibSrc}";
            GLIBC_TUNABLES = "glibc.rtld.optional_static_tls=10000";
            TYPST_FONT_PATHS = docsFontPath;
            # `typst.withPackages` injects this only into its executable wrapper.
            # The persistent Rust renderer and typst-py need the same package tree directly.
            TYPST_PACKAGE_CACHE_PATH = "${docsTypst}/lib/typst/packages";

            CC = nixCc;
            CXX = nixCxx;
            "${cargoLinkerVar}" = nixCc;

            LD_LIBRARY_PATH = runtimeLibPath;
            DYLD_LIBRARY_PATH = runtimeLibPath;

            # shellHook = ''
            #   export CC="${nixCc}"
            #   export CXX="${nixCxx}"
            #   export ${cargoLinkerVar}="${nixCc}"
            # '';

            packages = devShellPackages ++ extraPackages;
          }
          // lib.optionalAttrs (rustFlags != null) {
            RUSTFLAGS = rustFlags;
          });
      in
        craneLib.mkCargoDerivation (
          cleanedArgs
          // {
            inherit doCheck;

            src = dummySrc;
            pnameSuffix = "-deps";
            pname = args.pname or crateName.pname;
            version = args.version or crateName.version;
            nativeBuildInputs =
              (cleanedArgs.nativeBuildInputs or [])
              ++ lib.optionals (args.stripWorkspaceArtifacts or false) [
                pkgs.gnutar
                pkgs.rsync
                pkgs.zstd
              ];

            cargoArtifacts = args.cargoArtifacts or null;
            doNotLinkInheritedArtifacts = args.doNotLinkInheritedArtifacts or true;
            cargoVendorDir = args.cargoVendorDir or (craneLib.vendorCargoDeps argsMaybeDummySrcOverride);

            env = (args.env or {}) // {
              CRANE_BUILD_DEPS_ONLY = ((args.env or {}).CRANE_BUILD_DEPS_ONLY or 1);
            };

            postPatch =
              (args.postPatch or "")
              + ''

                ${normalizeWorkspaceHackBuildScriptTimestampScript}
              '';

            preBuild =
              (args.preBuild or "")
              + lib.optionalString (preBuildWorkspaceArtifactStripPackages != []) (stripSelectedWorkspaceCargoArtifactsScript preBuildWorkspaceArtifactStripPackages);

            buildPhaseCargoCommand =
              args.buildPhaseCargoCommand or ''
                ${cargoCheckCommand} ${cargoExtraArgs} ${cargoCheckExtraArgs}
                ${cargoBuildCommand} ${cargoExtraArgs}
              '';

            checkPhaseCargoCommand =
              args.checkPhaseCargoCommand or ''
                ${cargoTestCommand} ${cargoExtraArgs} ${cargoTestExtraArgs}
              '';

            postBuild =
              (args.postBuild or "")
              + stripWorkspaceArtifactsScriptText;

            preFixup =
              (args.preFixup or "")
              + lib.optionalString (args.stripWorkspaceArtifacts or false) ''
                if [ -f "$out/target.tar.zst" ] && { [ -e "$out/target.tar.zst.prev" ] || [ -L "$out/target.tar.zst.prev" ]; }; then
                  tmp="$(mktemp -d)"
                  mkdir -p "$tmp/target"

                  make_target_writable() {
                    chmod -R u+w "$tmp/target"
                  }

                  unpack_artifact() {
                    local artifact="$1"

                    if [ -d "$artifact" ] && [ -f "$artifact/target.tar.zst" ]; then
                      artifact="$artifact/target.tar.zst"
                    elif [ -d "$artifact" ] && [ -d "$artifact/target" ]; then
                      artifact="$artifact/target"
                    fi

                    if [ -f "$artifact" ]; then
                      if [ -e "$artifact.prev" ] || [ -L "$artifact.prev" ]; then
                        unpack_artifact "$(realpath "$artifact.prev")"
                      fi
                      make_target_writable
                      zstd -d "$artifact" --stdout | tar --no-same-permissions -x -C "$tmp/target"
                      make_target_writable
                    elif [ -d "$artifact" ]; then
                      make_target_writable
                      rsync -a --chmod=u+w "$artifact/" "$tmp/target/"
                      make_target_writable
                    else
                      echo "unsupported cargo artifact path: $artifact" >&2
                      exit 1
                    fi
                  }

                  unpack_artifact "$(realpath "$out/target.tar.zst.prev")"
                  unpack_artifact "$out/target.tar.zst"

                  (
                    cd "$tmp"
                    ${stripWorkspaceArtifactsScriptText}
                  )

                  # Match Crane's artifact and Nix source timestamps so Cargo
                  # does not treat the compacted target as stale.
                  tar --sort=name --mtime=@1 --owner=0 --group=0 --numeric-owner -C "$tmp/target" -cf - . \
                    | zstd -T0 --stdout > "$out/target.tar.zst.tmp"
                  mv "$out/target.tar.zst.tmp" "$out/target.tar.zst"
                  rm -f "$out/target.tar.zst.prev"
                  rm -rf "$tmp"
                fi
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
        workspaceHackBuildArtifacts;

      workspaceHackDependencyArtifacts = buildDepsOnlyWithArtifacts (ciArgs
        // {
          pname = "gammaloop-crate-${workspaceHackPackage}";
          src = workspacePackageSrcFor workspaceHackPackage;
          cargoExtraArgs = cargoPackageCiArgsFor workspaceHackPackage;
          buildPhaseCargoCommand = "cargoWithProfile build ${cargoPackageCiArgsFor workspaceHackPackage}";
          checkPhaseCargoCommand = "";
          doCheck = false;
          preBuildWorkspaceArtifactStripPackages = [workspaceHackPackage];
          stripWorkspaceArtifacts = true;
          extraDummyScript = workspaceAllDummyCargoTargetsScript;
        });

      workspaceHackBuildArtifacts = craneLib.cargoBuild (ciArgs
        // {
          cargoArtifacts = workspaceHackDependencyArtifacts;
          doNotLinkInheritedArtifacts = true;
          pname = "gammaloop-crate-${workspaceHackPackage}";
          src = workspacePackageSrcFor workspaceHackPackage;
          cargoExtraArgs = cargoPackageCiArgsFor workspaceHackPackage;
          postPatch = workspaceMissingCargoTargetsScript;
        });

      cranePackageDependencyModeDependencyArtifacts = lib.fix (self:
        lib.genAttrs workspacePackages (package:
          if package == workspaceHackPackage
          then workspaceHackDependencyArtifacts
          else if !workspacePackageHasLibTarget package
          then null
          else let
            dependencyArtifactFor = dependency:
              if dependency == workspaceHackPackage
              then workspaceHackBuildArtifacts
              else cranePackageDependencyModeArtifacts.${dependency};
            dependencySourcePackages =
              lib.filter (sourcePackage: sourcePackage != package) (workspaceNormalSourcePackageNamesFor package);
            preservedWorkspaceArtifactPackages =
              lib.filter workspacePackageIsProcMacro dependencySourcePackages;
          in
            buildDepsOnlyWithArtifacts ((builtins.removeAttrs ciArgs ["src"])
              // {
                cargoArtifacts = mergeCargoArtifactsOrNull "gammaloop-crate-${package}-dependency-deps-inputs" (
                  [cargoArtifacts]
                  ++ rootPackageDependencyArtifactsFor package
                  ++ map dependencyArtifactFor dependencySourcePackages
                );
                pname = "gammaloop-crate-${package}-dependency-deps";
                src = workspacePackageSrcForSourcePackages {
                  sourcePackages = dependencySourcePackages;
                };
                buildPhaseCargoCommand = "cargoWithProfile build ${cargoPackageDependencyModeArgsFor package}";
                checkPhaseCargoCommand = "";
                doCheck = false;
                stripWorkspaceArtifacts = true;
                inherit preservedWorkspaceArtifactPackages;
                extraDummyScript = ''
                  ${workspaceDependencyDummyCargoTargetsScriptFor package}
                  ${workspaceConsumerSourceScriptFor package dependencySourcePackages}
                '';
                postPatch = workspaceMissingCargoTargetsScript;
              })));

      cranePackageDependencyModeArtifacts = lib.fix (self:
        lib.genAttrs workspacePackages (package:
          if package == workspaceHackPackage
          then workspaceHackBuildArtifacts
          else if !workspacePackageHasLibTarget package
          then null
          else let
            dependencyArtifactFor = dependency:
              if dependency == workspaceHackPackage
              then workspaceHackBuildArtifacts
              else self.${dependency};
            sourcePackages = workspaceNormalSourcePackageNamesFor package;
            preservedWorkspaceArtifactPackages = sourcePackages;
          in
            buildDepsOnlyWithArtifacts ((builtins.removeAttrs ciArgs ["src"])
              // {
                cargoArtifacts = mergeCargoArtifactsOrNull "gammaloop-crate-${package}-dependency-inputs" (
                  [cranePackageDependencyModeDependencyArtifacts.${package}]
                  ++ map dependencyArtifactFor (lib.filter (sourcePackage: sourcePackage != package) sourcePackages)
                );
                pname = "gammaloop-crate-${package}-dependency";
                src = workspacePackageSrcForSourcePackages {
                  inherit sourcePackages;
                };
                buildPhaseCargoCommand = "cargoWithProfile build ${cargoPackageDependencyModeArgsFor package}";
                checkPhaseCargoCommand = "";
                doCheck = false;
                preBuildWorkspaceArtifactStripPackages = [package];
                stripWorkspaceArtifacts = true;
                inherit preservedWorkspaceArtifactPackages;
                extraDummyScript = ''
                  ${workspacePackageSourceRestoreInDummySrcScriptFor sourcePackages}
                  ${workspaceConsumerSourceScriptFor package sourcePackages}
                '';
                postPatch = workspaceMissingCargoTargetsScript;
              })));

      # Public crate-deps outputs are the reusable workspace crate artifacts.
      # Build them through the consumer-anchor path so each package's real
      # artifact is preserved for downstream package and test derivations.
      cranePackageDependencyArtifacts = cranePackageDependencyModeArtifacts;

      cranePackageBuildArtifacts = lib.fix (self:
        lib.genAttrs workspacePackages (package:
          if package == workspaceHackPackage
          then workspaceHackBuildArtifacts
          else
            craneLib.cargoBuild (ciArgs
              // {
                cargoArtifacts =
                  if cranePackageDependencyArtifacts.${package} == null
                  then cargoArtifacts
                  else cranePackageDependencyArtifacts.${package};
                pname = "gammaloop-crate-${package}";
                src = workspacePackageSrcFor package;
                cargoExtraArgs =
                  cargoPackageCiArgsFor package
                  + lib.optionalString (package == "gammaloop-api") " --lib --bins";
                preBuild = lib.optionalString (package == "gammaloop-api") drawingTypstBundleAssets;
                postPatch = workspaceMissingCargoTargetsScript;
              })));

      gammaloopApiPackageArtifacts = mergeCargoArtifacts "gammaloop-api-package-artifacts" [
        cranePackageBuildArtifacts."gammaloop-api"
      ];

      workspaceBaseTestContexts =
        map (package: workspaceTestContextFor {packages = [package];}) workspacePackages;
      nextestPackageTestContexts = lib.concatMap (target:
        map (package:
          workspaceTestContextFor {
            packages = [package];
            extraFeatures = target.extraFeatures or {};
          })
        target.packages)
      checkedNextestPackageGroups;
      workspaceTestContexts = lib.listToAttrs (map (context: {
          name = context.key;
          value = context;
        })
        (workspaceBaseTestContexts ++ nextestPackageTestContexts));

      craneTestDependencyArtifacts = lib.fix (self:
        lib.mapAttrs (_: context:
          if context.componentPackages == [workspaceHackPackage]
          then workspaceHackBuildArtifacts
          else let
            sourcePackages = lib.filter workspacePackageHasLibTarget context.sourcePackages;
            dependencyContexts = workspaceTestDependencyContextsFor context;
            procMacroPackages = lib.filter workspacePackageIsProcMacro context.componentPackages;
            hostAnchorPackage = "${testBinaryFeatureAnchorPackageFor context}-host";
            hostAnchorPackageDir = "crates/${hostAnchorPackage}";
            anchorConsumerDependencies = {
              ${testBinaryFeatureAnchorPackageFor context}.path = "../${testBinaryFeatureAnchorPackageFor context}";
            } // lib.listToAttrs (map (package: let
                features = lib.filter (feature: feature != "default" && !lib.hasInfix "/" feature) context.features.${package};
              in {
                name = package;
                value = {
                  path = workspacePrebuildDependencyPathFor package;
                } // lib.optionalAttrs (features != []) {inherit features;};
              })
              procMacroPackages);
            hostAnchorCargoToml = (pkgs.formats.toml {}).generate "${hostAnchorPackage}-Cargo.toml" {
              package = {
                name = hostAnchorPackage;
                version = "0.1.0";
                edition = "2024";
                publish = false;
              };
              lib = {
                path = "src/lib.rs";
                "proc-macro" = true;
              };
              dependencies = anchorConsumerDependencies;
            };
            targetFacadeFor = package: let
              manifest = workspaceManifestFor workspaceMemberPackageDirs.${package};
              packageName = "${testBinaryFeatureAnchorPackageFor context}-target-${package}";
              packageDir = "crates/${packageName}";
              localFeatures = lib.filter (feature: feature != "default" && !lib.hasInfix "/" feature) context.features.${package};
              resolvedDependencies = testBinaryFeatureAnchorDependenciesFor context;
              rawDependencies = (manifest.dependencies or {}) // (manifest."dev-dependencies" or {});
              dependencies = lib.mapAttrs (dependency: rawDependency: let
                  dependencyAttrs =
                    if builtins.isAttrs rawDependency
                    then rawDependency
                    else {version = rawDependency;};
                  features = sortedUnique (
                    (dependencyAttrs.features or [])
                    ++ (resolvedDependencies.${dependency}.features or [])
                  );
                in
                  dependencyAttrs // lib.optionalAttrs (features != []) {inherit features;})
                rawDependencies
                // anchorConsumerDependencies;
              cargoToml = (pkgs.formats.toml {}).generate "${packageName}-Cargo.toml" {
                package = {
                  name = packageName;
                  version = "0.1.0";
                  edition = "2024";
                  publish = false;
                };
                lib.path = "src/lib.rs";
                inherit dependencies;
                features = manifest.features or {};
              };
            in {
              inherit cargoToml localFeatures packageDir packageName;
            };
            targetFacades = map targetFacadeFor procMacroPackages;
            targetFacadeFeatures = lib.concatMap (facade:
              map (feature: "${facade.packageName}/${feature}") facade.localFeatures)
            targetFacades;
            anchorConsumerArgs = lib.concatStringsSep " " (
              ["--offline" "-p ${lib.escapeShellArg hostAnchorPackage}"]
              ++ map (facade: "-p ${lib.escapeShellArg facade.packageName}") targetFacades
              ++ lib.optional (targetFacadeFeatures != []) (cargoFeatureArgs targetFacadeFeatures)
            );
            hasProcMacro = procMacroPackages != [];
            buildPhaseCargoCommand = ''
              cargoWithProfile build ${cargoTestDependencyArgsFor context} --lib
              cargoWithProfile build ${cargoTestBinaryDependencyArgsFor context} --lib
            '';
            postPatch = ''
              ${workspaceMissingCargoTargetsScript}
              ${testDependencyFeatureAnchorSourceScriptFor context ""}
              ${testBinaryFeatureAnchorSourceScriptFor context ""}
            '';
          in
            buildDepsOnlyWithArtifacts ((builtins.removeAttrs ciArgs ["src"])
              // {
                cargoArtifacts = mergeCargoArtifacts "gammaloop-crate-test-dependencies-${context.key}-inputs" (
                  [
                    cargoArtifacts
                    workspaceHackBuildArtifacts
                  ]
                  ++ map (dependencyContext: self.${dependencyContext.key}) dependencyContexts
                );
                pname = "gammaloop-crate-test-dependencies-${context.key}";
                dummySrc = workspacePackageSrcForSourcePackages {
                  inherit sourcePackages;
                };
                # Rebuild workspace library crates against the final merge's dependency metadata.
                inherit buildPhaseCargoCommand postPatch;
                checkPhaseCargoCommand = "";
                doCheck = false;
              }
              // lib.optionalAttrs hasProcMacro {
                buildPhaseCargoCommand = buildPhaseCargoCommand + ''
                  cargoWithProfile build ${anchorConsumerArgs} --lib
                '';
                postPatch = postPatch + ''
                  install -D -m 0644 ${hostAnchorCargoToml} ${hostAnchorPackageDir}/Cargo.toml
                  install -D -m 0644 ${dummyProcMacroCargoTarget} ${hostAnchorPackageDir}/src/lib.rs
                  touch -d @0 ${hostAnchorPackageDir}/Cargo.toml ${hostAnchorPackageDir}/src/lib.rs
                  ${lib.concatMapStringsSep "\n" (facade: ''
                    install -D -m 0644 ${facade.cargoToml} ${facade.packageDir}/Cargo.toml
                    install -D -m 0644 ${dummyCargoTarget} ${facade.packageDir}/src/lib.rs
                    touch -d @0 ${facade.packageDir}/Cargo.toml ${facade.packageDir}/src/lib.rs
                  '')
                  targetFacades}
                '';
              }
              // lib.optionalAttrs context.usesPythonModule {
                nativeBuildInputs = (ciArgs.nativeBuildInputs or []) ++ [nextestPython];
                PYO3_PYTHON = "${nextestPython}/bin/python3";
                PYTHON = "${nextestPython}/bin/python3";
                PYTHONPATH = "${gammaloop-python-module}/${pythonSitePackages}:${nextestPython}/${pythonSitePackages}";
              }))
        workspaceTestContexts);

      craneTestBinaryArtifactFor = context: package:
        buildDepsOnlyWithArtifacts ((builtins.removeAttrs ciArgs ["src"])
          // {
            cargoArtifacts = mergeCargoArtifacts "gammaloop-crate-test-binaries-${package}-${context.key}-inputs" [
              cargoArtifacts
              craneTestDependencyArtifacts.${context.key}
            ];
            pname = "gammaloop-crate-test-binaries-${package}-${context.key}";
            dummySrc = workspacePackageSrcForSourcePackages {
              sourcePackages = context.sourcePackages;
              packageSourcePackages = [package];
              testSourcePackages = [package];
            };
            buildPhaseCargoCommand = "cargoWithProfile test --no-run ${cargoTestBinaryArgsFor context package}";
            checkPhaseCargoCommand = "";
            doCheck = false;
            postPatch = ''
              ${workspaceMissingCargoTargetsScript}
              ${testBinaryFeatureAnchorSourceScriptFor context ""}
              ${testBinaryFeatureAnchorDevDependencyScriptFor context package ""}
            '';
          }
          // lib.optionalAttrs context.usesPythonModule {
            nativeBuildInputs = (ciArgs.nativeBuildInputs or []) ++ [nextestPython];
            PYO3_PYTHON = "${nextestPython}/bin/python3";
            PYTHON = "${nextestPython}/bin/python3";
            PYTHONPATH = "${gammaloop-python-module}/${pythonSitePackages}:${nextestPython}/${pythonSitePackages}";
          });

      craneTestBinaryArtifacts = lib.genAttrs workspacePackages (package:
        if package == workspaceHackPackage
        then workspaceHackBuildArtifacts
        else craneTestBinaryArtifactFor (workspaceTestContextFor {packages = [package];}) package);

      workspaceCargoCheck = craneLib.mkCargoDerivation (ciArgs
        // {
          cargoArtifacts = cargoArtifacts;
          pname = "gammaloop-workspace-check";
          src = workspaceTestSrc;
          doNotLinkInheritedArtifacts = true;
          doInstallCargoArtifacts = false;
          buildPhaseCargoCommand = ''
            mkdir -p "$out"
            if [ -d target ]; then
              chmod -R u+w target
            fi
            cargoWithProfile check ${ciArgs.cargoExtraArgs} --all-targets
          '';
          checkPhaseCargoCommand = "";
          doCheck = false;
          installPhaseCommand = "";
        });

      workspaceClippyCheck = craneLib.mkCargoDerivation (ciArgs
        // {
          cargoArtifacts = cargoArtifacts;
          pname = "gammaloop-workspace-clippy";
          src = workspaceTestSrc;
          doNotLinkInheritedArtifacts = true;
          doInstallCargoArtifacts = false;
          buildPhaseCargoCommand = ''
            mkdir -p "$out"
            if [ -d target ]; then
              chmod -R u+w target
            fi
            cargoWithProfile clippy ${ciArgs.cargoExtraArgs} --all-targets --no-deps -- --deny warnings
          '';
          checkPhaseCargoCommand = "";
          doCheck = false;
          installPhaseCommand = "";
        });

      workspaceDocCheck = craneLib.mkCargoDerivation (ciArgs
        // {
          cargoArtifacts = cargoArtifacts;
          pname = "gammaloop-workspace-doc";
          src = workspaceTestSrc;
          doNotLinkInheritedArtifacts = true;
          doInstallCargoArtifacts = false;
          buildPhaseCargoCommand = ''
            mkdir -p "$out"
            if [ -d target ]; then
              chmod -R u+w target
            fi
            cargoWithProfile doc ${ciArgs.cargoExtraArgs} --no-deps
          '';
          checkPhaseCargoCommand = "";
          doCheck = false;
          installPhaseCommand = "";
        });

      workspaceDoctestCheck = craneLib.mkCargoDerivation (ciArgs
        // {
          cargoArtifacts = cargoArtifacts;
          pname = "gammaloop-workspace-doctest";
          src = workspaceTestSrc;
          doNotLinkInheritedArtifacts = true;
          doInstallCargoArtifacts = false;
          buildPhaseCargoCommand = ''
            mkdir -p "$out"
            if [ -d target ]; then
              chmod -R u+w target
            fi
            cargoWithProfile test --doc ${ciArgs.cargoExtraArgs}
          '';
          checkPhaseCargoCommand = "";
          doCheck = false;
          installPhaseCommand = "";
          preBuild = licensePreCheck;
          SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
        });

      cranePythonDependencyArtifacts = buildDepsOnlyWithArtifacts ((builtins.removeAttrs ciArgs ["src"])
        // {
          cargoArtifacts = mergeCargoArtifacts "gammaloop-python-deps-inputs" (
            [cargoArtifacts]
            ++ map (
              dependency:
                if dependency == workspaceHackPackage
                then workspaceHackBuildArtifacts
                else cranePackageDependencyModeArtifacts.${dependency}
            ) cranePythonDependencySourcePackageNames
          );
          pname = "gammaloop-api-python";
          src = workspacePackageSrcForSourcePackages {
            sourcePackages = cranePythonDependencySourcePackageNames;
          };
          buildPhaseCargoCommand = "cargoWithProfile build ${cranePythonCargoArgs}";
          checkPhaseCargoCommand = "";
          doCheck = false;
          preBuildWorkspaceArtifactStripPackages = ["gammaloop-api"];
          stripWorkspaceArtifacts = true;
          preservedWorkspaceArtifactPackages = cranePythonDependencySourcePackageNames;
          extraDummyScript = ''
            ${workspaceDependencyDummyCargoTargetsScriptFor "gammaloop-api"}
            ${workspacePackageSourceRestoreInDummySrcScriptFor cranePythonExtraDependencySourcePackageNames}
          '';
          postPatch = workspaceMissingCargoTargetsScript;
        });

      cranePythonBuildArtifacts = craneLib.cargoBuild (ciArgs
        // {
          cargoArtifacts = cranePythonDependencyArtifacts;
          pname = "gammaloop-api-python-build";
          src = cranePythonSrc;
          cargoExtraArgs = cranePythonCargoArgs;
          postPatch = workspaceMissingCargoTargetsScript;
        });

      cranePackageOutputs = lib.listToAttrs (map (package: {
          name = "crate-${package}";
          value = cranePackageBuildArtifacts.${package};
        })
        workspacePackages);

      cranePackageDependencyOutputs = lib.listToAttrs (map (package: {
          name = "crate-deps-${package}";
          value =
            nixCiArtifactBarrier
            "crate-deps-${package}"
            cranePackageDependencyArtifacts.${package};
        })
        (lib.filter (package: cranePackageDependencyArtifacts.${package} != null)
          workspacePackages));

      craneTestBinaryPackageOutputs = lib.listToAttrs (map (package: {
          name = "crate-test-binaries-${package}";
          value =
            nixCiArtifactBarrier
            "crate-test-binaries-${package}"
            craneTestBinaryArtifacts.${package};
        })
        workspacePackages);

      craneTestDependencyOutputs = lib.listToAttrs (map (representative: let
          context = workspaceTestContextFor {packages = [representative];};
        in {
          name = "crate-test-dependencies-${representative}";
          value =
            nixCiArtifactBarrier
            "crate-test-dependencies-${representative}"
            craneTestDependencyArtifacts.${context.key};
        })
        workspaceTestDependencyComponentRepresentatives);

      workspaceBuildArtifacts = craneLib.cargoBuild (ciArgs
        // {
          inherit cargoArtifacts;
          pname = "gammaloop-workspace-build-artifacts";
          src = workspaceNonIntegrationTestSrc;
          cargoExtraArgs = "${cargoPackagesArgsFor (lib.subtractLists ["gammaloop-integration-tests"] workspaceMemberPackages) craneTestFeaturesFor} --tests";
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
          name = "feynkit";
          packages = [
            "feynkit"
            "feynkit-cff"
            "feynkit-generator"
            "feynkit-graph"
            "feynkit-kinematics"
            "feynkit-model"
            "feynkit-py"
            "feynkit-tensor"
            "feynkit-ufo"
          ];
        }
        {
          name = "python-api";
          packages = ["gammaloop-integration-tests"];
          runtimeTestSourcePackages = [];
          filter = "package(gammaloop-integration-tests) & binary(test_python_api)";
          extraFeatures."gammaloop-integration-tests" = ["python-api-tests"];
        }
        {
          name = "clinnet";
          packages = ["clinnet"];
        }
        {
          name = "linnet";
          packages = [
            "kurvst"
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
            "symbolica-utils"
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
        "spynso3"
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
      nextestContextFor = target:
        workspaceTestContextFor {
          inherit (target) packages;
          extraFeatures = target.extraFeatures or {};
        };
      nextestPackageContextFor = target: package:
        workspaceTestContextFor {
          packages = [package];
          extraFeatures = target.extraFeatures or {};
        };
      nextestUsesPythonModule = target: (nextestContextFor target).usesPythonModule;
      nextestSourcePackagesFor = target: (nextestContextFor target).sourcePackages;
      nextestSrcFor = target:
        workspacePackageSrcForSourcePackages {
          sourcePackages = nextestSourcePackagesFor target;
          packageSourcePackages = target.packages;
          testSourcePackages = target.packages;
          extraFilesets = [./.config/nextest.toml];
        };
      nextestRuntimeSrcFor = target:
        workspacePackageSrcForSourcePackages {
          sourcePackages = nextestSourcePackagesFor target;
          packageSourcePackages = target.packages;
          testSourcePackages = target.packages;
          runtimeTestSourcePackages = target.runtimeTestSourcePackages or target.packages;
          extraFilesets = [./.config/nextest.toml];
        };

      nextestFeatureArgsFor = target:
        cargoQualifiedFeatureArgsFor target.packages (package:
          (nextestPackageContextFor target package).features.${package});

      nextestCargoArgsFor = target:
        lib.concatStringsSep " " (
          [
            "--offline"
          ]
          ++ map (package: "-p ${lib.escapeShellArg package}") target.packages
          ++ lib.optional (nextestFeatureArgsFor target != "") (nextestFeatureArgsFor target)
        );

      nextestArchiveNameFor = target: package: "gammaloop-nextest-${target.name}-${package}.tar.zst";

      nextestPackageTestBinaryArtifactFor = target: package: let
        context = nextestPackageContextFor target package;
        baseContext = workspaceTestContextFor {packages = [package];};
      in
        if context.key == baseContext.key
        then craneTestBinaryArtifacts.${package}
        else craneTestBinaryArtifactFor context package;

      nextestPackageArchiveFor = target: package: let
        packageTarget = target // {packages = [package];};
        context = nextestPackageContextFor target package;
        packageCargoArtifacts = nextestPackageTestBinaryArtifactFor target package;
      in
        craneLib.mkCargoDerivation (ciArgs
          // {
            pname = "gammaloop-nextest-binaries-${target.name}-${package}";
            src = nextestSrcFor packageTarget;
            cargoArtifacts = packageCargoArtifacts;
            doCheck = false;
            doInstallCargoArtifacts = false;
            nativeBuildInputs =
              (ciArgs.nativeBuildInputs or [])
              ++ [pkgs.cargo-nextest pkgs.form]
              ++ lib.optionals (nextestUsesPythonModule target) [nextestPython];
            postPatch = ''
              ${workspaceMissingCargoTargetsScript}
              ${testBinaryFeatureAnchorSourceScriptFor context ""}
              ${testBinaryFeatureAnchorDevDependencyScriptFor context package ""}

              # Crane only follows file-valued incremental artifact links. The
              # test-binary delta points directly to its materialized input
              # directory, so inherit that base before Crane overlays the
              # package-specific, writable delta in its post-patch hook.
              if [ -d ${packageCargoArtifacts}/target.tar.zst.prev ]; then
                inheritCargoArtifacts "$(realpath ${packageCargoArtifacts}/target.tar.zst.prev)" target
              fi
            '';
            buildPhaseCargoCommand = ''
              mkdir -p "$out"
              if [ -d target ]; then
                chmod -R u+w target
                find target -name .cargo-lock -delete
              fi
              cargo nextest --version
              cargo nextest archive \
                --cargo-profile ${ciCargoProfile} \
                ${nextestCargoArgsFor packageTarget} \
                ${nextestFilterFor packageTarget} \
                --profile ${nextestProfile} \
                --archive-file "$out/${nextestArchiveNameFor target package}"
            '';
            checkPhaseCargoCommand = "";
            installPhaseCommand = "";
          } // lib.optionalAttrs (nextestUsesPythonModule target) {
            PYO3_PYTHON = "${nextestPython}/bin/python3";
            PYTHON = "${nextestPython}/bin/python3";
            PYTHONPATH = "${gammaloop-python-module}/${pythonSitePackages}:${nextestPython}/${pythonSitePackages}";
          });

      nextestArchiveFor = target: let
        packageArchives = lib.genAttrs target.packages (nextestPackageArchiveFor target);
      in
        pkgs.linkFarm "gammaloop-nextest-binaries-${target.name}" (map (package: {
            name = nextestArchiveNameFor target package;
            path = "${packageArchives.${package}}/${nextestArchiveNameFor target package}";
          })
          target.packages);

      nextestBinarySets = lib.listToAttrs (map (target: {
          name = "gammaloop-nextest-binaries-${target.name}";
          value = nextestArchiveFor target;
        })
        checkedNextestPackageGroups);

      nextestContextualTestOutputs = lib.listToAttrs (lib.concatMap (target:
          lib.concatMap (package: let
              context = nextestPackageContextFor target package;
            in [
              {
                name = "crate-test-dependencies-${target.name}-${package}";
                value =
                  nixCiArtifactBarrier
                  "crate-test-dependencies-${target.name}-${package}"
                  craneTestDependencyArtifacts.${context.key};
              }
              {
                name = "crate-test-binaries-${target.name}-${package}";
                value =
                  nixCiArtifactBarrier
                  "crate-test-binaries-${target.name}-${package}"
                  (nextestPackageTestBinaryArtifactFor target package);
              }
            ])
          target.packages)
        (lib.filter (target: target ? extraFeatures) checkedNextestPackageGroups));

      nextestBinarySetForTarget = target: nextestBinarySets."gammaloop-nextest-binaries-${target.name}";

      nextestBinarySetAggregate = pkgs.linkFarm "gammaloop-nextest-binaries" (map (target: {
          name = target.name;
          path = nextestBinarySetForTarget target;
        })
        checkedNextestPackageGroups);

      nextestCheckFor = target:
        pkgs.runCommand "gammaloop-nextest-${target.name}" {
          nativeBuildInputs = [
            ciToolchain
            pkgs.cargo-nextest
            pkgs.form
            pkgs.gcc
            nextestFailureSummary
          ] ++ lib.optionals (nextestUsesPythonModule target) [nextestPython];
          CC = nixCc;
          CXX = nixCxx;
          "${cargoLinkerVar}" = nixCc;
          LD_LIBRARY_PATH = runtimeLibPath;
          DYLD_LIBRARY_PATH = runtimeLibPath;
          NEXTEST_SHOW_PROGRESS = "counter";
          RUST_BACKTRACE = "1";
          RUST_LIB_BACKTRACE = "1";
          SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
        } (''
          if [ -z "''${SYMBOLICA_LICENSE:-}" ]; then
            echo "Missing SYMBOLICA_LICENSE environment variable" >&2
            exit 1
          fi

          mkdir -p /build/source
          cp -R ${nextestRuntimeSrcFor target}/. /build/source/
          chmod -R u+w /build/source
          cd /build/source
          # Workspace-root discovery checks these directories even for test
          # targets that do not consume any files from them.
          mkdir -p tests/resources examples/cli
          ${workspaceMissingCargoTargetsScript}
        '' + lib.optionalString (nextestUsesPythonModule target) ''
          export PYO3_PYTHON=${nextestPython}/bin/python3
          export PYTHON=${nextestPython}/bin/python3
          export PYTHONPATH=${gammaloop-python-module}/${pythonSitePackages}:${nextestPython}/${pythonSitePackages}
        '' + ''
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

          mkdir -p target/nextest
          set +e
          status=0
          ${lib.concatMapStringsSep "\n" (package: ''
              rm -f ${lib.escapeShellArg nextestJunitPath}
              cargo nextest run \
                --archive-file ${nextestBinarySetForTarget target}/${nextestArchiveNameFor target package} \
                --workspace-remap . \
                ${nextestBaseExtraArgs}
              package_status=$?
              nextest-failure-summary ${lib.escapeShellArg nextestJunitPath} || true
              if [ "$package_status" -ne 0 ] && [ "$status" -eq 0 ]; then
                status="$package_status"
              fi
            '')
            target.packages}
          mkdir -p "$out"
          exit "$status"
        '');

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
            runtimeInputs = [cargoArtifacts];
          }
          {
            runnerAttr = "nix-ci-check-gammaloop-nextest";
            checkAttr = "gammaloop-nextest";
            runtimeInputs = [
              nextestBinarySetAggregate
              gammaloop-python-module
            ];
          }
        ]
        ++ map (target: {
          runnerAttr = "nix-ci-check-gammaloop-nextest-${target.name}";
          checkAttr = "gammaloop-nextest-${target.name}";
          runtimeInputs =
            [(nextestBinarySetForTarget target)]
            ++ lib.optionals (nextestUsesPythonModule target) [gammaloop-python-module];
        })
        checkedNextestPackageGroups;

      impureCheckRunnerPackages = lib.listToAttrs (map (target: {
          name = target.runnerAttr;
          value = pkgs.writeShellApplication {
            name = target.runnerAttr;
            # Retain the pure test inputs even when NixCI skips this runner as
            # cached, so the in-repo test does not rebuild them.
            runtimeInputs = [pkgs.nix] ++ target.runtimeInputs;
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

      allChecks =
        {
          # Keep existing check names for CI compatibility.
          gammaloop = gammaloop-cli;

          gammaloop-check = workspaceCargoCheck;

          gammaloop-clippy = workspaceClippyCheck;

          gammaloop-doc = workspaceDocCheck;

          gammaloop-doctest = workspaceDoctestCheck;

          gammaloop-fmt = craneLib.cargoFmt {
            src = workspaceFmtSrc;
            pname = "gammaloop-workspace";
            inherit (apiMeta) version;
          };

          gammaloop-guppy-workspace-graph = guppyWorkspaceGraphCheck;

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
        // nextestChecks
        // {
          gammaloop-nextest = nextestAggregate;
        };

      # Hestia builds the binary producers before evaluating this consumer
      # matrix. Omit those producers and the aggregate nextest check so each
      # nextest execution remains an independent row without racing its leaves.
      hestiaChecks = builtins.removeAttrs allChecks (
        [
          "gammaloop-nextest"
          "gammaloop-nextest-binaries"
        ]
        ++ map (target: "gammaloop-nextest-binaries-${target.name}") checkedNextestPackageGroups
      );
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
          "alphal00p-docs-cargo-artifacts" = alphal00pDocsCargoArtifacts;
          "alphal00p-docs-pages" = alphal00pDocsPages;
          "alphal00p-docs-snapshot-fixture" = alphal00pDocsSnapshotFixture;
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
        // lib.optionalAttrs (!pkgs.stdenv.hostPlatform.isDarwin) {
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
        default = mkDevShell craneLib [clinnet-cli] "-C linker=${nixCc}";
        full = mkDevShell craneLib [clinnet-cli rscls] "-C linker=${nixCc}";
        clinnet = mkDevShell craneLib [clinnet-cli] "-C linker=${nixCc}";
        # The native development shell deliberately does not pull every target
        # standard library. Use this shell for Tydenso and other Wasm builds.
        wasm = mkDevShell wasmCraneLib [clinnet-cli pkgs.binaryen] null;
      };
    });
}
