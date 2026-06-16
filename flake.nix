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
      workspacePrebuildPackage = "gammaloop-ci-prebuild";
      workspacePrebuildPackageDir = "crates/${workspacePrebuildPackage}";

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

      dummyProcMacroCargoTarget = pkgs.writeText "crane-dummy-proc-macro-cargo-target.rs" ''
        #![allow(clippy::all)]
        #![allow(dead_code)]
      '';

      normalizeWorkspaceHackBuildScriptTimestampScriptFor = prefix: ''
        if [ -e ${prefix}crates/${workspaceHackPackage}/build.rs ]; then
          if [ -e ${prefix}crates/${workspaceHackPackage}/Cargo.toml ]; then
            touch -d @0 ${prefix}crates/${workspaceHackPackage}/Cargo.toml
          fi
          if [ -e ${prefix}crates/${workspaceHackPackage}/src/lib.rs ]; then
            touch -d @0 ${prefix}crates/${workspaceHackPackage}/src/lib.rs
          fi
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

      workspaceTestSupportComponentRepresentatives =
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

      workspaceTestComponentDependencyClosureFor = representative:
        sortedUnique (lib.filter (dependencyRepresentative: dependencyRepresentative != representative) (
          map (entry: entry.key) (builtins.genericClosure {
            startSet = [{key = representative;}];
            operator = entry: map (dependencyRepresentative: {key = dependencyRepresentative;}) (
              workspaceTestComponentDependencyRepresentativesFor entry.key
            );
          })
        ));

      workspaceTestSyntheticConsumerExcludedPackages =
        workspaceFeatureUnificationExcludedPackages ++ [
          "gammaloop-api"
          "gammaloop-integration-tests"
          "gammalooprs"
        ];

      workspaceTestComponentSyntheticConsumerPackagesFor = representative:
        sortedUnique (lib.filter (
            package: let
              consumerRepresentative = workspaceTestComponentRepresentativeFor package;
              consumerDependencyRepresentatives =
                workspaceTestComponentDependencyRepresentativesFor consumerRepresentative;
              dependencyIsCoveredByIntermediate =
                builtins.any (
                  dependencyRepresentative:
                    dependencyRepresentative != representative
                    && builtins.elem representative (
                      workspaceTestComponentDependencyClosureFor dependencyRepresentative
                    )
                )
                consumerDependencyRepresentatives;
            in
              !(builtins.elem package workspaceTestComponentMembers.${representative})
              && !(builtins.elem package workspaceTestSyntheticConsumerExcludedPackages)
              && builtins.any (
                dependency: workspaceTestComponentRepresentativeFor dependency == representative
              )
              (workspaceTestDependencyNamesFor package)
              && !dependencyIsCoveredByIntermediate
          )
          workspaceMemberPackages);

      workspaceTestComponentSyntheticConsumerSourcePackageNamesFor = representative:
        sortedUnique (
          lib.concatMap workspaceTestSourcePackageNamesFor
          (workspaceTestComponentSyntheticConsumerPackagesFor representative)
        );

      workspaceSourcePackageNamesFor = package: dependencies:
        sortedUnique ([package] ++ dependencies);

      workspaceNormalSourcePackageNamesFor = package:
        workspaceSourcePackageNamesFor package (workspaceResolvedDependencyNamesFor package);

      workspaceTestSourcePackageNamesFor = package:
        workspaceSourcePackageNamesFor package (workspaceResolvedTestDependencyNamesFor package);

      workspaceDirectNormalSourcePackageNamesFor = package:
        workspaceSourcePackageNamesFor package (workspaceDependencyNamesFor package);

      workspaceDirectTestSourcePackageNamesFor = package:
        workspaceSourcePackageNamesFor package (workspaceTestDependencyNamesFor package);

      workspaceTestComponentSourcePackageNamesFor = representative:
        sortedUnique (lib.concatMap workspaceTestSourcePackageNamesFor workspaceTestComponentMembers.${representative});

      workspacePackageExtraSourceRoots = {
        "gammaloop-api" = [
          "assets"
        ];
        gammalooprs = [
          "assets"
        ];
        "gammaloop-integration-tests" = [
          "assets"
          "examples/cli"
        ];
      };

      workspacePackageExtraSourceRootsForSourcePackages = sourcePackages:
        sortedUnique (lib.concatMap (sourcePackage: workspacePackageExtraSourceRoots.${sourcePackage} or []) sourcePackages);

      workspacePackageExtraFilesetsForSourcePackages = sourcePackages:
        map (sourceRoot: workspaceRoot + "/${sourceRoot}") (workspacePackageExtraSourceRootsForSourcePackages sourcePackages);

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

      workspacePackageSrcForSourcePackages = package: sourcePackages:
        lib.fileset.toSource {
          root = workspaceRoot;
          fileset = lib.fileset.unions (
            workspaceDependencyManifestFiles
            ++ workspaceDependencyBuildScripts
            ++ map (
              sourcePackage:
                workspaceRoot + "/${workspaceMemberPackageDirs.${sourcePackage}}"
            ) sourcePackages
            ++ [
              ./tests/resources/fjcore
            ]
            ++ (workspacePackageExtraFilesetsForSourcePackages (sortedUnique ([package] ++ sourcePackages)))
          );
        };

      workspacePackageSrcFor = package:
        workspacePackageSrcForSourcePackages package (workspaceNormalSourcePackageNamesFor package);

      workspaceTestPackageSrcFor = package:
        workspacePackageSrcForSourcePackages package (workspaceTestSourcePackageNamesFor package);

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

      workspacePackageHasBinTarget = package: let
        manifest = workspaceManifestFor workspaceMemberPackageDirs.${package};
        packageDir = workspaceMemberPackageDirs.${package};
      in
        (manifest.bin or []) != [] || builtins.pathExists (workspaceRoot + "/${packageDir}/src/main.rs");

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

      workspaceDummyCargoTargetPathsForSourcePackages = sourcePackages:
        lib.filter (
          path:
            !builtins.any (
              sourcePackage:
                lib.hasPrefix "${workspaceMemberPackageDirs.${sourcePackage}}/" path
            ) sourcePackages
        )
        workspaceCargoTargetRelPaths;

      workspaceDummyCargoTargetsScriptForSourcePackages = sourcePackages:
        ''
          ${lib.concatMapStringsSep "\n" (path: ''
            if [ ! -e ${lib.escapeShellArg path} ]; then
              install -D -m 0644 ${workspaceDummyCargoTargetForPath path} ${lib.escapeShellArg path}
            fi
          '')
          (workspaceDummyCargoTargetPathsForSourcePackages sourcePackages)}

          ${normalizeWorkspaceHackBuildScriptTimestampScript}
        '';

      workspaceDummyCargoTargetsScriptFor = package:
        workspaceDummyCargoTargetsScriptForSourcePackages (workspaceNormalSourcePackageNamesFor package);

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
            source = workspaceRoot + "/${packageDir}";
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
        "linnet-py"
        "spynso3"
      ];

      workspaceIncomingNormalDependencyFeaturesFor = dependency:
        sortedUnique (lib.concatMap (
            package: workspaceGraph.normal_dependency_features.${package}.${dependency} or []
          ) (lib.subtractLists workspaceFeatureUnificationExcludedPackages workspaceMemberPackages));

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
                map (
                  feature:
                    if lib.hasInfix "/" feature
                    then feature
                    else "${package}/${feature}"
                ) (featuresFor package)
            )
            packages;
        in
          cargoFeatureArgs (sortedUnique features);

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
      testSupportFeatureAnchorPackageFor = representative: "gammaloop-ci-test-support-${representative}-features";
      testSupportFeatureAnchorPackageDirFor = representative: "crates/${testSupportFeatureAnchorPackageFor representative}";
      testSupportFeatureAnchorDependencyPackagesFor = sourcePackages:
        lib.filter workspacePackageHasLibTarget sourcePackages;
      testSupportFeatureAnchorCargoTomlFor = representative: sourcePackages: pkgs.writeText "${testSupportFeatureAnchorPackageFor representative}-Cargo.toml" ''
        [package]
        name = "${testSupportFeatureAnchorPackageFor representative}"
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
        '') (testSupportFeatureAnchorDependencyPackagesFor sourcePackages)}
      '';
      testSupportFeatureAnchorSourceScriptFor = representative: sourcePackages: prefix: ''
        install -D -m 0644 ${testSupportFeatureAnchorCargoTomlFor representative sourcePackages} "${prefix}${testSupportFeatureAnchorPackageDirFor representative}/Cargo.toml"
        install -D -m 0644 ${dummyCargoTarget} "${prefix}${testSupportFeatureAnchorPackageDirFor representative}/src/lib.rs"
      '';
      cargoTestSupportArgsFor = representative: packages: sourcePackages: let
        featureAnchorPackage = testSupportFeatureAnchorPackageFor representative;
        anchorPackages =
          lib.optionals (
            builtins.any (
              featurePackage:
                featurePackage != workspaceHackPackage
                && builtins.elem featurePackage workspaceGraph.symbolica_test_packages
            )
            sourcePackages
          ) [workspaceHackPackage];
        selectedPackages = sortedUnique (packages ++ anchorPackages ++ [featureAnchorPackage]);
        featureArgs = cargoQualifiedFeatureArgsFor (sortedUnique (sourcePackages ++ anchorPackages)) craneTestFeaturesFor;
      in
        lib.concatStringsSep " " (
          [
            "--offline"
          ]
          ++ map (selectedPackage: "-p ${lib.escapeShellArg selectedPackage}") selectedPackages
          ++ lib.optional (featureArgs != "") featureArgs
        );
      cranePythonExtraFeatureSets = {
        "gammaloop-api" = ["python_abi" "pyo3-extension-module"];
      };
      cranePythonFeaturesFor = package:
        sortedUnique (craneCiFeaturesFor package ++ (cranePythonExtraFeatureSets.${package} or []));
      cranePythonCargoArgs = let
        featurePackages = workspaceNormalSourcePackageNamesFor "gammaloop-api";
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
          cargoArtifacts = mergeCargoArtifacts "gammaloop-python-package-inputs" [
            cranePythonBuildArtifacts
          ];
          doNotLinkInheritedArtifacts = true;
          pname = "gammaloop-api-python";
          src = workspacePackageSrcFor "gammaloop-api";
          cargoExtraArgs = cranePythonCargoArgs;
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
          cargoExtraArgs = cargoPackageCiArgsFor "clinnet";
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
      cargoArtifacts = buildDepsOnlyWithArtifacts (ciArgs
        // {
          pname = "gammaloop-workspace-deps";
          src = workspaceDependencySrc;
          cargoArtifacts = workspaceHackBuildArtifacts;
          cargoExtraArgs = workspacePrebuildCargoArgs;
          buildPhaseCargoCommand = "cargoWithProfile build ${workspacePrebuildCargoArgs}";
          checkPhaseCargoCommand = "";
          doCheck = false;
          stripWorkspaceArtifacts = true;
          extraDummyScript =
            ''
              ${workspacePrebuildSourceScript}

              ${lib.concatMapStringsSep "\n" (path: ''
                install -D -m 0644 ${dummyCargoTarget} "$out/${path}"
              '')
              autoCargoTargetPaths}
            '';
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
        cat > "$tmp/ci-features.json" <<'EOF'
        ${guppyFeatureMapFor craneCiFeaturesFor}
        EOF
        cat > "$tmp/test-features.json" <<'EOF'
        ${guppyFeatureMapFor craneTestFeaturesFor}
        EOF

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

          def workspace_dep_features($path_map; $kinds):
            [
              .dependencies[]
              | select((.kind // "normal") as $kind | $kinds | index($kind))
              | select(.path != null)
              | ($path_map[.path] // empty) as $package
              | {
                  key: $package,
                  value: (.features // [])
                }
            ]
            | group_by(.key)
            | map({
                key: .[0].key,
                value: (map(.value[]) | unique | sort)
              })
            | from_entries;

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
              ),
              normal_dependency_features: (
                $packages
                | map({key: .name, value: workspace_dep_features($path_map; ["normal", "build"])})
                | from_entries
              ),
              test_dependency_features: (
                $packages
                | map({key: .name, value: workspace_dep_features($path_map; ["normal", "build", "dev"])})
                | from_entries
              )
            }
        ' "$tmp/metadata.json" > "$tmp/base.json"

          jq -r '.packages[]' "$tmp/base.json" > "$tmp/workspace-packages"

        guppy_feature_args() {
          local package="$1"
          local include_dev="$2"
          local feature_file="$tmp/ci-features.json"
          local features

          if [[ "$include_dev" == "true" ]]; then
            feature_file="$tmp/test-features.json"
          fi

          features="$(jq -r --arg package "$package" '.[$package] // [] | join(",")' "$feature_file")"
          if [[ -n "$features" ]]; then
            printf '%s\n' --features "$features"
          fi
        }

        guppy_names() {
          local package="$1"
          local kind="$2"
          local include_dev="$3"
          local feature_args=()
          mapfile -t feature_args < <(guppy_feature_args "$package" "$include_dev")

          if [[ "$include_dev" == "true" ]]; then
            cargo-guppy resolve-cargo \
              --manifest-path "$manifest_path" \
              --resolver-version v2 \
              --target-platform any \
              --host-platform any \
              -k "$kind" \
              -p "$package" \
              "''${feature_args[@]}" \
              --include-dev
          else
            cargo-guppy resolve-cargo \
              --manifest-path "$manifest_path" \
              --resolver-version v2 \
              --target-platform any \
              --host-platform any \
              -k "$kind" \
              -p "$package" \
              "''${feature_args[@]}"
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
        pkgs.runCommand name {
          nativeBuildInputs = [pkgs.rsync pkgs.zstd pkgs.gnutar];
        } ''
          mkdir -p "$out/target"

          make_target_writable() {
            chmod -R u+w "$out/target"
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
              zstd -d "$artifact" --stdout | tar --no-same-permissions -x -C "$out/target"
              make_target_writable
            elif [ -d "$artifact" ]; then
              make_target_writable
              rsync -a --chmod=u+w "$artifact/" "$out/target/"
              make_target_writable
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

      workspacePackageCargoArtifactNames = package: let
        manifest = workspaceManifestFor workspaceMemberPackageDirs.${package};
        libManifest = manifest.lib or {};
        targetNameForCargo = name: lib.replaceStrings ["-"] ["_"] name;
      in
        [(targetNameForCargo (libManifest.name or manifest.package.name))];

      stripSelectedWorkspaceCargoArtifactsScript = packagesToStrip: let
        stripPackageScript = package: let
          artifactNames = workspacePackageCargoArtifactNames package;
          fingerprintNames = sortedUnique ([package] ++ artifactNames);
        in ''
          ${lib.concatMapStringsSep "\n" (name: ''
              rm -f target/${ciCargoProfile}/deps/${lib.escapeShellArg name}
              rm -f target/${ciCargoProfile}/deps/${lib.escapeShellArg name}.d
              rm -f target/${ciCargoProfile}/deps/${lib.escapeShellArg name}-*
              rm -f target/${ciCargoProfile}/deps/lib${lib.escapeShellArg name}.*
              rm -f target/${ciCargoProfile}/deps/lib${lib.escapeShellArg name}-*
              rm -f target/${ciCargoProfile}/${lib.escapeShellArg name}
              rm -f target/${ciCargoProfile}/${lib.escapeShellArg name}.d
              rm -f target/${ciCargoProfile}/lib${lib.escapeShellArg name}.*
              rm -f target/${ciCargoProfile}/lib${lib.escapeShellArg name}-*
            '')
            artifactNames}
          ${lib.concatMapStringsSep "\n" (name: ''
              rm -rf target/${ciCargoProfile}/.fingerprint/${lib.escapeShellArg name}-[0-9a-f]*
              rm -rf target/${ciCargoProfile}/build/${lib.escapeShellArg name}-[0-9a-f]*
            '')
            fingerprintNames}
        '';
      in ''
        # Cargo feature anchors can compile dummy workspace units. Keep only
        # artifacts that are intended to be reused by later derivations.
        ${lib.concatMapStringsSep "\n" stripPackageScript packagesToStrip}
      '';

      stripWorkspaceCargoArtifactsScript = {
        preserveProcMacros ? false,
        preservedPackages ? [],
      }: let
        packagesToStrip =
          lib.filter (
            package:
              !(builtins.elem package preservedPackages)
              && !(preserveProcMacros && workspacePackageIsProcMacro package)
          )
          workspacePackages;
      in
        stripSelectedWorkspaceCargoArtifactsScript packagesToStrip;

      # Crane's buildDepsOnly deliberately clears cargoArtifacts. This local
      # variant keeps the same dummy-source behavior while allowing per-crate
      # dependency archives to build incrementally from earlier archives.
      buildDepsOnlyWithArtifacts = {
        cargoBuildCommand ? "cargoWithProfile build",
        cargoCheckCommand ? "cargoWithProfile check",
        cargoExtraArgs ? "--locked",
        cargoTestCommand ? "cargoWithProfile test",
        cargoTestExtraArgs ? "--no-run",
        preBuildWorkspaceArtifactStripPackages ? [],
        preservedWorkspaceArtifactPackages ? [],
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
          "preBuildWorkspaceArtifactStripPackages"
          "preserveWorkspaceProcMacroArtifacts"
          "preservedWorkspaceArtifactPackages"
          "stripWorkspaceArtifacts"
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
        stripWorkspaceArtifactsScriptText =
          lib.optionalString (args.stripWorkspaceArtifacts or false) (stripWorkspaceCargoArtifactsScript {
            preserveProcMacros = args.preserveWorkspaceProcMacroArtifacts or false;
            preservedPackages = preservedWorkspaceArtifactPackages;
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

                  tar --sort=name --mtime=@0 --owner=0 --group=0 --numeric-owner -C "$tmp/target" -cf - . \
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
          postPatch = workspaceDummyCargoTargetsScriptFor workspaceHackPackage;
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
              lib.filter workspacePackageIsProcMacro (workspaceResolvedDependencyNamesFor package);
          in
            buildDepsOnlyWithArtifacts ((builtins.removeAttrs ciArgs ["src"])
              // {
                cargoArtifacts = mergeCargoArtifactsOrNull "gammaloop-crate-${package}-dependency-deps-inputs" (
                  [cargoArtifacts]
                  ++ rootPackageDependencyArtifactsFor package
                  ++ map dependencyArtifactFor dependencySourcePackages
                );
                pname = "gammaloop-crate-${package}-dependency-deps";
                src = workspacePackageSrcForSourcePackages package dependencySourcePackages;
                buildPhaseCargoCommand = "cargoWithProfile build ${cargoPackageDependencyModeArgsFor package}";
                checkPhaseCargoCommand = "";
                doCheck = false;
                stripWorkspaceArtifacts = true;
                inherit preservedWorkspaceArtifactPackages;
                extraDummyScript = ''
                  ${workspaceDependencyDummyCargoTargetsScriptFor package}
                  ${workspaceConsumerSourceScriptFor package dependencySourcePackages}
                '';
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
            preservedWorkspaceArtifactPackages =
              sortedUnique ([package] ++ workspaceResolvedDependencyNamesFor package);
          in
            buildDepsOnlyWithArtifacts ((builtins.removeAttrs ciArgs ["src"])
              // {
                cargoArtifacts = mergeCargoArtifactsOrNull "gammaloop-crate-${package}-dependency-inputs" (
                  [cranePackageDependencyModeDependencyArtifacts.${package}]
                  ++ map dependencyArtifactFor (lib.filter (sourcePackage: sourcePackage != package) sourcePackages)
                );
                pname = "gammaloop-crate-${package}-dependency";
                src = workspacePackageSrcForSourcePackages package sourcePackages;
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
                postPatch = workspaceDummyCargoTargetsScriptForSourcePackages sourcePackages;
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
                cargoArtifacts = mergeCargoArtifacts "gammaloop-crate-${package}-deps" (
                  [
                    cargoArtifacts
                    cranePackageDependencyArtifacts.${package}
                  ]
                );
                doNotLinkInheritedArtifacts = true;
                pname = "gammaloop-crate-${package}";
                src = workspacePackageSrcFor package;
                cargoExtraArgs =
                  cargoPackageCiArgsFor package
                  + lib.optionalString (package == "gammaloop-api") " --lib --bins";
                postPatch = workspaceDummyCargoTargetsScriptFor package;
              })));

      gammaloopApiPackageArtifacts = mergeCargoArtifacts "gammaloop-api-package-artifacts" [
        cranePackageBuildArtifacts."gammaloop-api"
      ];

      craneTestSupportArtifacts = lib.fix (self:
        lib.genAttrs workspaceTestComponentRepresentatives (representative:
          if representative == workspaceHackPackage
          then workspaceHackBuildArtifacts
          else let
            componentPackages = workspaceTestComponentMembers.${representative};
            dependencyComponents = workspaceTestComponentDependencyRepresentativesFor representative;
            sourcePackages = workspaceTestComponentSourcePackageNamesFor representative;
            syntheticConsumerPackages = workspaceTestComponentSyntheticConsumerPackagesFor representative;
            syntheticConsumerSourcePackages = workspaceTestComponentSyntheticConsumerSourcePackageNamesFor representative;
            selectedPackages = sortedUnique (componentPackages ++ syntheticConsumerPackages);
            selectedBinaryPackages = lib.filter workspacePackageHasBinTarget componentPackages;
            featureSourcePackages = sortedUnique (sourcePackages ++ syntheticConsumerSourcePackages);
            supportCargoCommand =
              lib.concatStringsSep "\n" (
                lib.optional (selectedBinaryPackages != [])
                  "cargoWithProfile build ${cargoTestSupportArgsFor representative selectedBinaryPackages featureSourcePackages} --bins"
                ++ [
                  "cargoWithProfile test --no-run ${cargoTestSupportArgsFor representative selectedPackages featureSourcePackages}"
                ]
              );
          in
            buildDepsOnlyWithArtifacts ((builtins.removeAttrs ciArgs ["src"])
              // {
                cargoArtifacts = mergeCargoArtifactsOrNull "gammaloop-crate-test-support-${representative}-inputs" (
                  [workspaceHackBuildArtifacts]
                  ++ map (dependencyComponent: self.${dependencyComponent}) dependencyComponents
                );
                pname = "gammaloop-crate-test-support-${representative}";
                dummySrc = workspacePackageSrcForSourcePackages representative sourcePackages;
                buildPhaseCargoCommand = supportCargoCommand;
                checkPhaseCargoCommand = "";
                doCheck = false;
                stripWorkspaceArtifacts = syntheticConsumerPackages != [];
                preservedWorkspaceArtifactPackages = sourcePackages;
                postPatch = ''
                  ${workspaceDummyCargoTargetsScriptForSourcePackages sourcePackages}
                  ${testSupportFeatureAnchorSourceScriptFor representative featureSourcePackages ""}
                '';
              })));

      craneTestBinaryArtifacts = lib.genAttrs workspacePackages (package:
        if package == workspaceHackPackage
        then workspaceHackBuildArtifacts
        else craneTestSupportArtifacts.${workspaceTestComponentRepresentativeFor package});

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
            ) (lib.filter (sourcePackage: sourcePackage != "gammaloop-api") (workspaceNormalSourcePackageNamesFor "gammaloop-api"))
          );
          pname = "gammaloop-api-python";
          src = workspacePackageSrcForSourcePackages "gammaloop-api" (
            lib.filter (sourcePackage: sourcePackage != "gammaloop-api") (workspaceNormalSourcePackageNamesFor "gammaloop-api")
          );
          buildPhaseCargoCommand = "cargoWithProfile build ${cranePythonCargoArgs}";
          checkPhaseCargoCommand = "";
          doCheck = false;
          preBuildWorkspaceArtifactStripPackages = ["gammaloop-api"];
          stripWorkspaceArtifacts = true;
          preservedWorkspaceArtifactPackages = workspaceResolvedDependencyNamesFor "gammaloop-api";
          extraDummyScript = workspaceDependencyDummyCargoTargetsScriptFor "gammaloop-api";
        });

      cranePythonBuildArtifacts = craneLib.cargoBuild (ciArgs
        // {
          cargoArtifacts = mergeCargoArtifacts "gammaloop-python-build-inputs" [
            cranePythonDependencyArtifacts
          ];
          doNotLinkInheritedArtifacts = true;
          pname = "gammaloop-api-python-build";
          src = workspacePackageSrcFor "gammaloop-api";
          cargoExtraArgs = cranePythonCargoArgs;
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

      craneTestSupportOutputs = lib.listToAttrs (map (representative: {
          name = "crate-test-support-${representative}";
          value = craneTestSupportArtifacts.${representative};
        })
        workspaceTestSupportComponentRepresentatives);

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
          name = "clinnet";
          packages = ["clinnet"];
        }
        {
          name = "linnet";
          packages = [
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
      nextestUsesIntegrationTests = target: builtins.elem "gammaloop-integration-tests" target.packages;
      nextestTargetExtraFeaturesFor = target: package:
        if (target ? extraFeatures) && builtins.hasAttr package target.extraFeatures
        then target.extraFeatures.${package}
        else [];
      nextestUsesPythonModule = target:
        builtins.any (
          package: builtins.elem "python-api-tests" (nextestTargetExtraFeaturesFor target package)
        ) target.packages;
      nextestSourcePackagesFor = target:
        sortedUnique (
          target.packages
          ++ lib.concatMap workspaceResolvedTestDependencyNamesFor target.packages
        );
      nextestSrcFor = target: let
        sourcePackages = nextestSourcePackagesFor target;
      in
        lib.fileset.toSource {
          root = workspaceRoot;
          fileset = lib.fileset.unions (
            workspaceDependencyManifestFiles
            ++ workspaceDependencyBuildScripts
            ++ map (
              sourcePackage:
                workspaceRoot + "/${workspaceMemberPackageDirs.${sourcePackage}}"
            )
            sourcePackages
            ++ [
              ./.config/nextest.toml
              ./tests/resources
              ./examples/api
              ./examples/cli
            ]
            ++ (workspacePackageExtraFilesetsForSourcePackages sourcePackages)
          );
        };

      nextestFeatureArgsFor = target:
        cargoQualifiedFeatureArgsFor (nextestFeaturePackagesFor target) (
          package:
            sortedUnique (
              craneCiFeaturesFor package
              ++ nextestTargetExtraFeaturesFor target package
            )
        );

      nextestCargoArgsFor = target:
        lib.concatStringsSep " " (
          [
            "--offline"
          ]
          ++ map (package: "-p ${lib.escapeShellArg package}") (nextestSelectedPackagesFor target)
          ++ lib.optional (nextestFeatureArgsFor target != "") (nextestFeatureArgsFor target)
        );

      nextestArchiveNameFor = target: "gammaloop-nextest-${target.name}.tar.zst";
      nextestArchiveInputPackagesFor = target:
        sortedUnique target.packages;
      nextestAnchorPackagesFor = target:
        lib.optionals (
          builtins.any (
            package: builtins.elem package workspaceGraph.symbolica_test_packages
          ) (nextestSourcePackagesFor target)
        ) [workspaceHackPackage];
      nextestSelectedPackagesFor = target:
        sortedUnique (target.packages ++ nextestAnchorPackagesFor target ++ [(nextestFeatureAnchorPackageFor target)]);
      nextestFeaturePackagesFor = target:
        sortedUnique (
          nextestAnchorPackagesFor target
          ++ nextestSourcePackagesFor target
        );
      nextestFeatureAnchorPackageFor = target: "gammaloop-ci-nextest-${target.name}-features";
      nextestFeatureAnchorPackageDirFor = target: "crates/${nextestFeatureAnchorPackageFor target}";
      nextestFeatureAnchorDependencyPackagesFor = target:
        lib.filter workspacePackageHasLibTarget (nextestFeaturePackagesFor target);
      nextestFeatureAnchorCargoTomlFor = target: pkgs.writeText "${nextestFeatureAnchorPackageFor target}-Cargo.toml" ''
        [package]
        name = "${nextestFeatureAnchorPackageFor target}"
        version = "0.1.0"
        edition = "2024"
        publish = false

        [lib]
        path = "src/lib.rs"

        [dependencies]
        ${lib.concatMapStringsSep "\n" (package: let
          features = lib.filter (feature: !lib.hasInfix "/" feature) (
            sortedUnique (craneCiFeaturesFor package ++ nextestTargetExtraFeaturesFor target package)
          );
          featureEntry = lib.optionalString (features != []) ", features = ${builtins.toJSON features}";
        in ''
          ${package} = { path = "${workspacePrebuildDependencyPathFor package}"${featureEntry} }
        '') (nextestFeatureAnchorDependencyPackagesFor target)}
      '';
      nextestFeatureAnchorSourceScriptFor = target: prefix: ''
        install -D -m 0644 ${nextestFeatureAnchorCargoTomlFor target} "${prefix}${nextestFeatureAnchorPackageDirFor target}/Cargo.toml"
        install -D -m 0644 ${dummyCargoTarget} "${prefix}${nextestFeatureAnchorPackageDirFor target}/src/lib.rs"
      '';
      testSupportArtifactForPackage = package:
        if package == workspaceHackPackage
        then workspaceHackBuildArtifacts
        else craneTestSupportArtifacts.${workspaceTestComponentRepresentativeFor package};

      nextestPackageTestBinaryArtifactFor = target: package: let
        packageTarget =
          target
          // {
            name = "${target.name}-${package}";
            packages = [package];
          };
        sourcePackages = nextestSourcePackagesFor packageTarget;
      in
        buildDepsOnlyWithArtifacts ((builtins.removeAttrs ciArgs ["src"])
          // {
            cargoArtifacts = mergeCargoArtifacts "gammaloop-nextest-binaries-${target.name}-${package}-test-inputs" (
              [cargoArtifacts]
              ++ [(testSupportArtifactForPackage package)]
            );
            pname = "gammaloop-nextest-binaries-${target.name}-${package}-test";
            dummySrc = nextestSrcFor packageTarget;
            buildPhaseCargoCommand = ''
              if [ -d target ]; then
                chmod -R u+w target
              fi
              cargoWithProfile test --no-run ${nextestCargoArgsFor packageTarget}
            '';
            checkPhaseCargoCommand = "";
            doCheck = false;
            postPatch = ''
              ${workspaceDummyCargoTargetsScriptForSourcePackages sourcePackages}
              ${nextestFeatureAnchorSourceScriptFor packageTarget ""}
            '';
          } // lib.optionalAttrs (nextestUsesPythonModule packageTarget) {
            nativeBuildInputs = (ciArgs.nativeBuildInputs or []) ++ [nextestPython];
            PYO3_PYTHON = "${nextestPython}/bin/python3";
            PYTHON = "${nextestPython}/bin/python3";
            PYTHONPATH = "${gammaloop-python-module}/${pythonSitePackages}:${nextestPython}/${pythonSitePackages}";
          });

      nextestArchiveCargoArtifactsFor = target:
        mergeCargoArtifacts "gammaloop-nextest-binaries-${target.name}-inputs" (
          [
            cargoArtifacts
          ]
          ++ map (nextestPackageTestBinaryArtifactFor target) target.packages
        );

      nextestArchiveFor = target:
        craneLib.mkCargoDerivation (ciArgs
          // {
            pname = "gammaloop-nextest-binaries-${target.name}";
            src = nextestSrcFor target;
            cargoArtifacts = nextestArchiveCargoArtifactsFor target;
            doNotLinkInheritedArtifacts = true;
            doCheck = false;
            doInstallCargoArtifacts = false;
            nativeBuildInputs =
              (ciArgs.nativeBuildInputs or [])
              ++ [pkgs.cargo-nextest pkgs.form]
              ++ lib.optionals (nextestUsesPythonModule target) [nextestPython];
            postPatch = ''
              ${workspaceDummyCargoTargetsScriptForSourcePackages (nextestSourcePackagesFor target)}
              ${nextestFeatureAnchorSourceScriptFor target ""}
            '';
            buildPhaseCargoCommand = ''
              mkdir -p "$out"
              if [ -d target ]; then
                chmod -R u+w target
              fi
              cargo nextest --version
              cargo nextest archive \
                --cargo-profile ${ciCargoProfile} \
                ${nextestCargoArgsFor target} \
                ${nextestFilterFor target} \
                --profile ${nextestProfile} \
                --archive-file "$out/${nextestArchiveNameFor target}"
            '';
            checkPhaseCargoCommand = "";
            installPhaseCommand = "";
          } // lib.optionalAttrs (nextestUsesPythonModule target) {
            PYO3_PYTHON = "${nextestPython}/bin/python3";
            PYTHON = "${nextestPython}/bin/python3";
            PYTHONPATH = "${gammaloop-python-module}/${pythonSitePackages}:${nextestPython}/${pythonSitePackages}";
          });

      nextestBinarySets = lib.listToAttrs (map (target: {
          name = "gammaloop-nextest-binaries-${target.name}";
          value = nextestArchiveFor target;
        })
        checkedNextestPackageGroups);

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
          cp -R ${nextestSrcFor target}/. /build/source/
          chmod -R u+w /build/source
          cd /build/source
          ${workspaceDummyCargoTargetsScriptForSourcePackages (nextestSourcePackagesFor target)}
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
          RUST_MIN_STACK=''${RUST_MIN_STACK:-33554432} cargo nextest run \
            --archive-file ${nextestBinarySetForTarget target}/${nextestArchiveNameFor target} \
            --workspace-remap . \
            ${nextestBaseExtraArgs}
          status=$?
          nextest-failure-summary ${lib.escapeShellArg nextestJunitPath} || true
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
          inherit
            cargoArtifacts
            gammaloopApiPackageArtifacts
            workspaceBuildArtifacts;
          "nix-ci-passed" = nixCiPassed;
        }
        // cranePackageDependencyOutputs
        // cranePackageOutputs
        // craneTestSupportOutputs
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
