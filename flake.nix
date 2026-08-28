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

  outputs =
    {
      self,
      nixpkgs,
      crane,
      fenix,
      flake-utils,
      ...
    }:
    let
      # Nixpkgs unstable no longer supports x86_64-darwin, so exposing that
      # system makes `nix flake show` fail before NixCI can list jobs.
      supportedSystems = [
        "x86_64-linux"
        "aarch64-linux"
        "aarch64-darwin"
      ];
    in
    flake-utils.lib.eachSystem supportedSystems (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        inherit (pkgs) lib;

        typst015 =
          assert lib.assertMsg (
            pkgs.typst.version == "0.15.0"
          ) "the documentation build requires Typst 0.15.0, but nixpkgs provides ${pkgs.typst.version}";
          pkgs.typst;

        docsTypst = typst015.withPackages (
          typstPackages: with typstPackages; [
            cetz_0_5_1
            mitex_0_2_6
            tidy_0_4_3
          ]
        );

        docsFontPath = "${pkgs.roboto}/share/fonts/truetype";

        nixCiBarrierRevision =
          if self ? dirtyRev then
            self.dirtyRev
          else if self ? rev then
            self.rev
          else if self ? narHash then
            self.narHash
          else
            "local";
        # NixCI memoizes successful top-level derivations across commits without
        # re-realizing their closures. Salt only this zero-copy scheduling
        # wrapper so each commit primes the stable artifact in the shared cache.
        nixCiArtifactBarrier =
          name: artifact:
          pkgs.runCommand "nix-ci-artifact-barrier-${name}"
            {
              NIX_CI_BARRIER_REVISION = nixCiBarrierRevision;
            }
            ''
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

        craneLib = baseCraneLib.overrideToolchain ciToolchain;

        wasmTarget = "wasm32-unknown-unknown";

        wasmToolchain = fenix.packages.${system}.combine [
          (stableToolchain.withComponents [
            "cargo"
            "rust-std"
            "rustc"
          ])
          fenix.packages.${system}.targets.${wasmTarget}.stable.rust-std
        ];

        wasmCraneLib = baseCraneLib.overrideToolchain wasmToolchain;

        workspaceRoot = ./.;
        workspaceGraph = builtins.fromJSON (builtins.readFile ./nix/ci-workspace-graph.json);
        workspaceHackPackage = "gammaloop-workspace-hack";
        workspacePrebuildPackage = "gammaloop-ci-prebuild";
        workspacePrebuildPackageDir = "crates/${workspacePrebuildPackage}";

        cargoSources = craneLib.fileset.commonCargoSources workspaceRoot;
        repositoryDocumentationSources = lib.fileset.fileFilter (
          file:
          let
            name = lib.toLower file.name;
          in
          lib.any (extension: lib.hasSuffix ".${extension}" name) [
            "typ"
            "md"
            "markdown"
            "mdown"
            "mkd"
            "mdx"
            "html"
            "htm"
            "xhtml"
            "shtml"
            "rst"
            "rest"
            "adoc"
            "asciidoc"
            "org"
          ]
          || file.name == "LICENSE"
        ) workspaceRoot;

        documentationDeveloperScopePaths = lib.unique (
          map (scope: scope.path) (
            lib.concatMap (section: lib.concatMap (note: note.scope or [ ]) (section.note or [ ])) (
              (builtins.fromTOML (builtins.readFile ./docs/developers.toml)).section or [ ]
            )
          )
        );
        # Developer notes verify these source digests during the hermetic site
        # build, including scopes such as flake.nix that are not Cargo or prose.
        documentationDeveloperScopeSources = map (
          path: workspaceRoot + "/${path}"
        ) documentationDeveloperScopePaths;

        cargoVendorDir = craneLib.vendorCargoDeps {
          cargoLock = ./Cargo.lock;
          overrideVendorGitCheckout =
            packages: drv:
            if lib.any (package: package.name == "symbolica") packages then
              drv.overrideAttrs (old: {
                postInstall = (old.postInstall or "") + ''
                  for crate in ${
                    lib.concatMapStringsSep " " (
                      package: lib.escapeShellArg "${package.name}-${package.version}"
                    ) packages
                  }; do
                    if [ -d "$out/$crate" ]; then
                      mkdir -p "$out/$crate/.git"
                      printf 'ref: refs/heads/nix-vendor\n' > "$out/$crate/.git/HEAD"
                    fi
                  done
                '';
              })
            else
              drv;
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

        documentationSrc = lib.fileset.toSource {
          root = workspaceRoot;
          fileset = lib.fileset.unions (
            documentationDeveloperScopeSources
            ++ [
              cargoSources
              nonCargoBuildSources
              repositoryDocumentationSources
              ./docs
              ./scripts/check-docs-html.py
              ./scripts/render-docs-svg-assets.sh
              ./scripts/update-docs-pages.sh
              ./examples/api/python
              ./examples/cli/aa_aa/2L/graphs/GL00.dot
              ./examples/cli/aa_aa/2L/graphs/GL08.dot
              ./examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL000.dot
              ./examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL150.dot
              ./examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL300.dot
              ./examples/cli/aa_aa/3L/aa_aa_generation_timings.csv
              ./examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_unfiltered_pre_network.sym
              ./examples/cli/BNL/profiling/bnl_scalar_alias_captures.ansi.txt
              ./examples/cli/gg_hhh/3L/3L_graph.dot
              ./tests/resources/graphs/double_triangle.dot
              ./tests/resources/graphs/gghhh.dot
              ./tests/resources/graphs/qqx_aaa_pentabox_user_numerator.dot
              ./tests/resources/graphs/uv_tests/ad_ad_1L_gluon.dot
              ./tests/resources/graphs/uv_tests/epem_a_bbx.dot
              ./tests/resources/graphs/epemttbar.dot
              ./pyproject.toml
              ./crates/linnet-py/pyproject.toml
              ./crates/linnet-py/linnet_py.pyi
              ./crates/linnet-py/examples/physics_render_settings.py
              ./crates/linnet-py/tests/test_basic.py
            ]
          );
        };

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
            ./docs/api/python
            ./docs/examples.toml
            ./docs/products
            ./crates/linnet-py/linnet_py.pyi
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
            ./crates/kurvst/typst/src
            ./crates/kurvst/typst/typst.toml
            ./crates/linnest/typst/src
            ./crates/linnest/typst/typst.toml
          ];
        };

        workspaceMemberDirs =
          let
            crateEntries = builtins.readDir ./crates;
            crateMemberDirs = map (name: "crates/${name}") (
              lib.filter (
                name:
                crateEntries.${name} == "directory"
                && builtins.pathExists (workspaceRoot + "/crates/${name}/Cargo.toml")
              ) (builtins.attrNames crateEntries)
            );
          in
          crateMemberDirs ++ [ "tests" ];

        workspaceManifestFor =
          member: builtins.fromTOML (builtins.readFile (workspaceRoot + "/${member}/Cargo.toml"));

        workspaceMemberPackageDirs = lib.listToAttrs (
          map (member: {
            name = (workspaceManifestFor member).package.name;
            value = member;
          }) workspaceMemberDirs
        );

        workspaceMemberPackages = builtins.attrNames workspaceMemberPackageDirs;

        autoCargoTargetDirs = lib.concatMap (
          member:
          lib.filter (dir: builtins.pathExists (workspaceRoot + "/${dir}")) [
            "${member}/benches"
            "${member}/examples"
            "${member}/src/bin"
            "${member}/tests"
          ]
        ) workspaceMemberDirs;

        autoCargoTargetPaths = lib.sort (left: right: left < right) (
          lib.concatMap (
            dir:
            let
              entries = builtins.readDir (workspaceRoot + "/${dir}");
            in
            map (name: "${dir}/${name}") (
              lib.filter (name: entries.${name} == "regular" && lib.hasSuffix ".rs" name && name != "mod.rs") (
                builtins.attrNames entries
              )
            )
          ) autoCargoTargetDirs
        );

        workspaceExplicitCargoTargetPaths = lib.concatMap (
          member:
          let
            manifest = workspaceManifestFor member;
            targetPaths =
              targets: lib.concatMap (target: lib.optional (target ? path) "${member}/${target.path}") targets;
          in
          lib.optional ((manifest ? lib) && (manifest.lib ? path)) "${member}/${manifest.lib.path}"
          ++ targetPaths (manifest.bin or [ ])
          ++ targetPaths (manifest.example or [ ])
          ++ targetPaths (manifest.test or [ ])
          ++ targetPaths (manifest.bench or [ ])
        ) workspaceMemberDirs;

        workspaceDefaultCargoTargetPaths = lib.concatMap (member: [
          "${member}/src/lib.rs"
          "${member}/src/main.rs"
        ]) workspaceMemberDirs;

        workspaceCargoTargetRelPaths = lib.filter (path: builtins.pathExists (workspaceRoot + "/${path}")) (
          lib.sort (left: right: left < right) (
            lib.unique (
              workspaceDefaultCargoTargetPaths ++ workspaceExplicitCargoTargetPaths ++ autoCargoTargetPaths
            )
          )
        );

        workspaceCargoTargetEntrypoints = map (
          path: workspaceRoot + "/${path}"
        ) workspaceCargoTargetRelPaths;

        workspaceDependencyManifestFiles = [
          ./Cargo.lock
          ./Cargo.toml
        ]
        ++ map (member: workspaceRoot + "/${member}/Cargo.toml") workspaceMemberDirs;

        workspaceDependencyBuildScripts = lib.filter builtins.pathExists (
          [ ./build.rs ] ++ map (member: workspaceRoot + "/${member}/build.rs") workspaceMemberDirs
        );

        workspacePackageBuildScriptsForSourcePackages =
          sourcePackages:
          lib.filter builtins.pathExists (
            map (package: workspaceRoot + "/${workspaceMemberPackageDirs.${package}}/build.rs") sourcePackages
          );

        workspaceDependencySrc = lib.fileset.toSource {
          root = workspaceRoot;
          fileset = lib.fileset.unions workspaceDependencyManifestFiles;
        };

        # commonCargoSources includes every TOML file under the workspace. Keep
        # documentation metadata and test fixtures out of this reusable layer.
        documentationWorkspaceBuildSrc = lib.fileset.toSource {
          root = workspaceRoot;
          fileset = lib.fileset.unions [
            (craneLib.fileset.cargoTomlAndLock workspaceRoot)
            (craneLib.fileset.rust workspaceRoot)
            ./.cargo/config.toml
            (lib.fileset.difference nonCargoBuildSources (
              lib.fileset.unions [
                ./assets/gammalooplogo.svg
                ./assets/gammalooplogo-dark.svg
                ./assets/gammalooplogo-light.svg
              ]
            ))
          ];
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

        normalizeWorkspaceHackBuildScriptTimestampScript = normalizeWorkspaceHackBuildScriptTimestampScriptFor "";

        normalizeWorkspaceHackBuildScriptTimestampInDummySrcScript = normalizeWorkspaceHackBuildScriptTimestampScriptFor "$out/";

        workspaceDependencyNamesFor = package: workspaceGraph.normal_dependencies.${package} or [ ];

        workspaceResolvedDependencyNamesFor =
          package:
          workspaceGraph.resolved_normal_dependencies.${package} or (workspaceDependencyNamesFor package);

        workspaceResolvedTestDependencyNamesFor =
          package:
          workspaceGraph.resolved_test_dependencies.${package}
            or (workspaceResolvedDependencyNamesFor package);

        workspaceTestDependencyNamesFor =
          package: workspaceGraph.test_dependencies.${package} or (workspaceDependencyNamesFor package);

        workspaceDependencyClosureFor =
          dependencyNamesFor: package:
          sortedUnique (
            map (entry: entry.key) (
              builtins.genericClosure {
                startSet = [ { key = package; } ];
                operator = entry: map (dependency: { key = dependency; }) (dependencyNamesFor entry.key);
              }
            )
          );

        workspaceTestDependencyClosureFor = workspaceDependencyClosureFor workspaceTestDependencyNamesFor;

        workspaceTestComponentMembersFor =
          package:
          sortedUnique (
            lib.filter (
              other:
              builtins.elem other (workspaceTestDependencyClosureFor package)
              && builtins.elem package (workspaceTestDependencyClosureFor other)
            ) workspaceMemberPackages
          );

        workspaceTestComponentRepresentativeFor =
          package: builtins.head (workspaceTestComponentMembersFor package);

        workspaceTestComponentRepresentatives = sortedUnique (
          map workspaceTestComponentRepresentativeFor workspaceMemberPackages
        );

        workspaceTestDependencyComponentRepresentatives = lib.filter (
          representative: representative != workspaceHackPackage
        ) workspaceTestComponentRepresentatives;

        workspaceTestComponentMembers = lib.listToAttrs (
          map (representative: {
            name = representative;
            value = workspaceTestComponentMembersFor representative;
          }) workspaceTestComponentRepresentatives
        );

        workspaceTestComponentDependencyRepresentativesFor =
          representative:
          sortedUnique (
            lib.filter (dependencyRepresentative: dependencyRepresentative != representative) (
              map workspaceTestComponentRepresentativeFor (
                lib.concatMap workspaceTestDependencyNamesFor workspaceTestComponentMembers.${representative}
              )
            )
          );

        workspaceSourcePackageNamesFor = package: dependencies: sortedUnique ([ package ] ++ dependencies);

        workspaceNormalSourcePackageNamesFor = workspaceDependencyClosureFor workspaceResolvedDependencyNamesFor;

        workspaceTestSourcePackageNamesFor =
          package:
          sortedUnique (
            lib.concatMap workspaceNormalSourcePackageNamesFor (
              workspaceSourcePackageNamesFor package (workspaceResolvedTestDependencyNamesFor package)
            )
          );

        workspacePackageProductionExtraSourceRoots = {
          "alphal00p-docs-catalogs" = [
            "crates/gammaloop-api/src/commands/evaluate_samples.rs"
            "crates/gammaloop-api/src/commands/generate.rs"
            "crates/gammaloop-api/src/commands/import/model.rs"
            "crates/gammaloop-api/src/commands/integrate.rs"
            "crates/gammaloop-api/src/integrand_info.rs"
            "crates/gammaloop-api/src/lib.rs"
            "crates/gammaloop-api/src/session.rs"
            "crates/gammaloop-api/src/state.rs"
            "crates/gammalooprs/src/integrands/evaluation.rs"
            "crates/gammalooprs/src/integrands/mod.rs"
            "crates/gammalooprs/src/lib.rs"
            "crates/gammalooprs/src/observables/mod.rs"
            "crates/gammalooprs/src/settings/mod.rs"
            "crates/idenso/src/color/macros.rs"
            "crates/idenso/src/cook.rs"
            "crates/idenso/src/dirac/macros.rs"
            "crates/idenso/src/epsilon.rs"
            "crates/idenso/src/lib.rs"
            "crates/idenso/src/representations.rs"
            "crates/idenso/src/selective_expand.rs"
            "crates/linnet/src/half_edge.rs"
            "crates/linnet/src/half_edge/builder.rs"
            "crates/linnet/src/half_edge/involution.rs"
            "crates/linnet/src/half_edge/nodestore.rs"
            "crates/linnet/src/half_edge/subgraph.rs"
            "crates/linnet/src/half_edge/subgraph/subset.rs"
            "crates/linnet/src/half_edge/tree.rs"
            "crates/linnet/src/parser/mod.rs"
            "crates/spenso-hep-lib/src/lib.rs"
            "crates/spenso-macros/src/lib.rs"
            "crates/spenso/src/contraction.rs"
            "crates/spenso/src/network/contract.rs"
            "crates/spenso/src/network/mod.rs"
            "crates/spenso/src/network/parsing/mod.rs"
            "crates/spenso/src/structure.rs"
            "crates/spenso/src/symbolic_parallelism.rs"
            "crates/spenso/src/tensors/data/dense.rs"
            "crates/spenso/src/tensors/data/sparse.rs"
            "crates/spenso/src/tensors/parametric.rs"
            "crates/vakint/src/lib.rs"
            "crates/vakint/src/utils.rs"
          ];
          "alphal00p-docs-examples" = [
            "crates/linnet-py/pyproject.toml"
            "docs/api/python"
            "docs/examples.toml"
            "docs/products"
            "pyproject.toml"
          ];
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
          "alphal00p-docs-macros" = [
            "crates/alphal00p-docs-macros/tests/ui"
          ];
          "alphal00p-docs-python-exporter" = [
            "crates/linnet-py/linnet_py.pyi"
            "docs/api/python"
          ];
          gammalooprs = [
            "tests/resources/graphs/scalar/dod2_bubble.dot"
          ];
        };

        workspacePackageRuntimeTestExtraSourceRoots = {
          "alphal00p-docs-builder" = [
            "assets/embedded/drawing/templates/impl/physics-edge-style.typ"
            "assets/embedded/drawing/templates/layout-core.typ"
            "assets/embedded/drawing/templates/layout.typ"
            "assets/embedded/drawing/templates/physics-edge-style.typ"
            "assets/gammalooplogo-dark.svg"
            "assets/gammalooplogo-light.svg"
            "crates/idenso/CHANGELOG.typ"
            "crates/linnet/CHANGELOG.typ"
            "crates/spenso/CHANGELOG.typ"
            "docs"
            "scripts/render-docs-svg-assets.sh"
          ];
          "alphal00p-docs-catalogs" = [
            "docs/api/python"
          ];
          clinnet = [
            "assets/embedded/drawing/templates/impl/physics-edge-style.typ"
            "assets/embedded/drawing/templates/layout-core.typ"
            "assets/embedded/drawing/templates/physics-edge-style.typ"
          ];
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
            "crates/gammalooprs/src/feyngen/test.rs"
            "crates/gammalooprs/src/graph/parse/tests.rs"
            "crates/gammalooprs/src/model/test_polarization_sums.rs"
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

        workspacePackageExtraSourceRootsForSourcePackages =
          sourcePackages:
          sortedUnique (
            lib.concatMap (
              sourcePackage: workspacePackageProductionExtraSourceRoots.${sourcePackage} or [ ]
            ) sourcePackages
          );

        workspacePackageTestCompileTimeExtraSourceRootsForSourcePackages =
          sourcePackages:
          sortedUnique (
            lib.concatMap (
              sourcePackage: workspacePackageTestCompileTimeExtraSourceRoots.${sourcePackage} or [ ]
            ) sourcePackages
          );

        workspacePackageRuntimeTestExtraSourceRootsForSourcePackages =
          sourcePackages:
          sortedUnique (
            lib.concatMap (
              sourcePackage: workspacePackageRuntimeTestExtraSourceRoots.${sourcePackage} or [ ]
            ) sourcePackages
          );

        workspacePackageExtraFilesetsForSourcePackages =
          sourcePackages:
          map (sourceRoot: workspaceRoot + "/${sourceRoot}") (
            workspacePackageExtraSourceRootsForSourcePackages sourcePackages
          );

        workspacePackageTestCompileTimeExtraFilesetsForSourcePackages =
          sourcePackages:
          map (sourceRoot: workspaceRoot + "/${sourceRoot}") (
            workspacePackageTestCompileTimeExtraSourceRootsForSourcePackages sourcePackages
          );

        workspacePackageRuntimeTestExtraFilesetsForSourcePackages =
          sourcePackages:
          map (sourceRoot: workspaceRoot + "/${sourceRoot}") (
            workspacePackageRuntimeTestExtraSourceRootsForSourcePackages sourcePackages
          );

        workspacePackageOwnTestSourceRootsForSourcePackages =
          sourcePackages:
          sortedUnique (
            lib.concatMap (
              sourcePackage: workspacePackageOwnTestSourceRoots.${sourcePackage} or [ ]
            ) sourcePackages
          );

        workspacePackageOwnTestFilesetsForSourcePackages =
          sourcePackages:
          map (sourceRoot: workspaceRoot + "/${sourceRoot}") (
            workspacePackageOwnTestSourceRootsForSourcePackages sourcePackages
          );

        workspacePackageExtraSourceRestoreInDummySrcScriptFor = sourcePackages: ''
          ${lib.concatMapStringsSep "\n" (
            sourceRoot:
            let
              sourceParentDir = builtins.dirOf sourceRoot;
              source = workspaceRoot + "/${sourceRoot}";
            in
            ''
              rm -rf "$out"/${lib.escapeShellArg sourceRoot}
              mkdir -p "$out"/${lib.escapeShellArg sourceParentDir}
              cp -R --no-preserve=ownership ${source} "$out"/${lib.escapeShellArg sourceRoot}
              chmod -R u+w "$out"/${lib.escapeShellArg sourceRoot}
            ''
          ) (workspacePackageExtraSourceRootsForSourcePackages sourcePackages)}
        '';

        workspacePackageSrcForSourcePackages =
          {
            sourcePackages,
            packageSourcePackages ? [ ],
            testSourcePackages ? [ ],
            runtimeTestSourcePackages ? [ ],
            extraFilesets ? [ ],
          }:
          let
            productionSourceFilesets = map (
              sourcePackage:
              let
                packageRoot = workspaceRoot + "/${workspaceMemberPackageDirs.${sourcePackage}}";
                libraryTargetPath = workspacePackageLibTargetRelPath sourcePackage;
                testRoots = lib.filter builtins.pathExists (
                  map (root: packageRoot + "/${root}") [
                    "benches"
                    "examples"
                    "tests"
                  ]
                );
                nonLibraryTargetEntrypoints = map (path: workspaceRoot + "/${path}") (
                  lib.filter (
                    path:
                    lib.hasPrefix "${workspaceMemberPackageDirs.${sourcePackage}}/" path && path != libraryTargetPath
                  ) workspaceCargoTargetRelPaths
                );
                ownTestSources = workspacePackageOwnTestFilesetsForSourcePackages [ sourcePackage ];
                rustSources = lib.fileset.fileFilter (file: file.hasExt "rs") packageRoot;
              in
              lib.fileset.difference rustSources (
                lib.fileset.unions (testRoots ++ nonLibraryTargetEntrypoints ++ ownTestSources)
              )
            ) sourcePackages;
            packageSourceFilesets = map (
              sourcePackage:
              let
                libraryTargetPath = workspacePackageLibTargetRelPath sourcePackage;
              in
              map (path: workspaceRoot + "/${path}") (
                lib.filter (
                  path:
                  lib.hasPrefix "${workspaceMemberPackageDirs.${sourcePackage}}/" path && path != libraryTargetPath
                ) workspaceCargoTargetRelPaths
              )
            ) packageSourcePackages;
            testSourceFilesets = map (
              sourcePackage:
              let
                packageRoot = workspaceRoot + "/${workspaceMemberPackageDirs.${sourcePackage}}";
                testRoots = lib.filter builtins.pathExists (
                  map (root: packageRoot + "/${root}") [
                    "benches"
                    "examples"
                    "tests"
                  ]
                );
              in
              map (testRoot: lib.fileset.fileFilter (file: file.hasExt "rs") testRoot) testRoots
            ) testSourcePackages;
            runtimeTestSourceFilesets = map (
              sourcePackage:
              let
                packageRoot = workspaceRoot + "/${workspaceMemberPackageDirs.${sourcePackage}}";
              in
              lib.fileset.fileFilter (file: file.hasExt "snap") packageRoot
            ) runtimeTestSourcePackages;
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

        workspacePackageSrcFor =
          package:
          workspacePackageSrcForSourcePackages {
            sourcePackages = workspaceNormalSourcePackageNamesFor package;
            packageSourcePackages = [ package ];
          };

        workspaceTestPackageSrcFor =
          package:
          workspacePackageSrcForSourcePackages {
            sourcePackages = workspaceTestSourcePackageNamesFor package;
            packageSourcePackages = [ package ];
            testSourcePackages = [ package ];
          };

        workspacePackageIsProcMacro =
          package:
          let
            manifest = workspaceManifestFor workspaceMemberPackageDirs.${package};
            libManifest = manifest.lib or { };
            crateTypes = libManifest."crate-type" or [ ];
          in
          (libManifest."proc-macro" or false)
          || (libManifest.proc_macro or false)
          || builtins.elem "proc-macro" crateTypes
          || builtins.elem "proc_macro" crateTypes;

        workspacePackageHasLibTarget =
          package:
          let
            manifest = workspaceManifestFor workspaceMemberPackageDirs.${package};
            packageDir = workspaceMemberPackageDirs.${package};
          in
          manifest ? lib || builtins.pathExists (workspaceRoot + "/${packageDir}/src/lib.rs");

        workspacePackageLibTargetRelPath =
          package:
          let
            manifest = workspaceManifestFor workspaceMemberPackageDirs.${package};
          in
          "${workspaceMemberPackageDirs.${package}}/${manifest.lib.path or "src/lib.rs"}";

        workspacePackageForCargoTargetPath =
          path:
          lib.findFirst (
            package: lib.hasPrefix "${workspaceMemberPackageDirs.${package}}/" path
          ) null workspaceMemberPackages;

        workspaceDummyCargoTargetForPath =
          path:
          let
            package = workspacePackageForCargoTargetPath path;
          in
          if
            package != null
            && workspacePackageIsProcMacro package
            && path == workspacePackageLibTargetRelPath package
          then
            dummyProcMacroCargoTarget
          else
            dummyCargoTarget;

        workspaceMissingCargoTargetsScript = ''
          ${lib.concatMapStringsSep "\n" (path: ''
            if [ ! -e ${lib.escapeShellArg path} ]; then
              install -D -m 0644 ${workspaceDummyCargoTargetForPath path} ${lib.escapeShellArg path}
            fi
          '') workspaceCargoTargetRelPaths}

          ${normalizeWorkspaceHackBuildScriptTimestampScript}
        '';

        workspaceAllDummyCargoTargetsScript = ''
          ${lib.concatMapStringsSep "\n" (path: ''
            install -D -m 0644 ${workspaceDummyCargoTargetForPath path} "$out"/${lib.escapeShellArg path}
          '') workspaceCargoTargetRelPaths}

          ${normalizeWorkspaceHackBuildScriptTimestampInDummySrcScript}
        '';

        workspaceAllDummyCargoTargetsPreservingProcMacrosScript = ''
          ${lib.concatMapStringsSep "\n" (
            path:
            let
              package = workspacePackageForCargoTargetPath path;
              keepRealProcMacro =
                package != null
                && workspacePackageIsProcMacro package
                && path == workspacePackageLibTargetRelPath package;
              source = workspaceRoot + "/${path}";
            in
            if keepRealProcMacro then
              ''
                install -D -m 0644 ${source} "$out"/${lib.escapeShellArg path}
              ''
            else
              ''
                install -D -m 0644 ${workspaceDummyCargoTargetForPath path} "$out"/${lib.escapeShellArg path}
              ''
          ) workspaceCargoTargetRelPaths}

          ${normalizeWorkspaceHackBuildScriptTimestampInDummySrcScript}
        '';

        workspacePackageSourceRestoreInDummySrcScriptFor = sourcePackages: ''
          ${lib.concatMapStringsSep "\n" (
            sourcePackage:
            let
              packageDir = workspaceMemberPackageDirs.${sourcePackage};
              packageParentDir = builtins.dirOf packageDir;
              sourceRoot = workspacePackageSrcForSourcePackages {
                sourcePackages = [ sourcePackage ];
              };
              source = sourceRoot + "/${packageDir}";
            in
            ''
              rm -rf "$out"/${lib.escapeShellArg packageDir}
              mkdir -p "$out"/${lib.escapeShellArg packageParentDir}
              cp -R --no-preserve=ownership ${source} "$out"/${lib.escapeShellArg packageDir}
              chmod -R u+w "$out"/${lib.escapeShellArg packageDir}
            ''
          ) sourcePackages}

          ${workspacePackageExtraSourceRestoreInDummySrcScriptFor sourcePackages}
        '';

        workspaceDependencyDummyCargoTargetsScriptFor =
          package:
          let
            dependencySourcePackages = lib.filter (sourcePackage: sourcePackage != package) (
              workspaceNormalSourcePackageNamesFor package
            );
            packageTargetPaths = lib.filter (
              path: lib.hasPrefix "${workspaceMemberPackageDirs.${package}}/" path
            ) workspaceCargoTargetRelPaths;
          in
          ''
            ${workspacePackageSourceRestoreInDummySrcScriptFor dependencySourcePackages}

            ${lib.concatMapStringsSep "\n" (path: ''
              install -D -m 0644 ${workspaceDummyCargoTargetForPath path} "$out"/${lib.escapeShellArg path}
            '') packageTargetPaths}

            ${normalizeWorkspaceHackBuildScriptTimestampInDummySrcScript}
          '';

        packageCargoFeatures =
          package: (workspaceManifestFor workspaceMemberPackageDirs.${package}).features or { };

        packageFeatureIf =
          package: feature: lib.optional (builtins.hasAttr feature (packageCargoFeatures package)) feature;

        craneCiCommonFeaturesFor =
          package:
          lib.concatMap (packageFeatureIf package) [
            "bincode"
            "serde"
          ];

        craneCiExtraFeatureSets = {
          "gammaloop-workspace-hack" = [ "symbolica/tracing_max_level_info" ];
          "gammaloop-tracing-filter" = [
            "clap"
            "symbolica"
          ];
          idenso = [ "reference-cases" ];
          linnet = [
            "drawing"
            "symbolica"
          ];
          spenso = [ "shadowing" ];
          "spenso-macros" = [ "shadowing" ];
        };

        workspaceFeatureUnificationExcludedPackages = [
          "alphal00p-docs-python-exporter"
          "linnet-py"
          "spynso3"
        ];

        workspaceIncomingNormalDependencyFeaturesFor =
          dependency:
          sortedUnique (
            lib.concatMap (package: workspaceGraph.normal_dependency_features.${package}.${dependency} or [ ]) (
              lib.subtractLists workspaceFeatureUnificationExcludedPackages workspaceMemberPackages
            )
          );

        workspaceIncomingTestDependencyFeaturesFor =
          packages: dependency:
          sortedUnique (
            lib.concatMap (
              package: workspaceGraph.test_dependency_features.${package}.${dependency} or [ ]
            ) packages
          );

        craneCiFeaturesFor =
          package:
          sortedUnique (
            craneCiCommonFeaturesFor package
            ++ (craneCiExtraFeatureSets.${package} or [ ])
            ++ (workspaceIncomingNormalDependencyFeaturesFor package)
          );

        craneTestExtraFeatureSets = {
          "spenso-macros" = [ "spenso/shadowing" ];
        };

        craneTestFeaturesFor =
          package: sortedUnique (craneCiFeaturesFor package ++ (craneTestExtraFeatureSets.${package} or [ ]));

        craneTestContextFeaturesFor =
          sourcePackages: package:
          sortedUnique (
            craneCiCommonFeaturesFor package
            ++ (craneCiExtraFeatureSets.${package} or [ ])
            ++ (craneTestExtraFeatureSets.${package} or [ ])
            ++ (workspaceIncomingTestDependencyFeaturesFor sourcePackages package)
          );

        workspaceTestContextFor =
          {
            packages,
            extraFeatures ? { },
          }:
          let
            componentRepresentatives = sortedUnique (map workspaceTestComponentRepresentativeFor packages);
            componentPackages = sortedUnique (
              lib.concatMap (
                representative: workspaceTestComponentMembers.${representative}
              ) componentRepresentatives
            );
            sourcePackages = sortedUnique (lib.concatMap workspaceTestSourcePackageNamesFor componentPackages);
            anchorPackages = lib.optionals (builtins.any (
              package:
              package != workspaceHackPackage && builtins.elem package workspaceGraph.symbolica_test_packages
            ) sourcePackages) [ workspaceHackPackage ];
            featurePackages = sortedUnique (sourcePackages ++ anchorPackages);
            features = lib.listToAttrs (
              map (package: {
                name = package;
                value = sortedUnique (
                  craneTestContextFeaturesFor sourcePackages package ++ (extraFeatures.${package} or [ ])
                );
              }) featurePackages
            );
            resolvedFeatureVector = map (package: {
              inherit package;
              features = features.${package};
            }) featurePackages;
            usesPythonModule = builtins.any (
              package: builtins.elem "python-api-tests" features.${package}
            ) featurePackages;
            compileEnvironment = {
              inherit (ciArgs) NO_SYMBOLICA_OEM_LICENSE;
              inherit (commonArgs) CC CXX RUSTFLAGS;
              PYO3_PYTHON = if usesPythonModule then "${nextestPython}/bin/python3" else ciArgs.PYO3_PYTHON;
            };
            key = builtins.substring 0 16 (
              builtins.hashString "sha256" (
                builtins.toJSON {
                  inherit compileEnvironment resolvedFeatureVector;
                  profile = ciCargoProfile;
                  target = nextestTargetTriple;
                }
              )
            );
          in
          {
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

        workspaceTestDependencyContextsFor =
          context:
          let
            dependencyRepresentatives = sortedUnique (
              lib.concatMap workspaceTestComponentDependencyRepresentativesFor context.componentRepresentatives
            );
            dependencyContexts = map (
              representative:
              workspaceTestContextFor {
                packages = workspaceTestComponentMembers.${representative};
                inherit (context) extraFeatures;
              }
            ) dependencyRepresentatives;
          in
          builtins.attrValues (
            builtins.removeAttrs (lib.listToAttrs (
              map (dependencyContext: {
                name = dependencyContext.key;
                value = dependencyContext;
              }) dependencyContexts
            )) [ context.key ]
          );

        cargoFeatureArgs =
          features:
          lib.optionalString (
            features != [ ]
          ) "--features ${lib.escapeShellArg (lib.concatStringsSep "," features)}";

        cargoPackageArgsFor =
          package: features:
          lib.concatStringsSep " " (
            [
              "--locked"
              "-p ${lib.escapeShellArg package}"
            ]
            ++ lib.optional (features != [ ]) (cargoFeatureArgs features)
          );

        cargoQualifiedFeaturesFor =
          packages: featuresFor:
          let
            features = lib.concatMap (
              package:
              map (feature: if lib.hasInfix "/" feature then feature else "${package}/${feature}") (
                featuresFor package
              )
            ) packages;
          in
          sortedUnique features;

        cargoQualifiedFeatureArgsFor =
          packages: featuresFor: cargoFeatureArgs (cargoQualifiedFeaturesFor packages featuresFor);

        cargoPackageArgsWithFeaturePackagesFor =
          {
            package,
            featurePackages,
            featuresFor,
            selectFeaturePackages ? true,
          }:
          let
            anchorPackages = lib.optionals (
              package != workspaceHackPackage && builtins.elem package workspaceGraph.symbolica_normal_packages
            ) [ workspaceHackPackage ];
            selectedFeaturePackages =
              if selectFeaturePackages then
                lib.filter (
                  featurePackage: (featuresFor featurePackage) != [ ] && !workspacePackageIsProcMacro featurePackage
                ) featurePackages
              else
                [ ];
            selectedPackages = sortedUnique ([ package ] ++ anchorPackages ++ selectedFeaturePackages);
            featureArgs = cargoQualifiedFeatureArgsFor (sortedUnique (
              featurePackages ++ anchorPackages
            )) featuresFor;
          in
          lib.concatStringsSep " " (
            [
              "--locked"
            ]
            ++ map (selectedPackage: "-p ${lib.escapeShellArg selectedPackage}") selectedPackages
            ++ lib.optional (featureArgs != "") featureArgs
          );

        cargoPackagesArgsFor =
          packages: featuresFor:
          let
            featureArgs = cargoQualifiedFeatureArgsFor packages featuresFor;
          in
          lib.concatStringsSep " " (
            [
              "--locked"
            ]
            ++ map (package: "-p ${lib.escapeShellArg package}") packages
            ++ lib.optional (featureArgs != "") featureArgs
          );

        craneWorkspacePrebuildFeatureArgs = cargoQualifiedFeatureArgsFor workspaceMemberPackages craneTestFeaturesFor;
        workspacePrebuildExcludedPackages = workspaceFeatureUnificationExcludedPackages;
        workspacePrebuildDependencyPackages = lib.filter (
          package:
          package != workspaceHackPackage
          && !(builtins.elem package workspacePrebuildExcludedPackages)
          && workspacePackageHasLibTarget package
        ) workspaceMemberPackages;
        workspacePrebuildDependencyPathFor =
          package:
          let
            packageDir = workspaceMemberPackageDirs.${package};
          in
          if lib.hasPrefix "crates/" packageDir then
            "../${lib.removePrefix "crates/" packageDir}"
          else
            "../../${packageDir}";
        workspacePrebuildCargoToml = pkgs.writeText "${workspacePrebuildPackage}-Cargo.toml" ''
          [package]
          name = "${workspacePrebuildPackage}"
          version = "0.1.0"
          edition = "2024"
          publish = false

          [lib]
          path = "src/lib.rs"

          [dependencies]
          ${lib.concatMapStringsSep "\n" (
            package:
            let
              features = lib.filter (feature: !lib.hasInfix "/" feature) (craneTestFeaturesFor package);
              featureEntry = lib.optionalString (features != [ ]) ", features = ${builtins.toJSON features}";
            in
            ''
              ${package} = { path = "${workspacePrebuildDependencyPathFor package}"${featureEntry} }
            ''
          ) workspacePrebuildDependencyPackages}
        '';
        workspacePrebuildSourceScript = ''
          install -D -m 0644 ${workspacePrebuildCargoToml} "$out/${workspacePrebuildPackageDir}/Cargo.toml"
          install -D -m 0644 ${dummyCargoTarget} "$out/${workspacePrebuildPackageDir}/src/lib.rs"
        '';
        workspacePrebuildCargoArgs = lib.concatStringsSep " " (
          [
            "--offline"
            "-p ${lib.escapeShellArg workspacePrebuildPackage}"
            "-p ${lib.escapeShellArg workspaceHackPackage}"
          ]
          ++ lib.optional (craneWorkspacePrebuildFeatureArgs != "") craneWorkspacePrebuildFeatureArgs
        );
        workspaceConsumerPackageFor = package: "gammaloop-ci-consumer-${package}";
        workspaceConsumerPackageDirFor = package: "crates/${workspaceConsumerPackageFor package}";
        workspaceConsumerCargoTomlFor =
          package: dependencyPackages:
          pkgs.writeText "${workspaceConsumerPackageFor package}-Cargo.toml" ''
            [package]
            name = "${workspaceConsumerPackageFor package}"
            version = "0.1.0"
            edition = "2024"
            publish = false

            [lib]
            path = "src/lib.rs"

            [dependencies]
            ${lib.concatMapStringsSep "\n" (
              dependencyPackage:
              let
                features = lib.filter (feature: !lib.hasInfix "/" feature) (craneCiFeaturesFor dependencyPackage);
                featureEntry = lib.optionalString (features != [ ]) ", features = ${builtins.toJSON features}";
              in
              ''
                ${dependencyPackage} = { path = "${workspacePrebuildDependencyPathFor dependencyPackage}"${featureEntry} }
              ''
            ) dependencyPackages}
          '';
        workspaceConsumerSourceScriptFor = package: dependencyPackages: ''
          install -D -m 0644 ${workspaceConsumerCargoTomlFor package dependencyPackages} "$out/${workspaceConsumerPackageDirFor package}/Cargo.toml"
          install -D -m 0644 ${dummyCargoTarget} "$out/${workspaceConsumerPackageDirFor package}/src/lib.rs"
        '';
        cargoPackageDependencyModeArgsFor =
          package:
          let
            sourcePackages = workspaceNormalSourcePackageNamesFor package;
            anchorPackages = lib.optionals (
              package != workspaceHackPackage && builtins.elem package workspaceGraph.symbolica_normal_packages
            ) [ workspaceHackPackage ];
            crossFeatures = sortedUnique (
              lib.concatMap (
                featurePackage: lib.filter (feature: lib.hasInfix "/" feature) (craneCiFeaturesFor featurePackage)
              ) (sortedUnique (sourcePackages ++ anchorPackages))
            );
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
        cargoPackageCiArgsFor =
          package:
          cargoPackageArgsWithFeaturePackagesFor {
            inherit package;
            featurePackages = workspaceNormalSourcePackageNamesFor package;
            featuresFor = craneCiFeaturesFor;
          };
        testDependencyFeatureAnchorPackageFor = context: "gammaloop-ci-test-dependencies-${context.key}";
        testDependencyFeatureAnchorPackageDirFor =
          context: "crates/${testDependencyFeatureAnchorPackageFor context}";
        testDependencyFeatureAnchorDependencyPackagesFor =
          context: lib.filter workspacePackageHasLibTarget context.sourcePackages;
        testDependencyFeatureAnchorCargoTomlFor =
          context:
          pkgs.writeText "${testDependencyFeatureAnchorPackageFor context}-Cargo.toml" ''
            [package]
            name = "${testDependencyFeatureAnchorPackageFor context}"
            version = "0.1.0"
            edition = "2024"
            publish = false

            [lib]
            path = "src/lib.rs"

            [dependencies]
            ${lib.concatMapStringsSep "\n" (
              package:
              let
                features = lib.filter (feature: !lib.hasInfix "/" feature) context.features.${package};
                featureEntry = lib.optionalString (features != [ ]) ", features = ${builtins.toJSON features}";
              in
              ''
                ${package} = { path = "${workspacePrebuildDependencyPathFor package}"${featureEntry} }
              ''
            ) (testDependencyFeatureAnchorDependencyPackagesFor context)}
          '';
        testDependencyFeatureAnchorSourceScriptFor = context: prefix: ''
          install -D -m 0644 ${testDependencyFeatureAnchorCargoTomlFor context} "${prefix}${testDependencyFeatureAnchorPackageDirFor context}/Cargo.toml"
          install -D -m 0644 ${dummyCargoTarget} "${prefix}${testDependencyFeatureAnchorPackageDirFor context}/src/lib.rs"
          touch -d @0 "${prefix}${testDependencyFeatureAnchorPackageDirFor context}/Cargo.toml" "${prefix}${testDependencyFeatureAnchorPackageDirFor context}/src/lib.rs"
        '';
        testBinaryFeatureAnchorPackageFor = context: "gammaloop-ci-test-binary-dependencies-${context.key}";
        testBinaryFeatureAnchorPackageDirFor =
          context: "crates/${testBinaryFeatureAnchorPackageFor context}";
        testBinaryFeatureAnchorDevDependenciesFor =
          context:
          lib.foldl' lib.recursiveUpdate { } (
            map (
              package:
              builtins.removeAttrs ((workspaceManifestFor workspaceMemberPackageDirs.${package})
                ."dev-dependencies" or { }
              ) context.componentPackages
            ) context.componentPackages
          );
        testBinaryFeatureAnchorCrossFeaturesFor =
          context:
          lib.foldl'
            (
              features: feature:
              let
                parts = lib.splitString "/" feature;
                dependency = builtins.head parts;
              in
              features
              // {
                ${dependency} = sortedUnique ((features.${dependency} or [ ]) ++ [ (builtins.elemAt parts 1) ]);
              }
            )
            { }
            (
              lib.concatMap (
                package: lib.filter (feature: lib.hasInfix "/" feature) context.features.${package}
              ) context.featurePackages
            );
        testBinaryFeatureAnchorDependenciesFor =
          context:
          let
            inheritedDependencies = lib.listToAttrs (
              map
                (package: {
                  name = package;
                  value.path = workspacePrebuildDependencyPathFor package;
                })
                (
                  lib.subtractLists context.componentPackages (
                    testDependencyFeatureAnchorDependencyPackagesFor context
                  )
                )
            );
            rawDependencies =
              inheritedDependencies
              // (testBinaryFeatureAnchorDevDependenciesFor context)
              // lib.optionalAttrs (context.anchorPackages != [ ]) {
                ${workspaceHackPackage} = {
                  path = "../${workspaceHackPackage}";
                };
              };
            crossFeatures = testBinaryFeatureAnchorCrossFeaturesFor context;
            dependencyNames = sortedUnique (
              builtins.attrNames rawDependencies ++ builtins.attrNames crossFeatures
            );
          in
          lib.listToAttrs (
            map (
              dependency:
              let
                rawDependency = rawDependencies.${dependency} or { workspace = true; };
                dependencyAttrs =
                  if builtins.isAttrs rawDependency then rawDependency else { version = rawDependency; };
                features = sortedUnique (
                  (dependencyAttrs.features or [ ])
                  ++ lib.filter (feature: !lib.hasInfix "/" feature) (context.features.${dependency} or [ ])
                  ++ (crossFeatures.${dependency} or [ ])
                );
              in
              {
                name = dependency;
                value = dependencyAttrs // lib.optionalAttrs (features != [ ]) { inherit features; };
              }
            ) dependencyNames
          );
        testBinaryFeatureAnchorCargoTomlFor =
          context:
          (pkgs.formats.toml { }).generate "${testBinaryFeatureAnchorPackageFor context}-Cargo.toml" {
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
        testBinaryFeatureAnchorDependencyPathFor =
          context: package:
          if lib.hasPrefix "crates/" workspaceMemberPackageDirs.${package} then
            "../${testBinaryFeatureAnchorPackageFor context}"
          else
            "../${testBinaryFeatureAnchorPackageDirFor context}";
        testBinaryFeatureAnchorDevDependencyFor =
          context: package:
          pkgs.writeText "${testBinaryFeatureAnchorPackageFor context}-${package}-dev-dependency.toml" ''

            [dev-dependencies.${testBinaryFeatureAnchorPackageFor context}]
            path = "${testBinaryFeatureAnchorDependencyPathFor context package}"
          '';
        testBinaryFeatureAnchorDevDependencyScriptFor = context: package: prefix: ''
          cat ${testBinaryFeatureAnchorDevDependencyFor context package} >> "${prefix}${
            workspaceMemberPackageDirs.${package}
          }/Cargo.toml"
          touch -d @0 "${prefix}${workspaceMemberPackageDirs.${package}}/Cargo.toml"
        '';
        cargoTestDependencyArgsFor =
          context:
          let
            selectedPackages = sortedUnique (
              context.anchorPackages
              ++ [
                (testDependencyFeatureAnchorPackageFor context)
              ]
            );
            featureArgs = cargoQualifiedFeatureArgsFor context.featurePackages (
              package: context.features.${package}
            );
          in
          lib.concatStringsSep " " (
            [
              "--offline"
            ]
            ++ map (selectedPackage: "-p ${lib.escapeShellArg selectedPackage}") selectedPackages
            ++ lib.optional (featureArgs != "") featureArgs
          );
        cargoTestContextArgsFor =
          context: packages:
          let
            featureArgs = cargoQualifiedFeatureArgsFor context.featurePackages (
              featurePackage: context.features.${featurePackage}
            );
          in
          lib.concatStringsSep " " (
            [ "--offline" ]
            ++ map (package: "-p ${lib.escapeShellArg package}") (sortedUnique packages)
            ++ lib.optional (featureArgs != "") featureArgs
          );
        cargoTestBinaryDependencyArgsFor =
          context:
          cargoTestContextArgsFor context (
            [ (testBinaryFeatureAnchorPackageFor context) ] ++ context.componentPackages
          );
        cargoTestBinaryArgsFor =
          context: package:
          let
            featureArgs = cargoQualifiedFeatureArgsFor [ package ] (
              featurePackage: context.features.${featurePackage}
            );
          in
          lib.concatStringsSep " " (
            [
              "--offline"
              "-p ${lib.escapeShellArg package}"
            ]
            ++ lib.optional (featureArgs != "") featureArgs
          );
        cranePythonExtraFeatureSets = {
          "gammaloop-api" = [
            "python_abi"
            "pyo3-extension-module"
          ];
        };
        cranePythonFeaturesFor =
          package:
          sortedUnique (craneCiFeaturesFor package ++ (cranePythonExtraFeatureSets.${package} or [ ]));
        cranePythonCargoArgs =
          let
            featurePackages = workspaceNormalSourcePackageNamesFor "gammaloop-api";
            selectedFeaturePackages = lib.filter (
              featurePackage:
              (cranePythonFeaturesFor featurePackage) != [ ] && !workspacePackageIsProcMacro featurePackage
            ) featurePackages;
            selectedPackages = sortedUnique (
              [
                "gammaloop-api"
                workspaceHackPackage
              ]
              ++ selectedFeaturePackages
            );
            featureArgs = cargoQualifiedFeatureArgsFor featurePackages cranePythonFeaturesFor;
          in
          lib.concatStringsSep " " (
            [
              "--locked"
              "--no-default-features"
              "--lib"
            ]
            ++ map (package: "-p ${lib.escapeShellArg package}") selectedPackages
            ++ lib.optional (featureArgs != "") featureArgs
          );

        guppyFeatureMapFor =
          featuresFor:
          builtins.toJSON (
            lib.listToAttrs (
              map (package: {
                name = package;
                value = featuresFor package;
              }) workspaceMemberPackages
            )
          );

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
        cargoLinkerVar = "CARGO_TARGET_${
          lib.toUpper (lib.replaceStrings [ "-" ] [ "_" ] rustTarget)
        }_LINKER";

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

          nativeBuildInputs = [
            pkgs.pkg-config
            pkgs.gcc
            pkgs.python313
            pkgs.gnum4
          ]
          ++ lib.optionals pkgs.stdenv.isDarwin [
            pkgs.darwin.cctools
          ];

          buildInputs = [
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
        docsCargoProfile = "docs";

        ciArgs = commonArgs // {
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

        gammaloop-cli =
          pkgs.runCommand "gammaloop-api-${apiMeta.version}"
            {
              nativeBuildInputs = [
                pkgs.coreutils
                pkgs.findutils
              ];
            }
            ''
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
        gammaloop-python-lib = craneLib.buildPackage (
          ciArgs
          // {
            cargoArtifacts = cranePythonBuildArtifacts;
            pname = "gammaloop-api-python";
            src = workspacePackageSrcFor "gammaloop-api";
            cargoExtraArgs = cranePythonCargoArgs;
            doCheck = false;
            postPatch = workspaceMissingCargoTargetsScript;
          }
        );
        gammaloop-python-lib-output = lib.getLib gammaloop-python-lib;
        pythonSitePackages = "${pkgs.python313.sitePackages}";
        gammaloop-python-module = pkgs.runCommand "gammaloop-python-module" { } ''
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
        clinnetArgs = ciArgs // {
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

        clinnetCargoArtifacts = craneLib.buildDepsOnly (
          clinnetArgs
          // {
            preBuild = drawingTypstBundleAssets;
          }
        );

        clinnet-cli = craneLib.buildPackage (
          clinnetArgs
          // {
            cargoArtifacts = clinnetCargoArtifacts;
            doNotLinkInheritedArtifacts = true;
            preBuild = drawingTypstBundleAssets;
          }
        );

        rscls = pkgs.rustPlatform.buildRustPackage rec {
          pname = "rscls";
          version = "0.2.3";
          src = pkgs.fetchCrate {
            inherit pname version;
            sha256 = "sha256-tahAhWCjhIVjbJ1NzrtiHBwGb/FBmUdK4XP9VlSPqh0=";
          };
          cargoHash = "sha256-JikjBTFeDh4XHBm57yiorsCwZhKikz0aiWNOTaMn0Vo=";
        };

        devShellPackages =
          with pkgs;
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
            docsTypst
            roboto
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

        mkDevShell =
          extraPackages:
          craneLib.devShell {
            # checks = self.checks.${system};

            RUST_SRC_PATH = "${pkgs.rustPlatform.rustLibSrc}";
            GLIBC_TUNABLES = "glibc.rtld.optional_static_tls=10000";
            TYPST_FONT_PATHS = docsFontPath;

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

        nextestProfile = "ci_gammaloop";
        nextestJunitPath = "target/nextest/${nextestProfile}/junit.xml";

        nextestFailureSummary = pkgs.writeShellApplication {
          name = "nextest-failure-summary";
          runtimeInputs = [ pkgs.python313 ];
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
        # Keep this input to manifests so source-only commits can reuse the
        # dependency artifact from the NixCI cache.
        cargoArtifacts = buildDepsOnlyWithArtifacts (
          ciArgs
          // {
            pname = "gammaloop-workspace-deps";
            src = workspaceDependencySrc;
            cargoArtifacts = workspaceHackBuildArtifacts;
            cargoExtraArgs = workspacePrebuildCargoArgs;
            buildPhaseCargoCommand = "cargoWithProfile build ${workspacePrebuildCargoArgs}";
            checkPhaseCargoCommand = "";
            doCheck = false;
            stripWorkspaceArtifacts = true;
            extraDummyScript = ''
              ${workspacePrebuildSourceScript}

              ${lib.concatMapStringsSep "\n" (path: ''
                install -D -m 0644 ${dummyCargoTarget} "$out/${path}"
              '') autoCargoTargetPaths}
            '';
          }
        );

        documentationRegistry = builtins.fromTOML (builtins.readFile ./docs/products/registry.toml);
        alphal00pDocsCargoTargetRoot = "target/alphal00p-docs-rustdoc";
        alphal00pDocsCargoTarget = "${alphal00pDocsCargoTargetRoot}/cargo-target-v1";
        alphal00pDocsCargoArgs = commonArgs // {
          buildType = docsCargoProfile;
          CARGO_PROFILE = docsCargoProfile;
          CARGO_TARGET_DIR = alphal00pDocsCargoTarget;
          PYO3_PYTHON = "${pkgs.python313}/bin/python3";
          PYTHONPATH = "${pkgs.python313}/lib/python3.13/site-packages";
          # Keep the compile-time Symbolica setting explicit and identical in
          # both the reusable producer and the documentation consumer.
          SYMBOLICA_OEM_LICENSE =
            (builtins.fromTOML (builtins.readFile ./.cargo/config.toml)).env.SYMBOLICA_OEM_LICENSE.value;
        };
        alphal00pDocsRealCargoArgs = alphal00pDocsCargoArgs // {
          postPatch = normalizeWorkspaceHackBuildScriptTimestampScript;
        };
        documentationRustComponents = lib.concatMap (
          product: product.rust_components or [ ]
        ) documentationRegistry.product;
        documentationCatalogFeatures = "gammaloop-reference,vakint-reference";
        documentationPythonExporterBuildCommands = lib.concatMapStringsSep "\n" (product: ''
          cargoWithProfile build --locked -p alphal00p-docs-python-exporter --features ${lib.escapeShellArg product.id}
        '') documentationRegistry.product;
        documentationRustdocBuildCommands = lib.concatMapStringsSep "\n" (
          component:
          let
            features = component.features or [ ];
            featureArgs = lib.optionalString (features != [ ]) (
              " --features ${lib.escapeShellArg (lib.concatStringsSep "," features)}"
            );
          in
          ''
            cargoWithProfile doc --locked --no-deps --no-default-features -p ${lib.escapeShellArg component.package}${featureArgs}
          ''
        ) documentationRustComponents;

        # Keep real workspace outputs keyed only by Rust and build inputs so
        # documentation-only edits reuse the complete Cargo build. End on the
        # first consumer context after the isolated exporter/Rustdoc matrix.
        alphal00pDocsCargoArtifacts = craneLib.mkCargoDerivation (
          alphal00pDocsRealCargoArgs
          // {
            cargoArtifacts = null;
            doInstallCargoArtifacts = true;
            pname = "alphal00p-docs-cargo";
            src = documentationWorkspaceBuildSrc;
            # Prime content-sensitive test dependency contexts without adding
            # authored documentation to this reusable source boundary.
            postPatch = normalizeWorkspaceHackBuildScriptTimestampScript + ''
              install -D -m 0644 ${dummyCargoTarget} crates/alphal00p-docs-examples/build.rs
              install -D -m 0644 ${dummyCargoTarget} crates/alphal00p-docs-examples/src/lib.rs
            '';
            buildPhaseCargoCommand = ''
              cargoWithProfile build --locked -p alphal00p-docs-catalogs --bin alphal00p-docs-catalogs
              ${documentationPythonExporterBuildCommands}
              cargoWithProfile test --locked --no-run -p alphal00p-docs-examples
              cargo clean --profile ${docsCargoProfile} -p alphal00p-docs-examples
              cargoWithProfile build --locked -p alphal00p-docs-builder
              cargoWithProfile build --locked -p linnet-py --lib --features extension-module,abi3-py310
              ${documentationRustdocBuildCommands}
              cargoWithProfile build --locked -p alphal00p-docs-catalogs \
                --features ${lib.escapeShellArg documentationCatalogFeatures} \
                --bin alphal00p-docs-gammaloop-reference \
                --bin alphal00p-docs-vakint-reference
            '';
            checkPhaseCargoCommand = "";
            doCheck = false;
            installPhaseCommand = "";
          }
        );

        guppyWorkspaceGraphNativeBuildInputs = [
          ciToolchain
          pkgs.cargo-guppy
          pkgs.coreutils
          pkgs.gawk
          pkgs.gnugrep
          pkgs.jq
        ];

        guppyWorkspaceGraphJson =
          pkgs.runCommand "gammaloop-guppy-workspace-graph.json"
            {
              nativeBuildInputs = guppyWorkspaceGraphNativeBuildInputs;
            }
            ''
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
                  if guppy_names "$package" all "$include_dev" | grep -E '^(symbolica|numerica|graphica)$' >/dev/null; then
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

        guppyWorkspaceGraphCheck =
          pkgs.runCommand "gammaloop-guppy-workspace-graph-check"
            {
              nativeBuildInputs = [
                pkgs.diffutils
              ];
            }
            ''
              diff -u ${./nix/ci-workspace-graph.json} ${guppyWorkspaceGraphJson}
              mkdir -p "$out"
            '';

        mergeCargoArtifacts =
          name: artifacts:
          let
            artifactList = lib.filter (artifact: artifact != null) artifacts;
          in
          pkgs.runCommand name
            {
              nativeBuildInputs = [
                pkgs.rsync
                pkgs.zstd
                pkgs.gnutar
              ];
            }
            ''
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

              # Crane deletes inherited Cargo locks before using an artifact tree.
              # Remove them here while the merged files are still ordinary files.
              find "$out/target" -name .cargo-lock -delete
            '';

        mergeCargoArtifactsOrNull =
          name: artifacts:
          let
            artifactList = lib.filter (artifact: artifact != null) artifacts;
          in
          if artifactList == [ ] then null else mergeCargoArtifacts name artifactList;

        workspacePackageCargoArtifactNames =
          package:
          let
            manifest = workspaceManifestFor workspaceMemberPackageDirs.${package};
            libManifest = manifest.lib or { };
            targetNameForCargo = name: lib.replaceStrings [ "-" ] [ "_" ] name;
          in
          [ (targetNameForCargo (libManifest.name or manifest.package.name)) ];

        stripSelectedWorkspaceCargoArtifactsScript =
          {
            cargoTargetDir ? "target",
            packagesToStrip,
          }:
          let
            targetDir = if cargoTargetDir == "target" then "target" else lib.escapeShellArg cargoTargetDir;
            stripPackageScript =
              package:
              let
                artifactNames = workspacePackageCargoArtifactNames package;
                fingerprintNames = sortedUnique ([ package ] ++ artifactNames);
              in
              ''
                ${lib.concatMapStringsSep "\n" (name: ''
                  rm -f ${targetDir}/${ciCargoProfile}/deps/${lib.escapeShellArg name}
                  rm -f ${targetDir}/${ciCargoProfile}/deps/${lib.escapeShellArg name}.d
                  rm -f ${targetDir}/${ciCargoProfile}/deps/${lib.escapeShellArg name}-*
                  rm -f ${targetDir}/${ciCargoProfile}/deps/lib${lib.escapeShellArg name}.*
                  rm -f ${targetDir}/${ciCargoProfile}/deps/lib${lib.escapeShellArg name}-*
                  rm -f ${targetDir}/${ciCargoProfile}/${lib.escapeShellArg name}
                  rm -f ${targetDir}/${ciCargoProfile}/${lib.escapeShellArg name}.d
                  rm -f ${targetDir}/${ciCargoProfile}/lib${lib.escapeShellArg name}.*
                  rm -f ${targetDir}/${ciCargoProfile}/lib${lib.escapeShellArg name}-*
                '') artifactNames}
                ${lib.concatMapStringsSep "\n" (name: ''
                  rm -rf ${targetDir}/${ciCargoProfile}/.fingerprint/${lib.escapeShellArg name}-[0-9a-f]*
                  rm -rf ${targetDir}/${ciCargoProfile}/build/${lib.escapeShellArg name}-[0-9a-f]*
                '') fingerprintNames}
              '';
          in
          ''
            # Cargo feature anchors can compile dummy workspace units. Keep only
            # artifacts that are intended to be reused by later derivations.
            ${lib.concatMapStringsSep "\n" stripPackageScript packagesToStrip}
          '';

        stripWorkspaceCargoArtifactsScript =
          {
            cargoTargetDir ? "target",
            preserveProcMacros ? false,
            preservedPackages ? [ ],
          }:
          let
            packagesToStrip = lib.filter (
              package:
              !(builtins.elem package preservedPackages)
              && !(preserveProcMacros && workspacePackageIsProcMacro package)
            ) workspacePackages;
          in
          stripSelectedWorkspaceCargoArtifactsScript {
            inherit cargoTargetDir packagesToStrip;
          };

        # Crane's buildDepsOnly deliberately clears cargoArtifacts. This local
        # variant keeps the same dummy-source behavior while allowing per-crate
        # dependency archives to build incrementally from earlier archives.
        buildDepsOnlyWithArtifacts =
          {
            cargoBuildCommand ? "cargoWithProfile build",
            cargoCheckCommand ? "cargoWithProfile check",
            cargoExtraArgs ? "--locked",
            cargoTestCommand ? "cargoWithProfile test",
            cargoTestExtraArgs ? "--no-run",
            preBuildWorkspaceArtifactStripPackages ? [ ],
            preservedWorkspaceArtifactPackages ? [ ],
            ...
          }@args:
          let
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
              if args ? dummySrc then
                lib.warnIf (
                  args ? src && args.src != null
                ) "buildDepsOnlyWithArtifacts will ignore `src` when `dummySrc` is specified" args.dummySrc
              else
                craneLib.mkDummySrc args;
            argsMaybeDummySrcOverride = if args ? dummySrc then args // { src = args.dummySrc; } else args;
            crateName = craneLib.crateNameFromCargoToml argsMaybeDummySrcOverride;
            stripWorkspaceArtifactsScriptText =
              lib.optionalString (args.stripWorkspaceArtifacts or false)
                (stripWorkspaceCargoArtifactsScript {
                  cargoTargetDir = args.CARGO_TARGET_DIR or "target";
                  preserveProcMacros = args.preserveWorkspaceProcMacroArtifacts or false;
                  preservedPackages = preservedWorkspaceArtifactPackages;
                });
            compactStripWorkspaceArtifactsScriptText =
              lib.optionalString (args.stripWorkspaceArtifacts or false)
                (stripWorkspaceCargoArtifactsScript {
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
                (cleanedArgs.nativeBuildInputs or [ ])
                ++ lib.optionals (args.stripWorkspaceArtifacts or false) [
                  pkgs.gnutar
                  pkgs.rsync
                  pkgs.zstd
                ];

              cargoArtifacts = args.cargoArtifacts or null;
              doNotLinkInheritedArtifacts = args.doNotLinkInheritedArtifacts or true;
              cargoVendorDir = args.cargoVendorDir or (craneLib.vendorCargoDeps argsMaybeDummySrcOverride);

              env = (args.env or { }) // {
                CRANE_BUILD_DEPS_ONLY = ((args.env or { }).CRANE_BUILD_DEPS_ONLY or 1);
              };

              postPatch = (args.postPatch or "") + ''

                ${normalizeWorkspaceHackBuildScriptTimestampScript}
              '';

              preBuild =
                (args.preBuild or "")
                +
                  lib.optionalString (preBuildWorkspaceArtifactStripPackages != [ ])
                    (stripSelectedWorkspaceCargoArtifactsScript {
                      cargoTargetDir = args.CARGO_TARGET_DIR or "target";
                      packagesToStrip = preBuildWorkspaceArtifactStripPackages;
                    });

              buildPhaseCargoCommand =
                args.buildPhaseCargoCommand or ''
                  ${cargoCheckCommand} ${cargoExtraArgs} ${cargoCheckExtraArgs}
                  ${cargoBuildCommand} ${cargoExtraArgs}
                '';

              checkPhaseCargoCommand =
                args.checkPhaseCargoCommand or ''
                  ${cargoTestCommand} ${cargoExtraArgs} ${cargoTestExtraArgs}
                '';

              postBuild = (args.postBuild or "") + stripWorkspaceArtifactsScriptText;

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
                      ${compactStripWorkspaceArtifactsScriptText}
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

        rootPackageDependencyArtifactsFor =
          package:
          lib.optional (
            package != workspaceHackPackage
            && (workspaceDependencyNamesFor package) == [ ]
            && builtins.elem package workspaceGraph.symbolica_normal_packages
          ) workspaceHackBuildArtifacts;

        workspaceHackDependencyArtifacts = buildDepsOnlyWithArtifacts (
          ciArgs
          // {
            pname = "gammaloop-crate-${workspaceHackPackage}";
            src = workspacePackageSrcFor workspaceHackPackage;
            cargoExtraArgs = cargoPackageCiArgsFor workspaceHackPackage;
            buildPhaseCargoCommand = "cargoWithProfile build ${cargoPackageCiArgsFor workspaceHackPackage}";
            checkPhaseCargoCommand = "";
            doCheck = false;
            preBuildWorkspaceArtifactStripPackages = [ workspaceHackPackage ];
            stripWorkspaceArtifacts = true;
            extraDummyScript = workspaceAllDummyCargoTargetsScript;
          }
        );

        workspaceHackBuildArtifacts = craneLib.cargoBuild (
          ciArgs
          // {
            cargoArtifacts = workspaceHackDependencyArtifacts;
            doNotLinkInheritedArtifacts = true;
            pname = "gammaloop-crate-${workspaceHackPackage}";
            src = workspacePackageSrcFor workspaceHackPackage;
            cargoExtraArgs = cargoPackageCiArgsFor workspaceHackPackage;
            postPatch = workspaceMissingCargoTargetsScript;
          }
        );

        cranePackageDependencyModeDependencyArtifacts = lib.fix (
          self:
          lib.genAttrs workspacePackages (
            package:
            if package == workspaceHackPackage then
              workspaceHackDependencyArtifacts
            else if !workspacePackageHasLibTarget package then
              null
            else
              let
                dependencyArtifactFor =
                  dependency:
                  if dependency == workspaceHackPackage then
                    workspaceHackBuildArtifacts
                  else
                    cranePackageDependencyModeArtifacts.${dependency};
                dependencySourcePackages = lib.filter (sourcePackage: sourcePackage != package) (
                  workspaceNormalSourcePackageNamesFor package
                );
                preservedWorkspaceArtifactPackages = lib.filter workspacePackageIsProcMacro dependencySourcePackages;
              in
              buildDepsOnlyWithArtifacts (
                (builtins.removeAttrs ciArgs [ "src" ])
                // {
                  cargoArtifacts = mergeCargoArtifactsOrNull "gammaloop-crate-${package}-dependency-deps-inputs" (
                    [ cargoArtifacts ]
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
                }
              )
          )
        );

        cranePackageDependencyModeArtifacts = lib.fix (
          self:
          lib.genAttrs workspacePackages (
            package:
            if package == workspaceHackPackage then
              workspaceHackBuildArtifacts
            else if !workspacePackageHasLibTarget package then
              null
            else
              let
                dependencyArtifactFor =
                  dependency:
                  if dependency == workspaceHackPackage then workspaceHackBuildArtifacts else self.${dependency};
                sourcePackages = workspaceNormalSourcePackageNamesFor package;
                preservedWorkspaceArtifactPackages = sourcePackages;
              in
              buildDepsOnlyWithArtifacts (
                (builtins.removeAttrs ciArgs [ "src" ])
                // {
                  cargoArtifacts = mergeCargoArtifactsOrNull "gammaloop-crate-${package}-dependency-inputs" (
                    [ cranePackageDependencyModeDependencyArtifacts.${package} ]
                    ++ map dependencyArtifactFor (lib.filter (sourcePackage: sourcePackage != package) sourcePackages)
                  );
                  pname = "gammaloop-crate-${package}-dependency";
                  src = workspacePackageSrcForSourcePackages {
                    inherit sourcePackages;
                  };
                  buildPhaseCargoCommand = "cargoWithProfile build ${cargoPackageDependencyModeArgsFor package}";
                  checkPhaseCargoCommand = "";
                  doCheck = false;
                  preBuildWorkspaceArtifactStripPackages = [ package ];
                  stripWorkspaceArtifacts = true;
                  inherit preservedWorkspaceArtifactPackages;
                  extraDummyScript = ''
                    ${workspacePackageSourceRestoreInDummySrcScriptFor sourcePackages}
                    ${workspaceConsumerSourceScriptFor package sourcePackages}
                  '';
                  postPatch = workspaceMissingCargoTargetsScript;
                }
              )
          )
        );

        # Public crate-deps outputs are the reusable workspace crate artifacts.
        # Build them through the consumer-anchor path so each package's real
        # artifact is preserved for downstream package and test derivations.
        cranePackageDependencyArtifacts = cranePackageDependencyModeArtifacts;

        cranePackageBuildArtifacts = lib.fix (
          self:
          lib.genAttrs workspacePackages (
            package:
            if package == workspaceHackPackage then
              workspaceHackBuildArtifacts
            else
              craneLib.cargoBuild (
                ciArgs
                // {
                  cargoArtifacts =
                    if cranePackageDependencyArtifacts.${package} == null then
                      cargoArtifacts
                    else
                      cranePackageDependencyArtifacts.${package};
                  pname = "gammaloop-crate-${package}";
                  src = workspacePackageSrcFor package;
                  cargoExtraArgs =
                    cargoPackageCiArgsFor package + lib.optionalString (package == "gammaloop-api") " --lib --bins";
                  preBuild = lib.optionalString (package == "gammaloop-api") drawingTypstBundleAssets;
                  postPatch = workspaceMissingCargoTargetsScript;
                }
              )
          )
        );

        gammaloopApiPackageArtifacts = mergeCargoArtifacts "gammaloop-api-package-artifacts" [
          cranePackageBuildArtifacts."gammaloop-api"
        ];

        workspaceBaseTestContexts = map (
          package: workspaceTestContextFor { packages = [ package ]; }
        ) workspacePackages;
        nextestPackageTestContexts = lib.concatMap (
          target:
          map (
            package:
            workspaceTestContextFor {
              packages = [ package ];
              extraFeatures = target.extraFeatures or { };
            }
          ) target.packages
        ) checkedNextestPackageGroups;
        workspaceTestContexts = lib.listToAttrs (
          map (context: {
            name = context.key;
            value = context;
          }) (workspaceBaseTestContexts ++ nextestPackageTestContexts)
        );

        craneTestDependencyArtifacts = lib.fix (
          self:
          lib.mapAttrs (
            _: context:
            if context.componentPackages == [ workspaceHackPackage ] then
              workspaceHackBuildArtifacts
            else
              let
                sourcePackages = lib.filter workspacePackageHasLibTarget context.sourcePackages;
                dependencyContexts = workspaceTestDependencyContextsFor context;
                procMacroPackages = lib.filter workspacePackageIsProcMacro context.componentPackages;
                hostAnchorPackage = "${testBinaryFeatureAnchorPackageFor context}-host";
                hostAnchorPackageDir = "crates/${hostAnchorPackage}";
                anchorConsumerDependencies = {
                  ${testBinaryFeatureAnchorPackageFor context}.path =
                    "../${testBinaryFeatureAnchorPackageFor context}";
                }
                // lib.listToAttrs (
                  map (
                    package:
                    let
                      features = lib.filter (
                        feature: feature != "default" && !lib.hasInfix "/" feature
                      ) context.features.${package};
                    in
                    {
                      name = package;
                      value = {
                        path = workspacePrebuildDependencyPathFor package;
                      }
                      // lib.optionalAttrs (features != [ ]) { inherit features; };
                    }
                  ) procMacroPackages
                );
                hostAnchorCargoToml = (pkgs.formats.toml { }).generate "${hostAnchorPackage}-Cargo.toml" {
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
                targetFacadeFor =
                  package:
                  let
                    manifest = workspaceManifestFor workspaceMemberPackageDirs.${package};
                    packageName = "${testBinaryFeatureAnchorPackageFor context}-target-${package}";
                    packageDir = "crates/${packageName}";
                    localFeatures = lib.filter (
                      feature: feature != "default" && !lib.hasInfix "/" feature
                    ) context.features.${package};
                    resolvedDependencies = testBinaryFeatureAnchorDependenciesFor context;
                    rawDependencies = (manifest.dependencies or { }) // (manifest."dev-dependencies" or { });
                    dependencies =
                      lib.mapAttrs (
                        dependency: rawDependency:
                        let
                          dependencyAttrs =
                            if builtins.isAttrs rawDependency then rawDependency else { version = rawDependency; };
                          features = sortedUnique (
                            (dependencyAttrs.features or [ ]) ++ (resolvedDependencies.${dependency}.features or [ ])
                          );
                        in
                        dependencyAttrs // lib.optionalAttrs (features != [ ]) { inherit features; }
                      ) rawDependencies
                      // anchorConsumerDependencies;
                    cargoToml = (pkgs.formats.toml { }).generate "${packageName}-Cargo.toml" {
                      package = {
                        name = packageName;
                        version = "0.1.0";
                        edition = "2024";
                        publish = false;
                      };
                      lib.path = "src/lib.rs";
                      inherit dependencies;
                      features = manifest.features or { };
                    };
                  in
                  {
                    inherit
                      cargoToml
                      localFeatures
                      packageDir
                      packageName
                      ;
                  };
                targetFacades = map targetFacadeFor procMacroPackages;
                targetFacadeFeatures = lib.concatMap (
                  facade: map (feature: "${facade.packageName}/${feature}") facade.localFeatures
                ) targetFacades;
                anchorConsumerArgs = lib.concatStringsSep " " (
                  [
                    "--offline"
                    "-p ${lib.escapeShellArg hostAnchorPackage}"
                  ]
                  ++ map (facade: "-p ${lib.escapeShellArg facade.packageName}") targetFacades
                  ++ lib.optional (targetFacadeFeatures != [ ]) (cargoFeatureArgs targetFacadeFeatures)
                );
                hasProcMacro = procMacroPackages != [ ];
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
              buildDepsOnlyWithArtifacts (
                (builtins.removeAttrs ciArgs [ "src" ])
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
                    '') targetFacades}
                  '';
                }
                // lib.optionalAttrs context.usesPythonModule {
                  nativeBuildInputs = (ciArgs.nativeBuildInputs or [ ]) ++ [ nextestPython ];
                  PYO3_PYTHON = "${nextestPython}/bin/python3";
                  PYTHON = "${nextestPython}/bin/python3";
                  PYTHONPATH = "${gammaloop-python-module}/${pythonSitePackages}:${nextestPython}/${pythonSitePackages}";
                }
              )
          ) workspaceTestContexts
        );

        craneTestBinaryArtifactFor =
          context: package:
          buildDepsOnlyWithArtifacts (
            (builtins.removeAttrs ciArgs [ "src" ])
            // {
              cargoArtifacts =
                mergeCargoArtifacts "gammaloop-crate-test-binaries-${package}-${context.key}-inputs"
                  [
                    cargoArtifacts
                    craneTestDependencyArtifacts.${context.key}
                  ];
              pname = "gammaloop-crate-test-binaries-${package}-${context.key}";
              dummySrc = workspacePackageSrcForSourcePackages {
                sourcePackages = context.sourcePackages;
                packageSourcePackages = [ package ];
                testSourcePackages = [ package ];
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
              nativeBuildInputs = (ciArgs.nativeBuildInputs or [ ]) ++ [ nextestPython ];
              PYO3_PYTHON = "${nextestPython}/bin/python3";
              PYTHON = "${nextestPython}/bin/python3";
              PYTHONPATH = "${gammaloop-python-module}/${pythonSitePackages}:${nextestPython}/${pythonSitePackages}";
            }
          );

        craneTestBinaryArtifacts = lib.genAttrs workspacePackages (
          package:
          if package == workspaceHackPackage then
            workspaceHackBuildArtifacts
          else
            craneTestBinaryArtifactFor (workspaceTestContextFor { packages = [ package ]; }) package
        );

        workspaceCargoCheck = craneLib.mkCargoDerivation (
          ciArgs
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
          }
        );

        workspaceClippyCheck = craneLib.mkCargoDerivation (
          ciArgs
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
          }
        );

        workspaceDocCheck = craneLib.mkCargoDerivation (
          ciArgs
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
          }
        );

        alphal00pDocsDerivationArgs = alphal00pDocsRealCargoArgs // {
          cargoArtifacts = alphal00pDocsCargoArtifacts;
          src = documentationSrc;
          nativeBuildInputs = (alphal00pDocsCargoArgs.nativeBuildInputs or [ ]) ++ [
            docsTypst
            pkgs.gitMinimal
            pkgs.maturin
            pkgs.roboto
            pkgs.uv
          ];
          TYPST_FONT_PATHS = docsFontPath;
          ALPHAL00P_DOCS_CARGO_PROFILE = docsCargoProfile;
          doNotLinkInheritedArtifacts = true;
          doInstallCargoArtifacts = false;
          checkPhaseCargoCommand = "";
          doCheck = false;
          installPhaseCommand = "";
        };

        alphal00pDocsPagesChannel =
          let
            requested = builtins.getEnv "ALPHAL00P_DOCS_CHANNEL";
            channel = if requested == "" then "latest" else requested;
          in
          assert lib.assertMsg (builtins.elem channel [
            "latest"
            "snapshot"
          ]) "ALPHAL00P_DOCS_CHANNEL must be latest or snapshot";
          channel;
        alphal00pDocsPagesSnapshotTag =
          let
            tag = builtins.getEnv "ALPHAL00P_DOCS_SNAPSHOT_TAG";
          in
          assert lib.assertMsg (
            (alphal00pDocsPagesChannel == "latest" && tag == "")
            || (alphal00pDocsPagesChannel == "snapshot" && tag != "")
          ) "ALPHAL00P_DOCS_SNAPSHOT_TAG must be set exactly for snapshot builds";
          tag;
        alphal00pDocsPagesGitCommit =
          let
            commit = builtins.getEnv "ALPHAL00P_DOCS_GIT_COMMIT";
          in
          if commit == "" then nixCiBarrierRevision else commit;
        alphal00pDocsPagesGitTimestamp =
          let
            timestamp = builtins.getEnv "ALPHAL00P_DOCS_GIT_TIMESTAMP";
          in
          if timestamp == "" then toString self.lastModified else timestamp;

        alphal00pDocsValidationCommands = ''
          cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-catalogs --features ${lib.escapeShellArg documentationCatalogFeatures} --bin alphal00p-docs-gammaloop-reference -- --check
          cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-catalogs --features ${lib.escapeShellArg documentationCatalogFeatures} --bin alphal00p-docs-vakint-reference -- --check
          cargo test --locked --profile ${docsCargoProfile} -p alphal00p-docs-examples
          cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features gammaloop -- gammaloop-python docs/api/python/gammaloop-python.pyi --check
          cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features linnet -- linnet-py docs/api/python/linnet-py.pyi --check
          cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features spenso -- spynso3 docs/api/python/spynso3.pyi --check
          cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features idenso -- idenso-community docs/api/python/idenso-community.pyi --check
          cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features vakint -- vakint-community docs/api/python/vakint-community.pyi --check
          cargo test --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features gammaloop gammaloop_runtime_surface_and_signatures_match_the_docs_stub
          linnet_python="$TMPDIR/alphal00p-docs-linnet-python"
          export UV_CACHE_DIR="$TMPDIR/alphal00p-docs-uv-cache"
          uv venv "$linnet_python" --python "$PYO3_PYTHON"
          VIRTUAL_ENV="$linnet_python" maturin develop \
            --uv \
            --offline \
            --locked \
            --profile ${docsCargoProfile} \
            --manifest-path crates/linnet-py/Cargo.toml \
            --features extension-module,abi3-py310
          "$linnet_python/bin/python" -m unittest crates/linnet-py/tests/test_basic.py
          cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-builder -- check
          svg_assets="$TMPDIR/alphal00p-svg-assets"
          bash scripts/render-docs-svg-assets.sh "$svg_assets"
          checked_assets=(
            docs/assets/about-*.svg
            docs/assets/graphs/portal-*.svg
            docs/assets/local-unitarity-*.svg
            docs/assets/spensologo.svg
            assets/gammalooplogo*.svg
          )
          generated_assets=(
            "$svg_assets"/docs/assets/about-*.svg
            "$svg_assets"/docs/assets/graphs/portal-*.svg
            "$svg_assets"/docs/assets/local-unitarity-*.svg
            "$svg_assets"/docs/assets/spensologo.svg
            "$svg_assets"/assets/gammalooplogo*.svg
          )
          test "''${#checked_assets[@]}" -eq 32
          test "''${#generated_assets[@]}" -eq 32
          for checked_asset in "''${checked_assets[@]}"; do
            cmp "$checked_asset" "$svg_assets/$checked_asset"
          done
        '';

        alphal00pDocsDeveloperAssertions = output: ''
          test -s "${output}/developers/architecture/documentation-improvement-plan/index.html"
          plan_page="${output}/developers/architecture/documentation-improvement-plan/index.html"
          grep -Fq 'DOC-010' "$plan_page"
          grep -Fq 'id="executive-assessment"' "$plan_page"
          grep -Fq 'href="#executive-assessment"' "$plan_page"
          grep -Fq 'documentation-improvement-plan.typ' "$plan_page"
          test "$(grep -o '<html' "$plan_page" | wc -l)" -eq 1
          test "$(grep -o '<body' "$plan_page" | wc -l)" -eq 1
          grep -Fq '"title": "Executive assessment"' "${output}/developers/search-index.json"
          grep -Fq '"href": "architecture/documentation-improvement-plan/#executive-assessment"' "${output}/developers/search-index.json"
          for product in gammaloop linnet spenso idenso vakint; do
            test -s "${output}/developers/architecture/$product-architecture/index.html"
          done
          test -s "${output}/developers/architecture/spenso-parsing-flow/index.html"
          test ! -e "${output}/developers/architecture/spenso-parsing-flow/diagram.html"
        '';

        alphal00pDocsCheck = craneLib.mkCargoDerivation (
          alphal00pDocsDerivationArgs
          // {
            pname = "alphal00p-docs";
            ALPHAL00P_DOCS_GIT_COMMIT = nixCiBarrierRevision;
            ALPHAL00P_DOCS_GIT_TIMESTAMP = toString self.lastModified;
            buildPhaseCargoCommand = ''
              docs_rustdoc=${lib.escapeShellArg alphal00pDocsCargoTargetRoot}

              ${alphal00pDocsValidationCommands}

              docs_first="$TMPDIR/alphal00p-docs-first"
              docs_second="$TMPDIR/alphal00p-docs-second"
              docs_snapshot="$TMPDIR/alphal00p-docs-snapshot"
              cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-builder -- \
                build \
                --product all \
                --channel latest \
                --output "$docs_first" \
                --rustdoc-target-root "$docs_rustdoc"
              cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-builder -- \
                build \
                --product all \
                --channel latest \
                --output "$docs_second" \
                --rustdoc-target-root "$docs_rustdoc"
              diff --recursive --brief "$docs_first" "$docs_second"
              python3 scripts/check-docs-html.py "$docs_first"

              cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-builder -- \
                build \
                --product all \
                --channel snapshot \
                --snapshot-tag v0.3.4 \
                --output "$docs_snapshot" \
                --rustdoc-target-root "$docs_rustdoc"
              mkdir -p "$docs_snapshot/developers" \
                "$docs_snapshot/.staging" \
                "$docs_snapshot/products/gammaloop/latest" \
                "$docs_snapshot/products/gammaloop/snapshots/legacy"
              printf 'removed developer route\n' > "$docs_snapshot/developers/removed-before-rebuild.txt"
              printf 'stale staging\n' > "$docs_snapshot/.staging/incomplete.txt"
              printf 'latest route\n' > "$docs_snapshot/products/gammaloop/latest/.note"
              printf 'product redirect\n' > "$docs_snapshot/products/gammaloop/index.html"
              printf 'historical snapshot\n' > "$docs_snapshot/products/gammaloop/snapshots/legacy/.note"
              cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-builder -- \
                build \
                --product all \
                --channel snapshot \
                --snapshot-tag v0.3.4 \
                --output "$docs_snapshot" \
                --rustdoc-target-root "$docs_rustdoc"
              test ! -e "$docs_snapshot/developers/removed-before-rebuild.txt"
              test ! -e "$docs_snapshot/.staging"
              grep -Fq 'latest route' "$docs_snapshot/products/gammaloop/latest/.note"
              grep -Fq 'product redirect' "$docs_snapshot/products/gammaloop/index.html"
              grep -Fq 'historical snapshot' "$docs_snapshot/products/gammaloop/snapshots/legacy/.note"

              docs_pages_test="$TMPDIR/alphal00p-docs-pages-test"
              mkdir -p "$docs_pages_test/products/gammaloop/snapshots/legacy"
              mkdir -p "$docs_pages_test/developers"
              mkdir -p "$docs_pages_test/.git"
              printf 'historical snapshot\n' > "$docs_pages_test/products/gammaloop/snapshots/legacy/.note"
              printf 'removed developer route\n' > "$docs_pages_test/developers/old.txt"
              mv "$docs_first/search-index.json" "$TMPDIR/federated-search-index.json"
              if bash scripts/update-docs-pages.sh latest "$docs_first" "$docs_pages_test" \
                2> "$TMPDIR/missing-federated-search-index.err"; then
                missing_search_status=0
              else
                missing_search_status=$?
              fi
              mv "$TMPDIR/federated-search-index.json" "$docs_first/search-index.json"
              test "$missing_search_status" -ne 0
              grep -Fq 'latest build has no federated search index' \
                "$TMPDIR/missing-federated-search-index.err"
              bash scripts/update-docs-pages.sh latest "$docs_first" "$docs_pages_test"
              cmp "$docs_pages_test/search-index.json" "$docs_first/search-index.json"
              test -s "$docs_pages_test/products/gammaloop/snapshots/legacy/.note"
              test ! -e "$docs_pages_test/developers/old.txt"
              test -s "$docs_pages_test/developers/architecture/gammaloop-architecture/index.html"
              test -s "$docs_pages_test/developers/architecture/documentation-improvement-plan/index.html"
              cp "$docs_pages_test/index.html" "$TMPDIR/portal-before-snapshot.html"
              cp "$docs_pages_test/developers/.note" "$TMPDIR/developers-before-snapshot.note"
              cp "$docs_pages_test/products/gammaloop/latest/.note" "$TMPDIR/latest-before-snapshot.note"
              bash scripts/update-docs-pages.sh snapshot "$docs_snapshot" "$docs_pages_test" v0.3.4
              cmp "$docs_pages_test/index.html" "$TMPDIR/portal-before-snapshot.html"
              cmp "$docs_pages_test/developers/.note" "$TMPDIR/developers-before-snapshot.note"
              cmp "$docs_pages_test/products/gammaloop/latest/.note" "$TMPDIR/latest-before-snapshot.note"
              for product in gammaloop linnet spenso idenso vakint; do
                test -s "$docs_pages_test/products/$product/snapshots/v0.3.4/.note"
              done

              mkdir -p "$out"
              cp -R "$docs_first"/. "$out"/

              test -s "$out/index.html"
              test -e "$out/.nojekyll"
              test -s "$out/assets/site.css"
              test -s "$out/assets/site.js"
              test -s "$out/assets/local-unitarity-light.svg"
              test -s "$out/assets/local-unitarity-dark.svg"
              test -s "$out/assets/gammalooplogo-light.svg"
              test -s "$out/assets/gammalooplogo-dark.svg"
              test -s "$out/developers/.note"
              test -s "$out/developers/index.html"
              test -s "$out/developers/search-index.json"
              test -s "$out/developers/assets/site.css"
              test -s "$out/developers/assets/site.js"
              test -s "$out/developers/architecture/gammaloop-architecture/index.html"
              ${alphal00pDocsDeveloperAssertions "$out"}
              for product in gammaloop linnet spenso idenso vakint; do
                product_root="$out/products/$product"
                test -s "$product_root/index.html"
                test -s "$product_root/latest/index.html"
                test -s "$product_root/latest/.note"
                test -s "$product_root/latest/manual.pdf"
                test -s "$product_root/latest/search-index.json"
                test -s "$product_root/latest/snapshot.json"
                test -s "$product_root/latest/tutorial/index.html"
                test -s "$product_root/latest/reference/interfaces/index.html"
                test -s "$product_root/latest/version-history/index.html"
                test -s "$product_root/latest/manual/interfaces/index.html"
                test -s "$product_root/latest/manual/releases/index.html"
                grep -Fq 'url=../../reference/interfaces/' \
                  "$product_root/latest/manual/interfaces/index.html"
                grep -Fq 'url=../../version-history/' \
                  "$product_root/latest/manual/releases/index.html"
                test -s "$product_root/latest/assets/site.css"
                test -s "$product_root/latest/assets/site.js"
                test -s "$product_root/latest/assets/local-unitarity-light.svg"
                test -s "$product_root/latest/assets/local-unitarity-dark.svg"
                test -s "$product_root/latest/assets/gammalooplogo-light.svg"
                test -s "$product_root/latest/assets/gammalooplogo-dark.svg"
                test -s "$product_root/latest/reference/rust/index.html"
                test -s "$product_root/latest/reference/python/index.html"
                test -n "$(find "$product_root/latest/reference/python" \
                  -mindepth 2 -maxdepth 2 -name index.html -print -quit)"
                ! grep -q "Rustdoc generation was skipped" \
                  "$product_root/latest/reference/rust/index.html"
              done
              test -s "$out/products/gammaloop/latest/reference/rust/gammalooprs/index.html"
              test -s "$out/products/gammaloop/latest/reference/rust/gammaloop_api/index.html"
              test -s "$out/products/linnet/latest/reference/rust/linnet/index.html"
              test -s "$out/products/spenso/latest/reference/rust/spenso/index.html"
              test -s "$out/products/spenso/latest/reference/rust/spenso_macros/index.html"
              test -s "$out/products/spenso/latest/reference/rust/spenso_hep_lib/index.html"
              test -s "$out/products/idenso/latest/reference/rust/idenso/index.html"
              test -s "$out/products/vakint/latest/reference/rust/vakint/index.html"
              test -s "$out/products/gammaloop/latest/reference/rust/theme.css"
              grep -Fq 'href="../theme.css"' \
                "$out/products/gammaloop/latest/reference/rust/gammalooprs/index.html"
              grep -Fq 'class="alphal00p-rustdoc-bar"' \
                "$out/products/gammaloop/latest/reference/rust/gammalooprs/index.html"
              test -s \
                "$out/products/gammaloop/latest/reference/rust/src/gammalooprs/lib.rs.html"
              grep -Fq 'href="../../theme.css"' \
                "$out/products/gammaloop/latest/reference/rust/src/gammalooprs/lib.rs.html"
              grep -Fq 'class="alphal00p-rustdoc-bar"' \
                "$out/products/gammaloop/latest/reference/rust/src/gammalooprs/lib.rs.html"
              test -s "$out/products/linnet/latest/reference/typst/index.html"
              for typst_page in graph layout drawing templates subgraph; do
                test -s \
                  "$out/products/linnet/latest/reference/typst/$typst_page/index.html"
              done
              grep -Fq '"title": "graph.build"' \
                "$out/products/linnet/latest/search-index.json"
              ! grep -Fq '"title": "Parameters"' \
                "$out/products/linnet/latest/search-index.json"
              grep -Fq '../../../../gammaloop/latest/guides/kurvst/' \
                "$out/products/linnet/latest/reference/typst/index.html"
            '';
          }
        );

        alphal00pDocsPages = craneLib.mkCargoDerivation (
          alphal00pDocsDerivationArgs
          // {
            pname = "alphal00p-docs-pages";
            ALPHAL00P_DOCS_CHANNEL = alphal00pDocsPagesChannel;
            ALPHAL00P_DOCS_GIT_COMMIT = alphal00pDocsPagesGitCommit;
            ALPHAL00P_DOCS_GIT_TIMESTAMP = alphal00pDocsPagesGitTimestamp;
            ALPHAL00P_DOCS_SNAPSHOT_TAG = alphal00pDocsPagesSnapshotTag;
            buildPhaseCargoCommand = ''
              ${alphal00pDocsValidationCommands}

              docs_args=(
                build
                --product all
                --channel ${lib.escapeShellArg alphal00pDocsPagesChannel}
                --output "$out"
                --rustdoc-target-root ${lib.escapeShellArg alphal00pDocsCargoTargetRoot}
              )
              ${lib.optionalString (alphal00pDocsPagesChannel == "snapshot") ''
                docs_args+=(--snapshot-tag ${lib.escapeShellArg alphal00pDocsPagesSnapshotTag})
              ''}
              cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-builder -- "''${docs_args[@]}"
              python3 scripts/check-docs-html.py "$out"

              test -s "$out/index.html"
              test -e "$out/.nojekyll"
              ${alphal00pDocsDeveloperAssertions "$out"}
            '';
          }
        );

        workspaceDoctestCheck = craneLib.mkCargoDerivation (
          ciArgs
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
              cargoWithProfile test --doc ${ciArgs.cargoExtraArgs} \
                --exclude alphal00p-docs-python-exporter
            '';
            checkPhaseCargoCommand = "";
            doCheck = false;
            installPhaseCommand = "";
            preBuild = licensePreCheck;
            SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
          }
        );

        cranePythonDependencyArtifacts = buildDepsOnlyWithArtifacts (
          (builtins.removeAttrs ciArgs [ "src" ])
          // {
            cargoArtifacts = mergeCargoArtifacts "gammaloop-python-deps-inputs" (
              [ cargoArtifacts ]
              ++
                map
                  (
                    dependency:
                    if dependency == workspaceHackPackage then
                      workspaceHackBuildArtifacts
                    else
                      cranePackageDependencyModeArtifacts.${dependency}
                  )
                  (
                    lib.filter (sourcePackage: sourcePackage != "gammaloop-api") (
                      workspaceNormalSourcePackageNamesFor "gammaloop-api"
                    )
                  )
            );
            pname = "gammaloop-api-python";
            src = workspacePackageSrcForSourcePackages {
              sourcePackages = lib.filter (sourcePackage: sourcePackage != "gammaloop-api") (
                workspaceNormalSourcePackageNamesFor "gammaloop-api"
              );
            };
            buildPhaseCargoCommand = "cargoWithProfile build ${cranePythonCargoArgs}";
            checkPhaseCargoCommand = "";
            doCheck = false;
            preBuildWorkspaceArtifactStripPackages = [ "gammaloop-api" ];
            stripWorkspaceArtifacts = true;
            preservedWorkspaceArtifactPackages = workspaceResolvedDependencyNamesFor "gammaloop-api";
            extraDummyScript = workspaceDependencyDummyCargoTargetsScriptFor "gammaloop-api";
            postPatch = workspaceMissingCargoTargetsScript;
          }
        );

        cranePythonBuildArtifacts = craneLib.cargoBuild (
          ciArgs
          // {
            cargoArtifacts = cranePythonDependencyArtifacts;
            pname = "gammaloop-api-python-build";
            src = workspacePackageSrcFor "gammaloop-api";
            cargoExtraArgs = cranePythonCargoArgs;
            postPatch = workspaceMissingCargoTargetsScript;
          }
        );

        cranePackageOutputs = lib.listToAttrs (
          map (package: {
            name = "crate-${package}";
            value = cranePackageBuildArtifacts.${package};
          }) workspacePackages
        );

        cranePackageDependencyOutputs = lib.listToAttrs (
          map (package: {
            name = "crate-deps-${package}";
            value = nixCiArtifactBarrier "crate-deps-${package}" cranePackageDependencyArtifacts.${package};
          }) (lib.filter (package: cranePackageDependencyArtifacts.${package} != null) workspacePackages)
        );

        craneTestBinaryPackageOutputs = lib.listToAttrs (
          map (package: {
            name = "crate-test-binaries-${package}";
            value = nixCiArtifactBarrier "crate-test-binaries-${package}" craneTestBinaryArtifacts.${package};
          }) workspacePackages
        );

        craneTestDependencyOutputs = lib.listToAttrs (
          map (
            representative:
            let
              context = workspaceTestContextFor { packages = [ representative ]; };
            in
            {
              name = "crate-test-dependencies-${representative}";
              value =
                nixCiArtifactBarrier "crate-test-dependencies-${representative}"
                  craneTestDependencyArtifacts.${context.key};
            }
          ) workspaceTestDependencyComponentRepresentatives
        );

        workspaceBuildArtifacts = craneLib.cargoBuild (
          ciArgs
          // {
            inherit cargoArtifacts;
            pname = "gammaloop-workspace-build-artifacts";
            src = workspaceNonIntegrationTestSrc;
            cargoExtraArgs = "${
              cargoPackagesArgsFor (lib.subtractLists [
                "gammaloop-integration-tests"
              ] workspaceMemberPackages) craneTestFeaturesFor
            } --tests";
          }
        );

        symbolicaCrateArgs =
          usesSymbolica:
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

        linnestWasmCargoArtifacts = wasmCraneLib.buildDepsOnly (
          linnestWasmArgs
          // {
            pname = "linnest-wasm-deps";
          }
        );

        linnest-wasm = wasmCraneLib.buildPackage (
          linnestWasmArgs
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
          }
        );

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
            name = "docs";
            packages = [
              "alphal00p-docs-builder"
              "alphal00p-docs-catalogs"
              "alphal00p-docs-examples"
              "alphal00p-docs-macros"
              "alphal00p-docs-python-exporter"
              "alphal00p-docs-schema"
            ];
            runtimeInputs = [
              nextestPython
              pkgs.gitMinimal
            ];
          }
          {
            name = "integration";
            packages = [ "gammaloop-integration-tests" ];
          }
          {
            name = "python-api";
            packages = [ "gammaloop-integration-tests" ];
            runtimeTestSourcePackages = [ ];
            filter = "package(gammaloop-integration-tests) & binary(test_python_api)";
            extraFeatures."gammaloop-integration-tests" = [ "python-api-tests" ];
          }
          {
            name = "clinnet";
            packages = [ "clinnet" ];
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
            packages = [ "vakint" ];
          }
        ];

        sortedUnique = list: lib.sort (left: right: left < right) (lib.unique list);

        nextestSplitPackages = sortedUnique (lib.concatMap (target: target.packages) nextestPackageGroups);
        workspacePackages = sortedUnique workspaceMemberPackages;
        nextestCoverageIgnoredPackages = [
          "gammaloop-workspace-hack"
          "spynso3"
        ];
        nextestCoveredWorkspacePackages = sortedUnique (
          lib.subtractLists nextestCoverageIgnoredPackages workspaceMemberPackages
        );
        missingNextestPackages = lib.subtractLists nextestSplitPackages nextestCoveredWorkspacePackages;
        extraNextestPackages = lib.subtractLists nextestCoveredWorkspacePackages nextestSplitPackages;

        checkedNextestPackageGroups =
          assert lib.asserts.assertMsg (missingNextestPackages == [ ] && extraNextestPackages == [ ])
            "nextest split package coverage mismatch: missing [${lib.concatStringsSep ", " missingNextestPackages}], extra [${lib.concatStringsSep ", " extraNextestPackages}]";
          nextestPackageGroups;

        nextestTargetTriple =
          pkgs.stdenv.hostPlatform.rust.rustcTargetSpec or pkgs.stdenv.hostPlatform.config;
        nextestRustLibDir = "${ciToolchain}/lib/rustlib/${nextestTargetTriple}/lib";

        nextestBaseExtraArgs = "--profile ${nextestProfile} --no-fail-fast --final-status-level fail --no-tests=pass";

        nextestPackageFilter =
          packages:
          "-E ${
            lib.escapeShellArg (lib.concatMapStringsSep " | " (package: "package(${package})") packages)
          }";
        nextestFilterFor =
          target:
          if target ? filter then
            "-E ${lib.escapeShellArg target.filter}"
          else
            nextestPackageFilter target.packages;
        nextestContextFor =
          target:
          workspaceTestContextFor {
            inherit (target) packages;
            extraFeatures = target.extraFeatures or { };
          };
        nextestPackageContextFor =
          target: package:
          workspaceTestContextFor {
            packages = [ package ];
            extraFeatures = target.extraFeatures or { };
          };
        nextestUsesPythonModule = target: (nextestContextFor target).usesPythonModule;
        nextestSourcePackagesFor = target: (nextestContextFor target).sourcePackages;
        nextestSrcFor =
          target:
          workspacePackageSrcForSourcePackages {
            sourcePackages = nextestSourcePackagesFor target;
            packageSourcePackages = target.packages;
            testSourcePackages = target.packages;
            extraFilesets = [ ./.config/nextest.toml ];
          };
        nextestRuntimeSrcFor =
          target:
          workspacePackageSrcForSourcePackages {
            sourcePackages = nextestSourcePackagesFor target;
            packageSourcePackages = target.packages;
            testSourcePackages = target.packages;
            runtimeTestSourcePackages = target.runtimeTestSourcePackages or target.packages;
            extraFilesets = [
              ./.config/nextest.toml
            ]
            ++ lib.optionals (target.name == "docs") documentationDeveloperScopeSources;
          };

        nextestFeatureArgsFor =
          target:
          cargoQualifiedFeatureArgsFor target.packages (
            package: (nextestPackageContextFor target package).features.${package}
          );

        nextestCargoArgsFor =
          target:
          lib.concatStringsSep " " (
            [
              "--offline"
            ]
            ++ map (package: "-p ${lib.escapeShellArg package}") target.packages
            ++ lib.optional (nextestFeatureArgsFor target != "") (nextestFeatureArgsFor target)
          );

        nextestArchiveNameFor = target: package: "gammaloop-nextest-${target.name}-${package}.tar.zst";

        nextestPackageTestBinaryArtifactFor =
          target: package:
          let
            context = nextestPackageContextFor target package;
            baseContext = workspaceTestContextFor { packages = [ package ]; };
          in
          if context.key == baseContext.key then
            craneTestBinaryArtifacts.${package}
          else
            craneTestBinaryArtifactFor context package;

        nextestPackageArchiveFor =
          target: package:
          let
            packageTarget = target // {
              packages = [ package ];
            };
            context = nextestPackageContextFor target package;
            packageCargoArtifacts = nextestPackageTestBinaryArtifactFor target package;
          in
          craneLib.mkCargoDerivation (
            ciArgs
            // {
              pname = "gammaloop-nextest-binaries-${target.name}-${package}";
              src = nextestSrcFor packageTarget;
              cargoArtifacts = packageCargoArtifacts;
              doCheck = false;
              doInstallCargoArtifacts = false;
              nativeBuildInputs =
                (ciArgs.nativeBuildInputs or [ ])
                ++ [
                  pkgs.cargo-nextest
                  pkgs.form
                ]
                ++ lib.optionals (nextestUsesPythonModule target) [ nextestPython ];
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
            }
            // lib.optionalAttrs (nextestUsesPythonModule target) {
              PYO3_PYTHON = "${nextestPython}/bin/python3";
              PYTHON = "${nextestPython}/bin/python3";
              PYTHONPATH = "${gammaloop-python-module}/${pythonSitePackages}:${nextestPython}/${pythonSitePackages}";
            }
          );

        nextestArchiveFor =
          target:
          let
            packageArchives = lib.genAttrs target.packages (nextestPackageArchiveFor target);
          in
          pkgs.linkFarm "gammaloop-nextest-binaries-${target.name}" (
            map (package: {
              name = nextestArchiveNameFor target package;
              path = "${packageArchives.${package}}/${nextestArchiveNameFor target package}";
            }) target.packages
          );

        nextestBinarySets = lib.listToAttrs (
          map (target: {
            name = "gammaloop-nextest-binaries-${target.name}";
            value = nextestArchiveFor target;
          }) checkedNextestPackageGroups
        );

        nextestContextualTestOutputs = lib.listToAttrs (
          lib.concatMap (
            target:
            lib.concatMap (
              package:
              let
                context = nextestPackageContextFor target package;
              in
              [
                {
                  name = "crate-test-dependencies-${target.name}-${package}";
                  value =
                    nixCiArtifactBarrier "crate-test-dependencies-${target.name}-${package}"
                      craneTestDependencyArtifacts.${context.key};
                }
                {
                  name = "crate-test-binaries-${target.name}-${package}";
                  value = nixCiArtifactBarrier "crate-test-binaries-${target.name}-${package}" (
                    nextestPackageTestBinaryArtifactFor target package
                  );
                }
              ]
            ) target.packages
          ) (lib.filter (target: target ? extraFeatures) checkedNextestPackageGroups)
        );

        nextestBinarySetForTarget = target: nextestBinarySets."gammaloop-nextest-binaries-${target.name}";

        nextestBinarySetAggregate = pkgs.linkFarm "gammaloop-nextest-binaries" (
          map (target: {
            name = target.name;
            path = nextestBinarySetForTarget target;
          }) checkedNextestPackageGroups
        );

        nextestCheckFor =
          target:
          pkgs.runCommand "gammaloop-nextest-${target.name}"
            {
              nativeBuildInputs = [
                ciToolchain
                pkgs.cargo-nextest
                pkgs.form
                pkgs.gcc
                nextestFailureSummary
              ]
              ++ (target.runtimeInputs or [ ])
              ++ lib.optionals (nextestUsesPythonModule target) [ nextestPython ];
              CC = nixCc;
              CXX = nixCxx;
              "${cargoLinkerVar}" = nixCc;
              LD_LIBRARY_PATH = runtimeLibPath;
              DYLD_LIBRARY_PATH = runtimeLibPath;
              NEXTEST_SHOW_PROGRESS = "counter";
              RUST_BACKTRACE = "1";
              RUST_LIB_BACKTRACE = "1";
              SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
            }
            (
              ''
                if [ -z "''${SYMBOLICA_LICENSE:-}" ]; then
                  echo "Missing SYMBOLICA_LICENSE environment variable" >&2
                  exit 1
                fi

                mkdir -p /build/source
                cp -R ${nextestRuntimeSrcFor target}/. /build/source/
                chmod -R u+w /build/source
                cd /build/source
                ${lib.optionalString (target.name == "docs") ''
                  export ALPHAL00P_DOCS_GIT_COMMIT=${lib.escapeShellArg nixCiBarrierRevision}
                  export ALPHAL00P_DOCS_GIT_TIMESTAMP=${lib.escapeShellArg (toString self.lastModified)}
                  mkdir -p "$NIX_BUILD_TOP/nextest-tmp"
                  export TMPDIR="$NIX_BUILD_TOP/nextest-tmp"
                ''}
                mkdir -p .cargo .cargo-home
                cp ${cargoVendorDir}/config.toml .cargo/config.toml
                export CARGO_HOME="$PWD/.cargo-home"
                # Workspace-root discovery checks these directories even for test
                # targets that do not consume any files from them.
                mkdir -p tests/resources examples/cli
                ${workspaceMissingCargoTargetsScript}
              ''
              + lib.optionalString (nextestUsesPythonModule target) ''
                export PYO3_PYTHON=${nextestPython}/bin/python3
                export PYTHON=${nextestPython}/bin/python3
                export PYTHONPATH=${gammaloop-python-module}/${pythonSitePackages}:${nextestPython}/${pythonSitePackages}
              ''
              + ''
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
                    --extract-to . \
                    --extract-overwrite \
                    --workspace-remap . \
                    ${nextestBaseExtraArgs}
                  package_status=$?
                  nextest-failure-summary ${lib.escapeShellArg nextestJunitPath} || true
                  if [ "$package_status" -ne 0 ] && [ "$status" -eq 0 ]; then
                    status="$package_status"
                  fi
                '') target.packages}
                mkdir -p "$out"
                exit "$status"
              ''
            );

        nextestRunChecks = lib.listToAttrs (
          map (target: {
            name = "gammaloop-nextest-${target.name}";
            value = nextestCheckFor target;
          }) checkedNextestPackageGroups
        );

        nextestChecks = {
          gammaloop-nextest-binaries = nextestBinarySetAggregate;
        }
        // nextestBinarySets
        // nextestRunChecks;

        nextestAggregate = pkgs.runCommand "gammaloop-nextest" { } ''
          ${lib.concatMapStringsSep "\n" (check: "test -d ${check}") (builtins.attrValues nextestRunChecks)}
          mkdir -p "$out"
        '';

        impureCheckRunnerTargets = [
          {
            runnerAttr = "nix-ci-check-gammaloop-doctest";
            checkAttr = "gammaloop-doctest";
            runtimeInputs = [ cargoArtifacts ];
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
          runtimeInputs = [
            (nextestBinarySetForTarget target)
          ]
          ++ lib.optionals (nextestUsesPythonModule target) [ gammaloop-python-module ];
        }) checkedNextestPackageGroups;

        impureCheckRunnerPackages = lib.listToAttrs (
          map (target: {
            name = target.runnerAttr;
            value = pkgs.writeShellApplication {
              name = target.runnerAttr;
              # Retain the pure test inputs even when NixCI skips this runner as
              # cached, so the in-repo test does not rebuild them.
              runtimeInputs = [ pkgs.nix ] ++ target.runtimeInputs;
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
          }) impureCheckRunnerTargets
        );

        nixCiPassed = pkgs.writeShellApplication {
          name = "nix-ci-passed";
          text = ''
            echo "All NixCI build and test jobs passed."
          '';
        };

        allChecks = {
          # Keep existing check names for CI compatibility.
          gammaloop = gammaloop-cli;

          gammaloop-check = workspaceCargoCheck;

          gammaloop-clippy = workspaceClippyCheck;

          gammaloop-doc = workspaceDocCheck;

          alphal00p-docs = alphal00pDocsCheck;

          gammaloop-doctest = workspaceDoctestCheck;

          gammaloop-fmt = craneLib.cargoFmt {
            src = workspaceFmtSrc;
            pname = "gammaloop-workspace";
            inherit (apiMeta) version;
          };

          gammaloop-guppy-workspace-graph = guppyWorkspaceGraphCheck;

          linnest-wasm =
            pkgs.runCommand "linnest-wasm-check"
              {
                nativeBuildInputs = [ pkgs.wasm-tools ];
              }
              ''
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
      in
      {
        checks = allChecks;

        hydraJobs = hestiaChecks;

        packages = {
          default = gammaloop-cli;
          clinnet = clinnet-cli;
          gammaloop = gammaloop-cli;
          inherit clinnet-cli;
          "gammaloop-python-module" = nixCiArtifactBarrier "gammaloop-python-module" gammaloop-python-module;
          "alphal00p-docs-cargo-artifacts" = alphal00pDocsCargoArtifacts;
          "alphal00p-docs-pages" = alphal00pDocsPages;
          "ci-workspace-graph-json" = guppyWorkspaceGraphJson;
          inherit linnest-wasm;
          linnestWasmCargoArtifacts = nixCiArtifactBarrier "linnest-wasm-cargo-artifacts" linnestWasmCargoArtifacts;
          "crane-ci-prebuild" = cargoArtifacts;
          cargoArtifacts = nixCiArtifactBarrier "cargo-artifacts" cargoArtifacts;
          gammaloopApiPackageArtifacts = nixCiArtifactBarrier "gammaloop-api-package-artifacts" gammaloopApiPackageArtifacts;
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
          gammaloop-llvm-coverage = craneLib.cargoLlvmCov (
            commonArgs
            // {
              src = workspaceTestSrc;
              inherit cargoArtifacts;
              nativeBuildInputs = (commonArgs.nativeBuildInputs or [ ]) ++ [ pkgs.form ];
            }
          );
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
          default = mkDevShell [ clinnet-cli ];
          full = mkDevShell [
            clinnet-cli
            rscls
          ];
          clinnet = mkDevShell [ clinnet-cli ];
        };
      }
    );
}
