let
  system = "x86_64-linux";
  workspaceGraph = builtins.fromJSON (builtins.readFile ./nix/ci-workspace-graph.json);
  unique = values:
    builtins.attrNames (builtins.listToAttrs (map (value: {
        name = value;
        value = true;
      })
      values));
  workspacePackages = workspaceGraph.packages;
  cratePackageDepsAttr = package: "packages.${system}.crate-deps-${package}";
  cratePackageAttr = package: "packages.${system}.crate-${package}";
  crateTestSupportAttr = representative: "packages.${system}.crate-test-support-${representative}";
  crateTestBinaryAttr = package: "packages.${system}.crate-test-binaries-${package}";
  workspaceHackPackage = "gammaloop-workspace-hack";
  # clinnet is binary-only: the flake exposes crate-clinnet, but no
  # crate-deps-clinnet artifact.
  workspacePackagesWithDependencyArtifacts = builtins.filter (package: package != "clinnet") workspacePackages;
  nonWorkspaceHackPackages = builtins.filter (package: package != workspaceHackPackage) workspacePackages;
  workspaceHackCacheAttr = cratePackageDepsAttr workspaceHackPackage;
  gammaloopApiPackageArtifactsAttr = "packages.${system}.gammaloopApiPackageArtifacts";
  workspacePackageGraphAttr = package: cratePackageAttr package;
  mergeDependencySets = sets: let
    attrs = unique (builtins.concatLists (map builtins.attrNames sets));
  in
    builtins.listToAttrs (map (attr: {
        name = attr;
        value = unique (builtins.concatLists (map (set: set.${attr} or []) sets));
      })
      attrs);
  workspaceDependencyNamesFor = package:
    workspaceGraph.resolved_normal_dependencies.${package} or [];
  workspaceTestDependencyNamesFor = package:
    workspaceGraph.test_dependencies.${package} or (workspaceGraph.normal_dependencies.${package} or []);
  workspaceDependencyClosureFor = dependencyNamesFor: package:
    unique (map (entry: entry.key) (builtins.genericClosure {
      startSet = [{key = package;}];
      operator = entry: map (dependency: {key = dependency;}) (dependencyNamesFor entry.key);
    }));
  workspaceTestDependencyClosureFor =
    workspaceDependencyClosureFor workspaceTestDependencyNamesFor;
  workspaceTestComponentMembersFor = package:
    unique (builtins.filter (
        other:
          builtins.elem other (workspaceTestDependencyClosureFor package)
          && builtins.elem package (workspaceTestDependencyClosureFor other)
      )
      workspacePackages);
  workspaceTestComponentRepresentativeFor = package:
    builtins.head (workspaceTestComponentMembersFor package);
  workspaceTestComponentRepresentatives =
    unique (map workspaceTestComponentRepresentativeFor workspacePackages);
  workspaceTestSupportComponentRepresentatives =
    builtins.filter (representative: representative != workspaceHackPackage) workspaceTestComponentRepresentatives;
  workspaceTestComponentMembers =
    builtins.listToAttrs (map (representative: {
        name = representative;
        value = workspaceTestComponentMembersFor representative;
      })
      workspaceTestComponentRepresentatives);
  workspaceTestComponentDependencyRepresentativesFor = representative:
    unique (builtins.filter (dependencyRepresentative: dependencyRepresentative != representative) (
      map workspaceTestComponentRepresentativeFor (
        builtins.concatLists (map workspaceTestDependencyNamesFor workspaceTestComponentMembers.${representative})
      )
    ));
  workspaceCratePackageDependencies = builtins.listToAttrs (
    builtins.filter (entry: entry.value != []) (map (package: {
        name = workspacePackageGraphAttr package;
        value = map workspacePackageGraphAttr (workspaceDependencyNamesFor package);
      })
      workspacePackages)
  );
  workspaceCratePackageDependencyEdges = builtins.concatLists (map (dependent:
      map (dependency: {
        inherit dependency dependent;
      })
      (workspaceCratePackageDependencies.${dependent} or []))
    (builtins.attrNames workspaceCratePackageDependencies));
  workspaceCratePackageCacheDependencies = builtins.listToAttrs (map (package: {
      name = workspacePackageGraphAttr package;
      value = [(cratePackageDepsAttr package)];
    })
    workspacePackagesWithDependencyArtifacts);
  workspaceCratePackageCacheArtifactDependencies = builtins.listToAttrs (
    builtins.filter (entry: entry.value != []) (map (package: {
        name = cratePackageDepsAttr package;
        value =
          map cratePackageDepsAttr (workspaceDependencyNamesFor package)
          ++ (
            if package != workspaceHackPackage && (workspaceDependencyNamesFor package) == [] && builtins.elem package workspaceGraph.symbolica_normal_packages
            then [workspaceHackCacheAttr]
            else []
          );
      })
      workspacePackages)
  );
  workspaceCratePackageCacheArtifactDependencyEdges = builtins.concatLists (map (dependent:
      map (dependency: {
        inherit dependency dependent;
      })
      (workspaceCratePackageCacheArtifactDependencies.${dependent} or []))
    (builtins.attrNames workspaceCratePackageCacheArtifactDependencies));
  workspaceTestSupportDependencies = builtins.listToAttrs (map (representative: {
      name = crateTestSupportAttr representative;
      value =
        [(workspacePackageGraphAttr workspaceHackPackage)]
        ++ map crateTestSupportAttr (
          builtins.filter (
            dependencyRepresentative: dependencyRepresentative != workspaceHackPackage
          )
          (workspaceTestComponentDependencyRepresentativesFor representative)
        );
    })
    workspaceTestSupportComponentRepresentatives);
  nextestPackageGroups = {
    core = [
      "gammaloop-api"
      "gammaloop-tracing-filter"
      "gammaloop-tracing-filter-macros"
      "gammalooprs"
    ];
    integration = ["gammaloop-integration-tests"];
    "python-api" = ["gammaloop-integration-tests"];
    clinnet = ["clinnet"];
    linnet = [
      "linnet"
      "linnet-py"
      "linnest"
    ];
    spenso = [
      "idenso"
      "spenso"
      "spenso-hep-lib"
      "spenso-macros"
    ];
    vakint = ["vakint"];
  };
  nextestArchiveAttr = target: "checks.${system}.gammaloop-nextest-binaries-${target}";
  nextestPackageTestSupportAttr = package:
    crateTestSupportAttr (workspaceTestComponentRepresentativeFor package);
  nextestArchiveDependenciesFor = target:
    ["packages.${system}.cargoArtifacts"]
    ++ unique (map nextestPackageTestSupportAttr nextestPackageGroups.${target})
    ++ (
      if target == "python-api"
      then ["packages.${system}.gammaloop-python-module"]
      else []
    );
  # The Hakari workspace-hack deps artifact is the root for the
  # Symbolica-containing cache DAG. Higher-level crate cache jobs reach it
  # through their Guppy-resolved workspace cache dependencies.
  nextestBinaryChecks = map nextestArchiveAttr (builtins.attrNames nextestPackageGroups);
  dependencies = mergeDependencySets [
    workspaceCratePackageCacheArtifactDependencies
    workspaceCratePackageCacheDependencies
    workspaceCratePackageDependencies
    workspaceTestSupportDependencies
    {
      "packages.${system}.gammaloop" = [
        gammaloopApiPackageArtifactsAttr
        "checks.${system}.gammaloop-fmt"
      ];
      "checks.${system}.gammaloop" = ["packages.${system}.gammaloop"];
      "packages.${system}.default" = ["packages.${system}.gammaloop"];
      "packages.${system}.gammaloop-python-module" =
        workspaceCratePackageDependencies.${cratePackageAttr "gammaloop-api"} or [];
      ${gammaloopApiPackageArtifactsAttr} = [(cratePackageAttr "gammaloop-api")];
      "packages.${system}.cargoArtifacts" = [workspaceHackCacheAttr];
      "checks.${system}.gammaloop-check" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-clippy" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-doc" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-doctest" = ["packages.${system}.cargoArtifacts"];
      "packages.${system}.workspaceBuildArtifacts" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-nextest-binaries-core" = nextestArchiveDependenciesFor "core";
      "checks.${system}.gammaloop-nextest-binaries-clinnet" = nextestArchiveDependenciesFor "clinnet";
      "checks.${system}.gammaloop-nextest-binaries-integration" = nextestArchiveDependenciesFor "integration";
      "checks.${system}.gammaloop-nextest-binaries-python-api" = nextestArchiveDependenciesFor "python-api";
      "checks.${system}.gammaloop-nextest-binaries-linnet" = nextestArchiveDependenciesFor "linnet";
      "checks.${system}.gammaloop-nextest-binaries-spenso" = nextestArchiveDependenciesFor "spenso";
      "checks.${system}.gammaloop-nextest-binaries-vakint" = nextestArchiveDependenciesFor "vakint";
      "checks.${system}.gammaloop-nextest-binaries" = nextestBinaryChecks;
      "packages.${system}.linnest-wasm" = ["packages.${system}.linnestWasmCargoArtifacts"];
      "checks.${system}.linnest-wasm" = ["packages.${system}.linnest-wasm"];
      "packages.${system}.gammaloop-llvm-coverage" = ["packages.${system}.gammaloop"];
      "packages.${system}.nix-ci-check-gammaloop-doctest" = ["packages.${system}.cargoArtifacts"];
      "packages.${system}.nix-ci-check-gammaloop-nextest" =
        nextestBinaryChecks
        ++ ["packages.${system}.gammaloop-python-module"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-clinnet" = ["checks.${system}.gammaloop-nextest-binaries-clinnet"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-core" = ["checks.${system}.gammaloop-nextest-binaries-core"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-integration" = [
        "checks.${system}.gammaloop-nextest-binaries-integration"
      ];
      "packages.${system}.nix-ci-check-gammaloop-nextest-python-api" = [
        "checks.${system}.gammaloop-nextest-binaries-python-api"
        "packages.${system}.gammaloop-python-module"
      ];
      "packages.${system}.nix-ci-check-gammaloop-nextest-linnet" = ["checks.${system}.gammaloop-nextest-binaries-linnet"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-spenso" = ["checks.${system}.gammaloop-nextest-binaries-spenso"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-vakint" = ["checks.${system}.gammaloop-nextest-binaries-vakint"];
    }
  ];
  missingWorkspaceCratePackageEdges =
    builtins.filter (
      edge: !(builtins.elem edge.dependency (dependencies.${edge.dependent} or []))
    )
    workspaceCratePackageDependencyEdges;
  missingWorkspaceCratePackageCacheArtifactEdges =
    builtins.filter (
      edge: !(builtins.elem edge.dependency (dependencies.${edge.dependent} or []))
    )
    workspaceCratePackageCacheArtifactDependencyEdges;
  reciprocalWorkspaceCratePackageEdges =
    builtins.filter (
      edge: builtins.elem edge.dependent (workspaceCratePackageDependencies.${edge.dependency} or [])
    )
    workspaceCratePackageDependencyEdges;
  formatDependencyEdge = edge: "${edge.dependency} -> ${edge.dependent}";
  selfDependencies =
    builtins.filter (
      attr: builtins.elem attr (dependencies.${attr} or [])
    )
    (builtins.attrNames dependencies);
  validatedDependencies =
    assert missingWorkspaceCratePackageEdges == []
    || builtins.throw "manual NixCI dependency graph is missing workspace crate package dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge missingWorkspaceCratePackageEdges)}";
    assert missingWorkspaceCratePackageCacheArtifactEdges == []
    || builtins.throw "manual NixCI dependency graph is missing workspace crate package cache dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge missingWorkspaceCratePackageCacheArtifactEdges)}";
    assert reciprocalWorkspaceCratePackageEdges == []
    || builtins.throw "manual NixCI dependency graph contains reciprocal workspace crate package dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge reciprocalWorkspaceCratePackageEdges)}";
    assert selfDependencies == []
    || builtins.throw "manual NixCI dependency graph contains self dependencies: ${builtins.concatStringsSep ", " selfDependencies}";
      dependencies;
  doNotBuild = unique (
    [
      "checks.${system}.gammaloop-doctest"
      "checks.${system}.gammaloop-nextest"
      "checks.${system}.gammaloop-nextest-binaries"
      "checks.${system}.gammaloop-nextest-clinnet"
      "checks.${system}.gammaloop-nextest-core"
      "checks.${system}.gammaloop-nextest-integration"
      "checks.${system}.gammaloop-nextest-python-api"
      "checks.${system}.gammaloop-nextest-linnet"
      "checks.${system}.gammaloop-nextest-spenso"
      "checks.${system}.gammaloop-nextest-vakint"
      "packages.${system}.default"
      "packages.${system}.crane-ci-prebuild"
      "packages.${system}.workspaceBuildArtifacts"
      "packages.${system}.gammaloop-llvm-coverage"
      "packages.${system}.nix-ci-check-gammaloop-nextest"
    ]
    ++ map crateTestBinaryAttr workspacePackages
    ++ map cratePackageDepsAttr (
      builtins.filter (package: package != workspaceHackPackage) workspacePackagesWithDependencyArtifacts
    )
    ++ map cratePackageAttr nonWorkspaceHackPackages
  );
  # NixCI only schedules jobs it actually builds, so a dependency edge that
  # references a doNotBuild job is rejected as pointing at a non-existent job.
  # The manual graph above is constructed over the full crate/artifact DAG
  # (which keeps the drift and cycle asserts meaningful); here we drop every
  # edge touching a doNotBuild job so only ordering between built jobs remains.
  doNotBuildSet = builtins.listToAttrs (map (job: {
      name = job;
      value = true;
    })
    doNotBuild);
  isBuiltJob = job: !(doNotBuildSet ? ${job});
  buildableDependencies = deps:
    builtins.listToAttrs (
      builtins.filter (entry: entry.value != []) (map (dependent: {
          name = dependent;
          value = builtins.filter isBuiltJob deps.${dependent};
        })
        (builtins.filter isBuiltJob (builtins.attrNames deps)))
    );
in {
  systems = [system];
  inherit doNotBuild;
  fail-fast = false;
  # Keep dependency discovery manual. With generated Rust outputs,
  # automatic discovery asks NixCI to compute derivation paths for many
  # package/check attrs during `show`, including attrs listed in doNotBuild.
  # The manual graph below uses the Hakari workspace-hack cache artifact as the
  # root for Symbolica-containing cache jobs and orders nextest archive jobs
  # after the shared test-support artifacts that the archives reuse. The graph
  # is constructed over the full crate/artifact DAG so the drift and cycle
  # asserts stay meaningful, then buildableDependencies drops every edge that
  # references a doNotBuild job, since NixCI rejects edges to jobs it does not
  # build. Ordinary crate package attrs are not CI roots, so test-binary
  # generation can start before unrelated final package outputs.
  # See https://nix-ci.com/documentation/automatic-dependency-discovery
  # and https://nix-ci.com/documentation/manually-specified-dependencies
  dependency-discovery.enable = false;
  dependencies = buildableDependencies validatedDependencies;
  test = {
    gammaloop-doctest = {
      package = "packages.${system}.nix-ci-check-gammaloop-doctest";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-core = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-core";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-clinnet = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-clinnet";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-integration = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-integration";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-python-api = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-python-api";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-linnet = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-linnet";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-spenso = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-spenso";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest-vakint = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest-vakint";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };
  };
  deploy = {
    ci-passed = {
      package = "packages.${system}.nix-ci-passed";
      system = system;
      branches = "all";
    };
  };
}
