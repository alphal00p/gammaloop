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
  workspaceHackCacheAttr = cratePackageDepsAttr workspaceHackPackage;
  workspacePackageGraphAttr = package:
    if package == "gammaloop-api"
    then "packages.${system}.gammaloop"
    else cratePackageAttr package;
  workspaceTestSupportGraphAttr = representative:
    if representative == workspaceHackPackage
    then workspaceHackCacheAttr
    else crateTestSupportAttr representative;
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
    workspacePackages);
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
  workspaceCrateTestSupportDependencies = builtins.listToAttrs (
    builtins.filter (entry: entry.value != []) (map (representative: {
        name = crateTestSupportAttr representative;
        value =
          [workspaceHackCacheAttr]
          ++ map workspaceTestSupportGraphAttr (workspaceTestComponentDependencyRepresentativesFor representative);
      })
      workspaceTestSupportComponentRepresentatives)
  );
  workspaceCrateTestSupportDependencyEdges = builtins.concatLists (map (dependent:
      map (dependency: {
        inherit dependency dependent;
      })
      (workspaceCrateTestSupportDependencies.${dependent} or []))
    (builtins.attrNames workspaceCrateTestSupportDependencies));
  workspaceCrateTestBinaryCacheDependencies = builtins.listToAttrs (map (package: {
      name = crateTestBinaryAttr package;
      value = [
        (workspaceTestSupportGraphAttr (workspaceTestComponentRepresentativeFor package))
        (workspacePackageGraphAttr package)
      ];
    })
    workspacePackages);
  # The Hakari workspace-hack deps artifact is the root for the
  # Symbolica-containing cache DAG. Higher-level crate cache jobs reach it
  # through their Guppy-resolved workspace cache dependencies.
  nextestBinaryChecks = [
    "checks.${system}.gammaloop-nextest-binaries-core"
    "checks.${system}.gammaloop-nextest-binaries-integration"
    "checks.${system}.gammaloop-nextest-binaries-python-api"
    "checks.${system}.gammaloop-nextest-binaries-linnet"
    "checks.${system}.gammaloop-nextest-binaries-spenso"
    "checks.${system}.gammaloop-nextest-binaries-vakint"
  ];
  dependencies = mergeDependencySets [
    workspaceCratePackageCacheArtifactDependencies
    workspaceCratePackageCacheDependencies
    workspaceCratePackageDependencies
    workspaceCrateTestSupportDependencies
    workspaceCrateTestBinaryCacheDependencies
    {
      "packages.${system}.gammaloop" = [
        "checks.${system}.gammaloop-fmt"
        "devShells.${system}.default"
      ];
      "checks.${system}.gammaloop" = ["packages.${system}.gammaloop"];
      "packages.${system}.default" = ["packages.${system}.gammaloop"];
      "packages.${system}.gammaloop-python-module" =
        workspaceCratePackageDependencies."packages.${system}.gammaloop" or [];
      "packages.${system}.cargoArtifacts" = [workspaceHackCacheAttr];
      "checks.${system}.gammaloop-clippy" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-doc" = ["packages.${system}.cargoArtifacts"];
      "packages.${system}.workspaceBuildArtifacts" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-nextest-binaries-core" = [
        "packages.${system}.crate-test-binaries-gammaloop-api"
        "packages.${system}.crate-test-binaries-gammaloop-tracing-filter"
        "packages.${system}.crate-test-binaries-gammaloop-tracing-filter-macros"
        "packages.${system}.crate-test-binaries-gammalooprs"
      ];
      "checks.${system}.gammaloop-nextest-binaries-integration" = ["packages.${system}.crate-test-binaries-gammaloop-integration-tests"];
      "checks.${system}.gammaloop-nextest-binaries-python-api" = ["packages.${system}.crate-test-binaries-gammaloop-integration-tests"];
      "checks.${system}.gammaloop-nextest-binaries-linnet" = [
        "packages.${system}.crate-test-binaries-clinnet"
        "packages.${system}.crate-test-binaries-linnet"
        "packages.${system}.crate-test-binaries-linnet-py"
        "packages.${system}.crate-test-binaries-linnest"
      ];
      "checks.${system}.gammaloop-nextest-binaries-spenso" = [
        "packages.${system}.crate-test-binaries-idenso"
        "packages.${system}.crate-test-binaries-spenso"
        "packages.${system}.crate-test-binaries-spenso-hep-lib"
        "packages.${system}.crate-test-binaries-spenso-macros"
        "packages.${system}.crate-test-binaries-spynso3"
      ];
      "checks.${system}.gammaloop-nextest-binaries-vakint" = ["packages.${system}.crate-test-binaries-vakint"];
      "checks.${system}.gammaloop-nextest-binaries" = nextestBinaryChecks;
      "packages.${system}.linnest-wasm" = ["packages.${system}.linnestWasmCargoArtifacts"];
      "checks.${system}.linnest-wasm" = ["packages.${system}.linnest-wasm"];
      "packages.${system}.gammaloop-llvm-coverage" = ["packages.${system}.gammaloop"];
      "packages.${system}.nix-ci-check-gammaloop-doctest" = ["packages.${system}.cargoArtifacts"];
      "packages.${system}.nix-ci-check-gammaloop-nextest" =
        nextestBinaryChecks
        ++ ["packages.${system}.gammaloop-python-module"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-core" = ["checks.${system}.gammaloop-nextest-binaries-core"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-integration" = [
        "checks.${system}.gammaloop-nextest-binaries-integration"
        "packages.${system}.gammaloop-python-module"
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
  missingWorkspaceCrateTestSupportEdges =
    builtins.filter (
      edge: !(builtins.elem edge.dependency (dependencies.${edge.dependent} or []))
    )
    workspaceCrateTestSupportDependencyEdges;
  reciprocalWorkspaceCratePackageEdges =
    builtins.filter (
      edge: builtins.elem edge.dependent (workspaceCratePackageDependencies.${edge.dependency} or [])
    )
    workspaceCratePackageDependencyEdges;
  reciprocalWorkspaceCrateTestSupportEdges =
    builtins.filter (
      edge: builtins.elem edge.dependent (workspaceCrateTestSupportDependencies.${edge.dependency} or [])
    )
    workspaceCrateTestSupportDependencyEdges;
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
    assert missingWorkspaceCrateTestSupportEdges == []
    || builtins.throw "manual NixCI dependency graph is missing workspace crate-test-support dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge missingWorkspaceCrateTestSupportEdges)}";
    assert reciprocalWorkspaceCratePackageEdges == []
    || builtins.throw "manual NixCI dependency graph contains reciprocal workspace crate package dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge reciprocalWorkspaceCratePackageEdges)}";
    assert reciprocalWorkspaceCrateTestSupportEdges == []
    || builtins.throw "manual NixCI dependency graph contains reciprocal workspace crate-test-support dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge reciprocalWorkspaceCrateTestSupportEdges)}";
    assert selfDependencies == []
    || builtins.throw "manual NixCI dependency graph contains self dependencies: ${builtins.concatStringsSep ", " selfDependencies}";
      dependencies;
in {
  systems = [system];
  doNotBuild =
    [
      # Alias of packages.${system}.gammaloop; the product attr carries the
      # gammaloop-api crate graph.
      (cratePackageAttr "gammaloop-api")
      "checks.${system}.gammaloop-doctest"
      "checks.${system}.gammaloop-nextest"
      "checks.${system}.gammaloop-nextest-binaries"
      "checks.${system}.gammaloop-nextest-core"
      "checks.${system}.gammaloop-nextest-integration"
      "checks.${system}.gammaloop-nextest-python-api"
      "checks.${system}.gammaloop-nextest-linnet"
      "checks.${system}.gammaloop-nextest-spenso"
      "checks.${system}.gammaloop-nextest-vakint"
      "packages.${system}.default"
      "packages.${system}.crane-ci-prebuild"
      "packages.${system}.cargoArtifacts"
      "packages.${system}.workspaceBuildArtifacts"
      "packages.${system}.gammaloop-llvm-coverage"
      "packages.${system}.nix-ci-check-gammaloop-nextest"
    ];
  fail-fast = false;
  # Keep dependency discovery manual. With generated Rust outputs,
  # automatic discovery asks NixCI to compute derivation paths for many
  # package/check attrs during `show`, including attrs listed in doNotBuild.
  # The manual graph below uses the Hakari workspace-hack cache artifact as the
  # root for Symbolica-containing cache jobs, orders crate package jobs by
  # Guppy's normal resolved package closures, and orders test-binary cache jobs
  # by SCC-compressed Guppy test dependencies. Packages in a test-dependency SCC
  # share one crate-test-support derivation, so dev-dependency cycles like
  # linnet/linnest and spenso-macros/spenso do not become reciprocal NixCI
  # edges.
  # See https://nix-ci.com/documentation/automatic-dependency-discovery
  # and https://nix-ci.com/documentation/manually-specified-dependencies
  dependency-discovery.enable = false;
  dependencies = validatedDependencies;
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
