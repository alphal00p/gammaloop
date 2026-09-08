{
  workspaceGraph,
  system ? "x86_64-linux",
}: let
  groups = [
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
  unique = values:
    builtins.attrNames (builtins.listToAttrs (map (value: {
        name = value;
        value = true;
      })
      values));
  workspacePackages = workspaceGraph.packages;
  cratePackageDepsAttr = package: "packages.${system}.crate-deps-${package}";
  cratePackageAttr = package: "packages.${system}.crate-${package}";
  crateTestDependencyAttr = representative: "packages.${system}.crate-test-dependencies-${representative}";
  crateTestBinaryAttr = package: "packages.${system}.crate-test-binaries-${package}";
  nextestContextualTestDependencyAttr = target: package: "packages.${system}.crate-test-dependencies-${target}-${package}";
  nextestContextualTestBinaryAttr = target: package: "packages.${system}.crate-test-binaries-${target}-${package}";
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
  workspaceTestDependencyComponentRepresentatives =
    builtins.filter (representative: representative != workspaceHackPackage) workspaceTestComponentRepresentatives;
  workspaceTestComponentMembers = builtins.listToAttrs (map (representative: {
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
          (
            if package == workspaceHackPackage
            then []
            else ["packages.${system}.cargoArtifacts"]
          )
          ++ map cratePackageDepsAttr (workspaceDependencyNamesFor package)
          ++ (
            if package != workspaceHackPackage && (workspaceDependencyNamesFor package) == [] && builtins.elem package workspaceGraph.symbolica_normal_packages
            then [workspaceHackCacheAttr]
            else []
          );
      })
      workspacePackagesWithDependencyArtifacts)
  );
  workspaceCratePackageCacheArtifactDependencyEdges = builtins.concatLists (map (dependent:
    map (dependency: {
      inherit dependency dependent;
    })
    (workspaceCratePackageCacheArtifactDependencies.${dependent} or []))
  (builtins.attrNames workspaceCratePackageCacheArtifactDependencies));
  workspaceTestDependencyArtifactDependencies = builtins.listToAttrs (map (representative: {
      name = crateTestDependencyAttr representative;
      value =
        [
          "packages.${system}.cargoArtifacts"
          workspaceHackCacheAttr
        ]
        ++ map crateTestDependencyAttr (
          builtins.filter (
            dependencyRepresentative: dependencyRepresentative != workspaceHackPackage
          )
          (workspaceTestComponentDependencyRepresentativesFor representative)
        );
    })
    workspaceTestDependencyComponentRepresentatives);
  workspaceTestBinaryArtifactDependencies = builtins.listToAttrs (map (package: {
      name = crateTestBinaryAttr package;
      value = [(crateTestDependencyAttr (workspaceTestComponentRepresentativeFor package))];
    })
    nonWorkspaceHackPackages);
  nextestPackageGroups = builtins.listToAttrs (map (group: {
      name = group.name;
      value = group.packages;
    })
    groups);
  nextestArchiveAttr = target: "checks.${system}.gammaloop-nextest-binaries-${target}";
  nextestPackageArtifactAttrFor = target: package:
    if target == "python-api"
    then nextestContextualTestBinaryAttr target package
    else crateTestBinaryAttr package;
  nextestArchiveDependenciesFor = target:
    ["packages.${system}.cargoArtifacts"]
    ++ unique (map (nextestPackageArtifactAttrFor target) nextestPackageGroups.${target})
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
    workspaceTestDependencyArtifactDependencies
    workspaceTestBinaryArtifactDependencies
    (builtins.listToAttrs (map (group: {
        name = nextestArchiveAttr group.name;
        value = nextestArchiveDependenciesFor group.name;
      })
      groups))
    (builtins.listToAttrs (map (group: {
        name = "packages.${system}.nix-ci-check-gammaloop-nextest-${group.name}";
        value = ["packages.${system}.ci-test-inputs"];
      })
      groups))
    {
      "packages.${system}.ci-test-inputs" = nextestBinaryChecks ++ ["packages.${system}.gammaloop-python-module"];
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
      ${nextestContextualTestDependencyAttr "python-api" "gammaloop-integration-tests"} =
        workspaceTestDependencyArtifactDependencies.${crateTestDependencyAttr (workspaceTestComponentRepresentativeFor "gammaloop-integration-tests")}
        ++ ["packages.${system}.gammaloop-python-module"];
      ${nextestContextualTestBinaryAttr "python-api" "gammaloop-integration-tests"} = [
        (nextestContextualTestDependencyAttr "python-api" "gammaloop-integration-tests")
        "packages.${system}.gammaloop-python-module"
      ];
      "checks.${system}.gammaloop-check" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-clippy" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-doc" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-doctest" = ["packages.${system}.cargoArtifacts"];
      "packages.${system}.workspaceBuildArtifacts" = ["packages.${system}.cargoArtifacts"];
      "checks.${system}.gammaloop-nextest-binaries" = nextestBinaryChecks;
      "packages.${system}.linnest-wasm" = ["packages.${system}.linnestWasmCargoArtifacts"];
      "checks.${system}.linnest-wasm" = ["packages.${system}.linnest-wasm"];
      "packages.${system}.gammaloop-llvm-coverage" = ["packages.${system}.gammaloop"];
      "packages.${system}.nix-ci-check-gammaloop-doctest" = ["packages.${system}.cargoArtifacts"];
      "packages.${system}.nix-ci-check-gammaloop-nextest" =
        nextestBinaryChecks
        ++ ["packages.${system}.gammaloop-python-module"];
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
  validatedDependencies = assert missingWorkspaceCratePackageEdges
  == []
  || builtins.throw "manual NixCI dependency graph is missing workspace crate package dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge missingWorkspaceCratePackageEdges)}";
  assert missingWorkspaceCratePackageCacheArtifactEdges
  == []
  || builtins.throw "manual NixCI dependency graph is missing workspace crate package cache dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge missingWorkspaceCratePackageCacheArtifactEdges)}";
  assert reciprocalWorkspaceCratePackageEdges
  == []
  || builtins.throw "manual NixCI dependency graph contains reciprocal workspace crate package dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge reciprocalWorkspaceCratePackageEdges)}";
  assert selfDependencies
  == []
  || builtins.throw "manual NixCI dependency graph contains self dependencies: ${builtins.concatStringsSep ", " selfDependencies}"; dependencies;
  primaryChecks =
    [
      "checks.${system}.gammaloop-clippy"
      "checks.${system}.gammaloop-fmt"
      "checks.${system}.gammaloop-guppy-workspace-graph"
      "packages.${system}.nix-ci-check-gammaloop-doctest"
      "packages.${system}.nix-ci-passed"
    ]
    ++ map (group: "packages.${system}.nix-ci-check-gammaloop-nextest-${group.name}") groups;
  onlyBuild = unique (primaryChecks
    ++ [
      "packages.${system}.ci-test-inputs"
      "packages.${system}.cargoArtifacts"
      workspaceHackCacheAttr
    ]);
  # NixCI only schedules jobs it actually builds, so a dependency edge that
  # references an unselected job is rejected as pointing at a non-existent job.
  # The manual graph above is constructed over the full crate/artifact DAG
  # (which keeps the drift and cycle asserts meaningful); here hidden paths are
  # contracted to their nearest built dependency so their ordering is retained.
  selectedJobs = builtins.listToAttrs (map (job: {
      name = job;
      value = true;
    })
    onlyBuild);
  isBuiltJob = job: selectedJobs ? ${job};
  builtDependencyFrontierFor = deps: dependent:
    unique (map (entry: entry.key) (builtins.filter (
        entry: isBuiltJob entry.key
      ) (builtins.genericClosure {
        startSet = map (dependency: {key = dependency;}) (deps.${dependent} or []);
        operator = entry:
          if isBuiltJob entry.key
          then []
          else map (dependency: {key = dependency;}) (deps.${entry.key} or []);
      })));
  buildableDependencies = deps:
    builtins.listToAttrs (
      builtins.filter (entry: entry.value != []) (map (dependent: {
          name = dependent;
          value = builtDependencyFrontierFor deps dependent;
        })
        (builtins.filter isBuiltJob (builtins.attrNames deps)))
    );
  projectedDependencies = buildableDependencies validatedDependencies;
  projectedDependencyClosureFor = dependent:
    map (entry: entry.key) (builtins.genericClosure {
      startSet = map (dependency: {key = dependency;}) (projectedDependencies.${dependent} or []);
      operator = entry:
        map (dependency: {key = dependency;}) (projectedDependencies.${entry.key} or []);
    });
  projectedDependencyCycles = builtins.filter (
    attr: builtins.elem attr (projectedDependencyClosureFor attr)
  ) (builtins.attrNames projectedDependencies);
  validatedProjectedDependencies = assert projectedDependencyCycles
  == []
  || builtins.throw "projected NixCI dependency graph contains cycles through: ${builtins.concatStringsSep ", " projectedDependencyCycles}"; projectedDependencies;
in {
  inherit groups;
  configuration = {
    systems = [system];
    inherit onlyBuild;
    fail-fast = false;
    fail-on-dangling-dependencies = true;
    # Trial synchronous discovery before jobs start. The earlier concern was
    # evaluating derivation paths for many generated package/check attrs during
    # `show`, including outputs not selected for CI; narrow onlyBuild selection
    # remains important. Keep the explicit graph and its barrier ordering too.
    # The manual graph below uses the Hakari workspace-hack cache artifact as the
    # root for Symbolica-containing cache jobs and orders nextest archive jobs
    # after the package-local test-binary artifacts that the archives reuse. The
    # exported artifact attrs are revision-scoped symlink barriers around stable
    # Cargo artifacts, so a memoized top-level result cannot release consumers
    # without one worker realizing and publishing the underlying closure. The
    # graph is constructed over the full crate/artifact DAG so the drift and
    # cycle asserts stay meaningful, then hidden paths are contracted to their
    # nearest built producer because NixCI rejects edges to jobs it does not
    # build. Ordinary crate package attrs are not CI roots, so test-binary
    # generation can start before unrelated final package outputs.
    # See https://nix-ci.com/documentation/automatic-dependency-discovery
    # and https://nix-ci.com/documentation/manually-specified-dependencies
    dependency-discovery = {
      enable = true;
      synchronous = true;
    };
    dependencies = validatedProjectedDependencies;
    test = builtins.listToAttrs (map (name: {
      inherit name;
      value = {
        package = "packages.${system}.nix-ci-check-${name}";
        inherit system;
        in-repo = true;
        secrets = ["SYMBOLICA_LICENSE"];
      };
    }) (["gammaloop-doctest"] ++ map (group: "gammaloop-nextest-${group.name}") groups));
    deploy = {
      ci-passed = {
        package = "packages.${system}.nix-ci-passed";
        system = system;
        branches = "all";
      };
    };
  };
}
