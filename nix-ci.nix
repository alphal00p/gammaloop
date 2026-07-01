let
  system = "x86_64-linux";
  unique = values:
    builtins.attrNames (builtins.listToAttrs (map (value: {
        name = value;
        value = true;
      })
      values));
  workspaceMemberDirs = let
    crateEntries = builtins.readDir ./crates;
    crateMemberDirs = map (name: "crates/${name}") (
      builtins.filter (
        name:
          crateEntries.${name}
          == "directory"
          && builtins.pathExists (./. + "/crates/${name}/Cargo.toml")
      ) (builtins.attrNames crateEntries)
    );
  in
    crateMemberDirs ++ ["tests"];
  workspaceMembers = builtins.listToAttrs (map (member: let
      manifest = builtins.fromTOML (builtins.readFile (./. + "/${member}/Cargo.toml"));
    in {
      name = manifest.package.name;
      value = {
        inherit manifest member;
      };
    })
    workspaceMemberDirs);
  cratePackageAttr = package: "packages.${system}.crate-${package}";
  crateTestBinaryAttr = package: "packages.${system}.crate-test-binaries-${package}";
  crate2nixCiPrebuildAttr = "packages.${system}.crate2nix-ci-prebuild";
  workspaceCratePackageAttrs = map cratePackageAttr (builtins.attrNames workspaceMembers);
  workspaceCrateTestBinaryAttrs = map crateTestBinaryAttr (builtins.attrNames workspaceMembers);
  mergeDependencySets = sets: let
    attrs = unique (builtins.concatLists (map builtins.attrNames sets));
  in
    builtins.listToAttrs (map (attr: {
        name = attr;
        value = unique (builtins.concatLists (map (set: set.${attr} or []) sets));
      })
      attrs);
  workspaceDependencySections = manifest: [
    (manifest.dependencies or {})
    (manifest."build-dependencies" or {})
  ];
  workspaceDependencyPackageName = name: spec:
    if builtins.isAttrs spec && spec ? package
    then spec.package
    else name;
  workspaceDependencyNamesInSection = section:
    builtins.filter (package: builtins.hasAttr package workspaceMembers) (
      map (name: workspaceDependencyPackageName name section.${name})
      (builtins.attrNames section)
    );
  workspaceDependencyNamesFor = manifest:
    unique (builtins.concatLists (map workspaceDependencyNamesInSection (workspaceDependencySections manifest)));
  workspaceCrateTestBinaryDependencies = builtins.listToAttrs (
    builtins.filter (entry: entry.value != []) (map (package: {
        name = crateTestBinaryAttr package;
        value = map crateTestBinaryAttr (workspaceDependencyNamesFor workspaceMembers.${package}.manifest);
      })
      (builtins.attrNames workspaceMembers))
  );
  workspaceCrateTestBinaryDependencyEdges = builtins.concatLists (map (dependent:
      map (dependency: {
        inherit dependency dependent;
      })
      (workspaceCrateTestBinaryDependencies.${dependent} or []))
    (builtins.attrNames workspaceCrateTestBinaryDependencies));
  crate2nixCiPrebuildDependentAttrs =
    workspaceCrateTestBinaryAttrs
    ++ [
      "packages.${system}.gammaloop"
      "packages.${system}.gammaloop-python-module"
    ];
  crate2nixCiPrebuildDependencies = builtins.listToAttrs (map (attr: {
      name = attr;
      value = [crate2nixCiPrebuildAttr];
    })
    crate2nixCiPrebuildDependentAttrs);
  nextestBinaryChecks = [
    "checks.${system}.gammaloop-nextest-binaries-core"
    "checks.${system}.gammaloop-nextest-binaries-integration"
    "checks.${system}.gammaloop-nextest-binaries-linnet"
    "checks.${system}.gammaloop-nextest-binaries-spenso"
    "checks.${system}.gammaloop-nextest-binaries-vakint"
  ];
  dependencies = mergeDependencySets [
    crate2nixCiPrebuildDependencies
    workspaceCrateTestBinaryDependencies
    {
      "packages.${system}.gammaloop" = [
        "checks.${system}.gammaloop-fmt"
        "devShells.${system}.default"
        crate2nixCiPrebuildAttr
      ];
      "checks.${system}.gammaloop" = ["packages.${system}.gammaloop"];
      "packages.${system}.default" = ["packages.${system}.gammaloop"];
      "packages.${system}.gammaloop-python-module" = [crate2nixCiPrebuildAttr];
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
      "packages.${system}.nix-ci-check-gammaloop-nextest-linnet" = ["checks.${system}.gammaloop-nextest-binaries-linnet"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-spenso" = ["checks.${system}.gammaloop-nextest-binaries-spenso"];
      "packages.${system}.nix-ci-check-gammaloop-nextest-vakint" = ["checks.${system}.gammaloop-nextest-binaries-vakint"];
    }
  ];
  missingCrate2nixCiPrebuildEdges =
    builtins.filter (
      attr: !(builtins.elem crate2nixCiPrebuildAttr (dependencies.${attr} or []))
    )
    crate2nixCiPrebuildDependentAttrs;
  missingWorkspaceCrateTestBinaryEdges =
    builtins.filter (
      edge: !(builtins.elem edge.dependency (dependencies.${edge.dependent} or []))
    )
    workspaceCrateTestBinaryDependencyEdges;
  reciprocalWorkspaceCrateTestBinaryEdges =
    builtins.filter (
      edge: builtins.elem edge.dependent (workspaceCrateTestBinaryDependencies.${edge.dependency} or [])
    )
    workspaceCrateTestBinaryDependencyEdges;
  formatDependencyEdge = edge: "${edge.dependency} -> ${edge.dependent}";
  selfDependencies =
    builtins.filter (
      attr: builtins.elem attr (dependencies.${attr} or [])
    )
    (builtins.attrNames dependencies);
  mentionedWorkspaceCratePackageAttrs =
    unique (builtins.filter (
        attr: builtins.elem attr workspaceCratePackageAttrs
      )
      ((builtins.attrNames dependencies) ++ builtins.concatLists (builtins.attrValues dependencies)));
  validatedDependencies =
    assert missingCrate2nixCiPrebuildEdges == []
    || builtins.throw "missing ${crate2nixCiPrebuildAttr} dependency edges for: ${builtins.concatStringsSep ", " missingCrate2nixCiPrebuildEdges}";
    assert missingWorkspaceCrateTestBinaryEdges == []
    || builtins.throw "manual NixCI dependency graph is missing workspace crate-test-binary dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge missingWorkspaceCrateTestBinaryEdges)}";
    assert reciprocalWorkspaceCrateTestBinaryEdges == []
    || builtins.throw "manual NixCI dependency graph contains reciprocal workspace crate-test-binary dependency edges: ${builtins.concatStringsSep ", " (map formatDependencyEdge reciprocalWorkspaceCrateTestBinaryEdges)}";
    assert selfDependencies == []
    || builtins.throw "manual NixCI dependency graph contains self dependencies: ${builtins.concatStringsSep ", " selfDependencies}";
    assert mentionedWorkspaceCratePackageAttrs == []
    || builtins.throw "manual NixCI dependency graph still mentions standalone workspace crate packages: ${builtins.concatStringsSep ", " mentionedWorkspaceCratePackageAttrs}";
      dependencies;
in {
  systems = [system];
  doNotBuild =
    workspaceCratePackageAttrs
    ++ [
      "checks.${system}.gammaloop-doctest"
      "checks.${system}.gammaloop-nextest"
      "checks.${system}.gammaloop-nextest-binaries"
      "checks.${system}.gammaloop-nextest-core"
      "checks.${system}.gammaloop-nextest-integration"
      "checks.${system}.gammaloop-nextest-linnet"
      "checks.${system}.gammaloop-nextest-spenso"
      "checks.${system}.gammaloop-nextest-vakint"
      "packages.${system}.default"
      "packages.${system}.cargoArtifacts"
      "packages.${system}.workspaceBuildArtifacts"
      "packages.${system}.gammaloop-llvm-coverage"
      "packages.${system}.nix-ci-check-gammaloop-nextest"
    ];
  fail-fast = false;
  # Keep dependency discovery manual. With the generated crate2nix outputs,
  # automatic discovery asks NixCI to compute derivation paths for many
  # package/check attrs during `show`, including attrs listed in doNotBuild.
  # The manual graph below prebuilds the crate2nix Symbolica closure before
  # crate jobs fan out, orders crate test-binary jobs by regular workspace
  # dependencies, and keeps impure Symbolica users behind test runners.
  # Test-only workspace dependencies are deliberately omitted: the flake has
  # synthetic crate2nix package IDs to break the linnet/linnest and
  # spenso-macros/spenso dev-dependency cycles, but those synthetic crates are
  # not exposed as NixCI jobs.
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
