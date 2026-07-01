let
  system = "x86_64-linux";
  mapAttrsToList = f: attrs: map (name: f name attrs.${name}) (builtins.attrNames attrs);
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
  workspaceDependencyName = name: dependency:
    if builtins.isAttrs dependency && dependency ? package
    then dependency.package
    else name;
  workspaceDependenciesForSections = sections: manifest:
    unique (builtins.filter (dependency: builtins.hasAttr dependency workspaceMembers) (
      builtins.concatLists (map (section: let
          dependencies = manifest.${section} or {};
        in
          map (name: workspaceDependencyName name dependencies.${name}) (builtins.attrNames dependencies))
        sections)
    ));
  workspaceCrateBuildDeps = builtins.mapAttrs (_: member:
    workspaceDependenciesForSections ["dependencies" "build-dependencies"] member.manifest)
  workspaceMembers;
  workspaceCrateTestDeps = builtins.mapAttrs (_: member:
    workspaceDependenciesForSections ["dependencies" "build-dependencies" "dev-dependencies"] member.manifest)
  workspaceMembers;
  cratePackageAttr = package: "packages.${system}.crate-${package}";
  crateTestBinaryAttr = package: "packages.${system}.crate-test-binaries-${package}";
  workspaceCratePackageDependencies =
    builtins.listToAttrs (mapAttrsToList (package: dependencies: {
        name = cratePackageAttr package;
        value = map cratePackageAttr dependencies;
      })
      workspaceCrateBuildDeps);
  workspaceCrateTestBinaryDependencies =
    builtins.listToAttrs (mapAttrsToList (package: dependencies: {
        name = crateTestBinaryAttr package;
        value = [(cratePackageAttr package)] ++ map cratePackageAttr dependencies;
      })
      workspaceCrateTestDeps);
  nextestBinaryChecks = [
    "checks.${system}.gammaloop-nextest-binaries-core"
    "checks.${system}.gammaloop-nextest-binaries-integration"
    "checks.${system}.gammaloop-nextest-binaries-linnet"
    "checks.${system}.gammaloop-nextest-binaries-spenso"
    "checks.${system}.gammaloop-nextest-binaries-vakint"
  ];
in {
  systems = [system];
  doNotBuild = [
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
  # The manual graph below keeps impure Symbolica users behind test runners
  # while still prebuilding the pure binaries those runners need.
  # See https://nix-ci.com/documentation/automatic-dependency-discovery
  # and https://nix-ci.com/documentation/manually-specified-dependencies
  dependency-discovery.enable = false;
  dependencies =
    workspaceCratePackageDependencies
    // workspaceCrateTestBinaryDependencies
    // {
      "packages.${system}.gammaloop" = [
        "checks.${system}.gammaloop-fmt"
        "devShells.${system}.default"
        "packages.${system}.crate-gammaloop-api"
      ];
      "checks.${system}.gammaloop" = ["packages.${system}.gammaloop"];
      "packages.${system}.default" = ["packages.${system}.gammaloop"];
      "packages.${system}.gammaloop-python-module" = ["packages.${system}.crate-gammaloop-api"];
      "packages.${system}.clinnet-cli" = ["packages.${system}.crate-clinnet"];
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
    };
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
