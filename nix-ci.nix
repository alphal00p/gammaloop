let
  system = "x86_64-linux";
in {
  systems = [system];
  doNotBuild = [
    "checks.${system}.gammaloop-doctest"
    "checks.${system}.gammaloop-nextest"
    "checks.${system}.gammaloop-nextest-binaries"
    "checks.${system}.gammaloop-nextest-clinnet"
    "checks.${system}.gammaloop-nextest-core"
    "checks.${system}.gammaloop-nextest-integration"
    "checks.${system}.gammaloop-nextest-linnet"
    "checks.${system}.gammaloop-nextest-python-api"
    "checks.${system}.gammaloop-nextest-spenso"
    "checks.${system}.gammaloop-nextest-vakint"
    "packages.${system}.crane-ci-prebuild"
    "packages.${system}.default"
    "packages.${system}.gammaloop-llvm-coverage"
    "packages.${system}.nix-ci-check-gammaloop-nextest"
    "packages.${system}.workspaceBuildArtifacts"
  ];
  fail-fast = false;
  # Automatic discovery derives the job graph from the derivations themselves,
  # so a dependency nobody wrote down here is still ordered correctly, and
  # synchronous means no job starts before the graph is known: nothing gets
  # scheduled that discovery would have ruled out as already building
  # elsewhere.
  # See https://nix-ci.com/documentation/automatic-dependency-discovery
  dependency-discovery = {
    enable = true;
    synchronous = true;
  };
  # Only the orderings discovery cannot see, because they are not derivation
  # dependencies at all. The impure test runners are small scripts that run
  # `nix build` themselves, so nothing in their derivations mentions the
  # archives and artifacts they will ask for, and without these edges each of
  # them would build its own. Formatting gates the package rather than
  # feeding it.
  # See https://nix-ci.com/documentation/manually-specified-dependencies
  dependencies = {
    "packages.${system}.gammaloop" = ["checks.${system}.gammaloop-fmt"];
    "packages.${system}.nix-ci-check-gammaloop-doctest" = ["packages.${system}.cargoArtifacts"];
    "packages.${system}.nix-ci-check-gammaloop-nextest-clinnet" = ["checks.${system}.gammaloop-nextest-binaries-clinnet"];
    "packages.${system}.nix-ci-check-gammaloop-nextest-core" = ["checks.${system}.gammaloop-nextest-binaries-core"];
    "packages.${system}.nix-ci-check-gammaloop-nextest-integration" = ["checks.${system}.gammaloop-nextest-binaries-integration"];
    "packages.${system}.nix-ci-check-gammaloop-nextest-linnet" = ["checks.${system}.gammaloop-nextest-binaries-linnet"];
    "packages.${system}.nix-ci-check-gammaloop-nextest-python-api" = [
      "checks.${system}.gammaloop-nextest-binaries-python-api"
      "packages.${system}.gammaloop-python-module"
    ];
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
