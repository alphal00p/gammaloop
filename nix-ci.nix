let
  system = "x86_64-linux";
in {
  systems = [system];
  doNotBuild = [
    "checks.${system}.gammaloop-doctest"
    "checks.${system}.gammaloop-nextest"
    "checks.${system}.gammaloop-nextest-core"
    "checks.${system}.gammaloop-nextest-integration"
    "checks.${system}.gammaloop-nextest-linnet"
    "checks.${system}.gammaloop-nextest-spenso"
    "checks.${system}.gammaloop-nextest-vakint"
    "packages.${system}.default"
    "packages.${system}.gammaloop-llvm-coverage"
  ];
  fail-fast = false;
  # We specify dependencies manually
  # See https://nix-ci.com/documentation/automatic-dependency-discovery
  # and https://nix-ci.com/documentation/manually-specified-dependencies
  dependency-discovery.enable = false;
  dependencies = {
    "packages.${system}.gammaloop" = [
      "packages.${system}.cargoArtifacts"
      "checks.${system}.gammaloop-fmt"
      "devShells.${system}.default"
    ];
    "checks.${system}.gammaloop" = ["packages.${system}.gammaloop"];
    "packages.${system}.default" = ["packages.${system}.gammaloop"];
    "checks.${system}.gammaloop-clippy" = ["packages.${system}.cargoArtifacts"];
    "checks.${system}.gammaloop-doc" = ["packages.${system}.cargoArtifacts"];
    "packages.${system}.linnest-wasm" = ["packages.${system}.linnestWasmCargoArtifacts"];
    "checks.${system}.linnest-wasm" = ["packages.${system}.linnest-wasm"];
    "packages.${system}.gammaloop-llvm-coverage" = ["packages.${system}.gammaloop"];
    "packages.${system}.nix-ci-check-gammaloop-doctest" = ["packages.${system}.cargoArtifacts"];
    "packages.${system}.nix-ci-check-gammaloop-nextest" = ["packages.${system}.cargoArtifacts"];
    "packages.${system}.nix-ci-check-gammaloop-nextest-core" = ["packages.${system}.cargoArtifacts"];
    "packages.${system}.nix-ci-check-gammaloop-nextest-integration" = ["packages.${system}.cargoArtifacts"];
    "packages.${system}.nix-ci-check-gammaloop-nextest-linnet" = ["packages.${system}.cargoArtifacts"];
    "packages.${system}.nix-ci-check-gammaloop-nextest-spenso" = ["packages.${system}.cargoArtifacts"];
    "packages.${system}.nix-ci-check-gammaloop-nextest-vakint" = ["packages.${system}.cargoArtifacts"];
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
