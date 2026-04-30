let
  system = "x86_64-linux";
in {
  systems = [system];
  doNotBuild = [
    "checks.${system}.gammaloop-doctest"
    "checks.${system}.gammaloop-nextest"
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
    "packages.${system}.gammaloop-llvm-coverage" = ["packages.${system}.gammaloop"];
    "packages.${system}.nix-ci-check-gammaloop-doctest" = ["packages.${system}.cargoArtifacts"];
    "packages.${system}.nix-ci-check-gammaloop-nextest" = ["packages.${system}.cargoArtifacts"];
  };
  test = {
    gammaloop-doctest = {
      package = "packages.${system}.nix-ci-check-gammaloop-doctest";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest = {
      package = "packages.${system}.nix-ci-check-gammaloop-nextest";
      system = system;
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };
  };
}
