let system = "x86_64-linux";
in {
  systems = [ system ];
  doNotBuild = [
    "checks.${system}.gammaloop-doctest"
    "checks.${system}.gammaloop-nextest"
  ];
  fail-fast = false;
  # We specify dependencies manually
  # See https://nix-ci.com/documentation/automatic-dependency-discovery
  # and https://nix-ci.com/documentation/manually-specified-dependencies
  dependency-discovery.enable = false;
  dependencies = {
    "packages.${system}.gammaloop" = [
      "packages.${system}.cargoArtifacts"
      "packages.${system}.gammaloop-fmt"
      "devShells.${system}.default"
    ];
    "packages.${system}.gammaloop-clippy" = [ "packages.${system}.gammaloop" ];
    "packages.${system}.gammaloop-doc" = [ "packages.${system}.gammaloop" ];
    "packages.${system}.gammaloop-llvm-coverage" = [ "packages.${system}.gammaloop" ];
    "packages.${system}.nix-ci-check-gammaloop-doctest" = [ "packages.${system}.gammaloop" ];
    "packages.${system}.nix-ci-check-gammaloop-nextest" = [ "packages.${system}.gammaloop" ];
  };
  test = {
    gammaloop-doctest = {
      package = "packages.x86_64-linux.nix-ci-check-gammaloop-doctest";
      system = "x86_64-linux";
      in-repo = true;
      secrets = [ "SYMBOLICA_LICENSE" ];
    };

    gammaloop-nextest = {
      package = "packages.x86_64-linux.nix-ci-check-gammaloop-nextest";
      system = "x86_64-linux";
      in-repo = true;
      secrets = [ "SYMBOLICA_LICENSE" ];
    };
  };
}
