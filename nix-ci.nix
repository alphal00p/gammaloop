{
  systems = ["x86_64-linux"];
  doNotBuild = [
    "checks.x86_64-linux.gammaloop-clippy"
    "checks.x86_64-linux.gammaloop-doctest"
    "checks.x86_64-linux.gammaloop-nextest"
  ];
  fail-fast = false;
  test = {
    gammaloop-clippy = {
      package = "packages.x86_64-linux.nix-ci-check-gammaloop-clippy";
      system = "x86_64-linux";
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-doctest = {
      package = "packages.x86_64-linux.nix-ci-check-gammaloop-doctest";
      system = "x86_64-linux";
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };

    gammaloop-nextest = {
      package = "packages.x86_64-linux.nix-ci-check-gammaloop-nextest";
      system = "x86_64-linux";
      in-repo = true;
      secrets = ["SYMBOLICA_LICENSE"];
    };
  };
}
