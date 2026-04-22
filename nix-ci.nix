let
  crates = [
    {
      attr = "clinnet";
      usesSymbolica = false;
    }
    {
      attr = "linnet";
      usesSymbolica = false;
    }
    {
      attr = "spenso-macros";
      usesSymbolica = true;
    }
    {
      attr = "vakint";
      usesSymbolica = true;
    }
    {
      attr = "linnest";
      usesSymbolica = false;
    }
    {
      attr = "linnet-py";
      usesSymbolica = false;
    }
    {
      attr = "spenso";
      usesSymbolica = true;
    }
    {
      attr = "idenso";
      usesSymbolica = true;
    }
    {
      attr = "spenso-hep-lib";
      usesSymbolica = true;
    }
    {
      attr = "spynso3";
      usesSymbolica = true;
    }
    {
      attr = "gammalooprs";
      usesSymbolica = true;
    }
    {
      attr = "gammaloop-api";
      usesSymbolica = true;
    }
  ];

  impureChecks =
    (map
      (crate: "checks.x86_64-linux.clippy-${crate.attr}")
      (builtins.filter (crate: crate.usesSymbolica) crates))
    ++ map
    (crate: "checks.x86_64-linux.nextest-${crate.attr}")
    (builtins.filter (crate: crate.usesSymbolica) crates)
    ++ [
      "checks.x86_64-linux.nextest-integration"
    ];

  symbolicaClippyChecks = builtins.listToAttrs (map
    (crate: {
      name = "clippy-${crate.attr}";
      value = {
        package = "packages.x86_64-linux.nix-ci-check-clippy-${crate.attr}";
        system = "x86_64-linux";
        in-repo = true;
        secrets = ["SYMBOLICA_LICENSE"];
      };
    })
    (builtins.filter (crate: crate.usesSymbolica) crates));

  symbolicaNextestChecks = builtins.listToAttrs (map
    (crate: {
      name = "nextest-${crate.attr}";
      value = {
        package = "packages.x86_64-linux.nix-ci-check-nextest-${crate.attr}";
        system = "x86_64-linux";
        in-repo = true;
        secrets = ["SYMBOLICA_LICENSE"];
      };
    })
    (builtins.filter (crate: crate.usesSymbolica) crates));
in {
  systems = ["x86_64-linux"];
  doNotBuild = impureChecks;
  fail-fast = false;
  test =
    symbolicaClippyChecks
    // symbolicaNextestChecks
    // {
      nextest-integration = {
        package = "packages.x86_64-linux.nix-ci-check-nextest-integration";
        system = "x86_64-linux";
        in-repo = true;
        secrets = ["SYMBOLICA_LICENSE"];
      };
    };
}
