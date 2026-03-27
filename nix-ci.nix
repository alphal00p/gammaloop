let
  ciPartitionCount = 6;

  crates = [
    {
      attr = "clinnet";
      usesSymbolica = false;
    }
    {
      attr = "gammaloop-api";
      usesSymbolica = true;
    }
    {
      attr = "gammalooprs";
      usesSymbolica = true;
    }
    {
      attr = "idenso";
      usesSymbolica = true;
    }
    {
      attr = "linnest";
      usesSymbolica = false;
    }
    {
      attr = "linnet";
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
      attr = "spenso-hep-lib";
      usesSymbolica = true;
    }
    {
      attr = "spenso-macros";
      usesSymbolica = true;
    }
    {
      attr = "spynso3";
      usesSymbolica = true;
    }
    {
      attr = "vakint";
      usesSymbolica = true;
    }
  ];

  impureChecks =
    [
      "checks.x86_64-linux.gammaloop-nextest"
    ]
    ++ map
    (partition: "checks.x86_64-linux.gammaloop-nextest-partition-${toString (partition + 1)}")
    (builtins.genList (partition: partition) ciPartitionCount)
    ++ map
    (crate: "checks.x86_64-linux.crate-${crate.attr}")
    (builtins.filter (crate: crate.usesSymbolica) crates);

  symbolicaTests = builtins.listToAttrs (map
    (crate: {
      name = "crate-${crate.attr}";
      value = {
        package = "packages.x86_64-linux.nix-ci-check-${crate.attr}";
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
  test = symbolicaTests;
}
