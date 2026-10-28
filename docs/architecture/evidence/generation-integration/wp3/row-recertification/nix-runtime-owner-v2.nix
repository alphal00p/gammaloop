let
  root = /shared/lcnbr/dev/gammaloop/generation-integration-plan;
  graph = builtins.fromJSON (builtins.readFile (root + "/nix/ci-workspace-graph.json"));
  owner = import (root + "/nix/ci.nix") { workspaceGraph = graph; system = "x86_64-linux"; };
  group = builtins.head (builtins.filter (g: g.name == "feynkit") owner.groups);
in {
  inherit (group) name packages;
  runFilter = group.runFilter or null;
  generatorExtraFeatures = owner.testFeatures.feynkit-generator or [];
}
