# Deliberately independent of GammaLoop's manifests and committed CI graph.
let
  root = builtins.getEnv "GAMMALOOP_WORKSPACE_ROOT";
  lock = builtins.fromJSON (builtins.readFile (root + "/flake.lock"));
  input = name: (builtins.fetchTree lock.nodes.${lock.nodes.${lock.root}.inputs.${name}}.locked).outPath;
  pkgs = import (input "nixpkgs") {
    system = builtins.currentSystem;
    config = {};
    overlays = [];
  };
  fenix = import (input "fenix") {inherit pkgs;};
  vendorBuilder = pkgs.writeText "dyndrv-vendor.nix" ''
    let
      pkgs = import ${input "nixpkgs"} {
        system = "${pkgs.stdenv.hostPlatform.system}";
        config = {};
        overlays = [];
      };
      crane = import ${input "crane"} {inherit pkgs;};
    in crane.vendorCargoDeps {
      cargoLock = builtins.path {
        path = /. + builtins.getEnv "CARGO_DYNDRV_CARGO_LOCK";
        name = "Cargo.lock";
      };
    }
  '';
  requestedPython = builtins.getEnv "GAMMALOOP_DYNDRV_PYTHON";
  pythonName =
    if requestedPython == ""
    then "python313"
    else requestedPython;
  pythonPackage = assert pkgs.lib.assertMsg (builtins.match "python3[0-9]+" pythonName != null && builtins.hasAttr pythonName pkgs)
  "Select a standard CPython package provided by the pinned nixpkgs";
    builtins.getAttr pythonName pkgs;
in {
  tools = import ./dyndrv.nix {
    inherit pkgs fenix pythonPackage;
    vendorExpression = "import ${vendorBuilder}";
    fork = input "cargo-dyndrv";
  };
  wheel-runtime =
    pkgs.runCommand "dyndrv-wheel-runtime" {
      src = builtins.path {
        path = /. + builtins.getEnv "GAMMALOOP_DYNDRV_WHEEL_PAYLOAD";
        name = "dyndrv-wheel-payload";
      };
      # Raw path imports have no scanned references. Declare the tool/runtime
      # closure so the output scanner can retain the symlink and ELF references.
      bootstrap = builtins.storePath (builtins.getEnv "GAMMALOOP_DYNDRV_BOOTSTRAP");
      runtimeInputs = map builtins.storePath (builtins.fromJSON (builtins.getEnv "GAMMALOOP_DYNDRV_WHEEL_REFERENCES"));
    } ''
      cp -r "$src" "$out"
    '';
}
