# Compile the same measurement-only source against each revision's normal API
# build. The revision owns all production code, dependencies, features and flags.
{
  revisionRoot,
  revision,
  system ? builtins.currentSystem,
}: let
  measured = builtins.getFlake "git+file://${toString revisionRoot}?rev=${revision}";
  pkgs = measured.inputs.nixpkgs.legacyPackages.${system};
  craneLib = measured.inputs.crane.mkLib pkgs;
  compiled = measured.packages.${system}.crate-gammaloop-api.overrideAttrs (old: let
    buildCommands = builtins.filter (pkgs.lib.hasPrefix "cargoWithProfile build ") (pkgs.lib.splitString "\n" old.buildPhase);
    testCommand =
      builtins.replaceStrings
      ["cargoWithProfile build " " --lib --bins"]
      ["cargoWithProfile test --example precise_sample " " -- --test-threads 1"]
      (builtins.head buildCommands);
  in
    # Crane has already rendered cargoExtraArgs into buildPhase. Preserve its
    # package/features/profile selection and add only the shared example target.
    assert builtins.length buildCommands == 1;
    assert pkgs.lib.hasSuffix " --lib --bins" (builtins.head buildCommands); {
      pname = "gammaloop-precise-performance-sample";
      PERFORMANCE_REVISION = revision;
      # Follow the repository's Nix test owner: require the supplied license,
      # with no restricted-instance fallback and no license value in logs.
      SYMBOLICA_LICENSE = builtins.getEnv "SYMBOLICA_LICENSE";
      doCheck = true;
      postPatch =
        (old.postPatch or "")
        + ''
          mkdir -p crates/gammaloop-api/examples
          cp ${./precise_sample.rs} crates/gammaloop-api/examples/precise_sample.rs
          mkdir -p tests/resources/graphs tests/resources/run_cards
            cp ${../../tests/resources/graphs/paper_figure_b1_double_triangle_soft_ir.dot} tests/resources/graphs/paper_figure_b1_double_triangle_soft_ir.dot
            cp ${../../tests/resources/run_cards/paper_figure_b1_double_triangle_soft_ir.toml} tests/resources/run_cards/paper_figure_b1_double_triangle_soft_ir.toml
        '';
      buildPhase =
        builtins.replaceStrings
        ["cargoWithProfile build "]
        ["cargoWithProfile build --example precise_sample "]
        old.buildPhase;
      checkPhase = ''
        runHook preCheck
        if [ -z "''${SYMBOLICA_LICENSE:-}" ]; then
          echo "Missing SYMBOLICA_LICENSE environment variable" >&2
          exit 1
        fi
        export RUST_MIN_STACK=67108864
        export RAYON_NUM_THREADS=1
        ${testCommand}
        runHook postCheck
      '';
    });
in
  assert builtins.match "[0-9a-f]{40}" revision != null;
    pkgs.runCommand "gammaloop-precise-performance-sample-${revision}" {
      nativeBuildInputs = [pkgs.gnutar pkgs.zstd craneLib.inheritCargoArtifactsHook];
      buildInputs = [pkgs.stdenv.cc.cc.lib pkgs.gcc.cc.lib];
    } ''
      inheritCargoArtifacts ${compiled} target
      install -D -m0755 target/ci-optim/examples/precise_sample "$out/bin/precise_sample"
    ''
