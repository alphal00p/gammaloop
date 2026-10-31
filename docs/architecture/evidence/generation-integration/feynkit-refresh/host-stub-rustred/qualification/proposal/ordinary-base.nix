let
  f = builtins.getFlake "git+file:///shared/lcnbr/dev/gammaloop/main?rev=2cfb52b28d99f92305a9d571bbe3911cf04b8a19&shallow=1&allRefs=1";
  pkgs = f.inputs.nixpkgs.legacyPackages.x86_64-linux;
  lib = pkgs.lib;
  p = f.packages.x86_64-linux;
  base = p.crate-feynkit-py;
  root = toString f.outPath;
  normal = toString base.src;
  extraRoots = [
    "examples/notebooks/symbolica-host"
    "crates/feynkit-py/python"
    "crates/linnet-py/python"
    "crates/symbolica-data"
    "crates/alphal00p-docs-python-exporter"
    "docs/api/python/feynkit-community.pyi"
    "README.md" "LICENSE" "LICENSE-APACHE" "LICENSE-MIT"
  ];
  selectedSource = builtins.path {
    name = "source";
    path = f.outPath;
    filter = path: kind:
      let rel = lib.removePrefix (root + "/") (toString path);
      in kind == "directory"
        || builtins.pathExists (normal + "/" + rel)
        || builtins.any (prefix: rel == prefix || lib.hasPrefix (prefix + "/") rel) extraRoots;
  };
  common = old: {
    src = selectedSource;
    doCheck = false;
    doInstallCargoArtifacts = false;
    doNotLinkInheritedArtifacts = true;
    preBuild = "";
    postBuild = "";
    postInstall = "";
    CARGO_BUILD_INCREMENTAL = "false";
    CARGO_TARGET_DIR = "target";
    PYO3_PYTHON = "${pkgs.python313}/bin/python3";
  };
  hostWheel = base.overrideAttrs (old: common old // {
    name = "generation-refresh-installed-host-wheel-0.1.0";
    pname = "generation-refresh-installed-host-wheel";
    cargoArtifacts = base.cargoArtifacts;
    nativeBuildInputs = old.nativeBuildInputs ++ [pkgs.maturin];
    buildPhase = ''
      runHook preBuild
      maturin build --offline --locked --profile ci-optim \
        --manifest-path examples/notebooks/symbolica-host/Cargo.toml \
        --interpreter ${pkgs.python313}/bin/python3 \
        --out wheels --compatibility linux --auditwheel skip
      runHook postBuild
    '';
    installPhase = ''
      runHook preInstall
      mkdir -p "$out/dist" "$out/lib"
      wheels=(wheels/*.whl)
      test "''${#wheels[@]}" -eq 1
      test -s "''${wheels[0]}"
      cp "''${wheels[0]}" "$out/dist/"
      ${pkgs.python313}/bin/python3 -B - "$out" <<'EXTRACT'
      import pathlib, sys, zipfile
      output = pathlib.Path(sys.argv[1])
      wheel = next((output / "dist").glob("*.whl"))
      count = 0
      with zipfile.ZipFile(wheel) as archive:
          for name in archive.namelist():
              if name.endswith((".so", ".dylib")):
                  relative = pathlib.PurePosixPath(name)
                  assert not relative.is_absolute() and ".." not in relative.parts
                  target = output / "lib" / relative
                  target.parent.mkdir(parents=True, exist_ok=True)
                  target.write_bytes(archive.read(name))
                  count += 1
      assert count > 0
      EXTRACT
      runHook postInstall
    '';
  });
  stubExporter = base.overrideAttrs (old: common old // {
    name = "generation-refresh-feynkit-stub-exporter-0.1.0";
    pname = "generation-refresh-feynkit-stub-exporter";
    cargoArtifacts = p.crate-test-binaries-feynkit-py.artifact;
    buildPhase = ''
      runHook preBuild
      cargo build --offline --locked --profile ci-optim \
        -p alphal00p-docs-python-exporter --features feynkit
      runHook postBuild
    '';
    installPhase = ''
      mkdir -p "$out/bin"
      install -m 0755 target/ci-optim/alphal00p-docs-python-exporter "$out/bin/"
    '';
  });
in { inherit hostWheel stubExporter selectedSource; }
