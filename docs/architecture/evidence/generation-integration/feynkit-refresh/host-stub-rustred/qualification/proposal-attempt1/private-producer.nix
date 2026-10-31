let
  f = builtins.getFlake "git+file:///shared/lcnbr/dev/gammaloop/main?rev=2cfb52b28d99f92305a9d571bbe3911cf04b8a19&shallow=1&allRefs=1";
  pkgs = f.inputs.nixpkgs.legacyPackages.x86_64-linux;
  base = f.packages.x86_64-linux.crate-feynkit-py;
  rustredSource = builtins.path {
    path = /nix/store/5hr22rycs6dpj77afgq0maqd8qbipk6a-rustred-source-7c1ed037;
    name = "rustred-source-7c1ed037";
    sha256 = "sha256-6Sgm1lGDDGCe2ch784rGp0krsqlOiZtKAEuafvNnKlA=";
  };
  source = builtins.path {
    path = /nix/store/h5g6njarsri3hacwzgzcrh4dzxm480ni-source;
    name = "source";
    sha256 = "sha256-xVFNiWJSJzU3+3M9Lk1SHEx3mesvPOmS/2/oL3EW6p0=";
  };
  hostWheel = base.overrideAttrs (old: {
    name = "generation-refresh-rustred-validation-host-wheel-0.1.0";
    pname = "generation-refresh-rustred-validation-host-wheel";
    src = source;
    # All external package versions/sources are the already frozen root vendor set.
    cargoVendorDir = base.cargoVendorDir;
    cargoArtifacts = base.cargoArtifacts;
    RUSTRED_SOURCE = rustredSource;
    # This private qualification covers the one-denominator QCD family only.
    RUSTRED_RUNTIME_ARITIES = "1";
    doCheck = false;
    doInstallCargoArtifacts = false;
    doNotLinkInheritedArtifacts = true;
    preBuild = "";
    postBuild = "";
    postInstall = "";
    CARGO_BUILD_INCREMENTAL = "false";
    CARGO_TARGET_DIR = "target";
    PYO3_PYTHON = "${pkgs.python313}/bin/python3";
    nativeBuildInputs = old.nativeBuildInputs ++ [pkgs.maturin];
    buildPhase = ''
      runHook preBuild
      rustc -Vv > producer-rustc.txt
      cargo metadata --offline --locked --format-version 1 \
        --manifest-path examples/notebooks/symbolica-host/Cargo.toml \
        > private-host-metadata.json
      ${pkgs.python313}/bin/python3 -B - <<'METADATA'
      import json, pathlib
      metadata = json.loads(pathlib.Path("private-host-metadata.json").read_text())
      packages = metadata["packages"]
      nodes = {n["id"]: n for n in metadata["resolve"]["nodes"]}
      for name in ("symbolica", "numerica", "graphica"):
          matches = [p for p in packages if p["name"] == name]
          assert len(matches) == 1, (name, matches)
          assert matches[0]["version"] == "3.0.1"
          assert matches[0]["source"].endswith("#8eb7c3729d6abefc916decc1bdf0c06536f8af80")
      for name in ("pyo3", "pyo3-ffi", "pyo3-build-config", "pyo3-macros", "pyo3-macros-backend"):
          matches = [p for p in packages if p["name"] == name]
          assert len(matches) == 1 and matches[0]["version"] == "0.28.3", name
      source_root = pathlib.Path.cwd().resolve()
      for package in packages:
          if package["name"].startswith("feynkit"):
              assert package["source"] is None, package
              assert pathlib.Path(package["manifest_path"]).is_relative_to(source_root / "crates")
      for name in ("rustred-feynkit", "rustred", "rustred-order"):
          package, = [p for p in packages if p["name"] == name]
          assert package["source"] is None
          assert pathlib.Path(package["manifest_path"]).is_relative_to(pathlib.Path("${rustredSource}"))
          assert nodes[package["id"]]["features"] == [], (name, nodes[package["id"]])
      assert not any(p["name"] in ("rustred-app", "rustred-python") and p["id"] in nodes for p in packages)
      METADATA
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
      assert count == 1
      EXTRACT
      mkdir -p "$out/provenance/cargo-fingerprints" "$out/provenance/depfiles"
      cp private-host-metadata.json producer-rustc.txt "$out/provenance/"
      cp examples/notebooks/symbolica-host/Cargo.toml "$out/provenance/host-Cargo.toml"
      cp examples/notebooks/symbolica-host/Cargo.lock "$out/provenance/host-Cargo.lock"
      printf '%s\n' "$RUSTRED_RUNTIME_ARITIES" > "$out/provenance/runtime-arities.txt"
      for fingerprints in target/ci-optim/.fingerprint/{alphal00p-notebook-host,rustred,rustred-order,rustred-feynkit,symbolica,numerica,graphica,pyo3,pyo3-ffi,pyo3-build-config}-*; do
        if [ -d "$fingerprints" ]; then cp -a "$fingerprints" "$out/provenance/cargo-fingerprints/"; fi
      done
      for depfile in target/ci-optim/deps/{core,rustred,rustred_order,rustred_feynkit,symbolica,numerica,graphica,pyo3,pyo3_ffi}-*.d target/ci-optim/libcore.d; do
        if [ -f "$depfile" ]; then cp "$depfile" "$out/provenance/depfiles/"; fi
      done
      runHook postInstall
    '';
  });
in { inherit hostWheel source rustredSource; }
