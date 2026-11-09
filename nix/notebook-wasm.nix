{
  pkgs,
  docsPkgs,
  craneLib,
  toolchain,
  workspace,
  workspaceRoot,
}: let
  inherit (pkgs) lib;
  hostRoot = "examples/notebooks/symbolica-host";
  target = "wasm32-unknown-emscripten";
  crossCrane = craneLib.overrideToolchain toolchain;
  python = docsPkgs.python314.withPackages (p: [p.wheel]);

  # Match cibuildwheel 4.2.0's cp314 Pyodide ABI. Both downloads are fixed
  # inputs; compilation and wheel packaging run offline in the Nix sandbox.
  sdk = pkgs.stdenv.mkDerivation {
    pname = "notebook-emscripten";
    version = "5.0.3";
    src = pkgs.fetchurl {
      url = "https://storage.googleapis.com/webassembly/emscripten-releases-builds/linux/065bfade0ee2c02138b83cb10049d8be95e146b6/wasm-binaries.tar.xz";
      hash = "sha256-eumW2ySRHtxTXMzwafuKPtY+jax6vljnWEMav+MlUWU=";
    };
    nativeBuildInputs = [pkgs.autoPatchelfHook pkgs.makeWrapper];
    buildInputs = [pkgs.stdenv.cc.cc.lib pkgs.zlib];
    dontBuild = true;
    dontStrip = true;
    installPhase = ''
      mkdir -p "$out"
      cp -R . "$out/"
      patchShebangs "$out/emscripten"
      cat > "$out/emscripten/.emscripten" <<EOF
      LLVM_ROOT = '$out/bin'
      BINARYEN_ROOT = '$out'
      NODE_JS = '${pkgs.nodejs_24}/bin/node'
      FROZEN_CACHE = False
      EOF
      for program in emcc em++ emar emranlib; do
        makeWrapper ${python}/bin/python "$out/bin/$program" \
          --add-flags "$out/emscripten/$program.py" \
          --set EM_CONFIG "$out/emscripten/.emscripten"
      done
    '';
    meta.platforms = ["x86_64-linux"];
  };
  crossPython =
    pkgs.runCommand "pyodide-cross-python-314.0.4" {
      src = pkgs.fetchurl {
        url = "https://github.com/pyodide/pyodide/releases/download/314.0.4/xbuildenv-314.0.4.tar.bz2";
        hash = "sha256-0+PP4JlGinRglklC3xQE7AyTc1uLp2iwoQu2YUtrtsw=";
      };
    } ''
      mkdir -p "$out"
      tar xf "$src" --strip-components=5 -C "$out" \
        xbuildenv/pyodide-root/cpython/installs/python-3.14.2
    '';
  sourcePackages = lib.unique (lib.concatMap
    (workspace.workspaceDependencyClosureFor workspace.workspaceDependencyNamesFor)
    ["feynkit-py" "linnet-py" "spynso3"]);
  notebookSource = extraFilesets:
    workspace.workspacePackageSrcForSourcePackages {
      inherit sourcePackages;
      packageSourcePackages = [];
      testSourcePackages = [];
      extraFilesets = [(workspaceRoot + "/${hostRoot}")] ++ extraFilesets;
    };
  cargoVendorDir = crossCrane.vendorCargoDeps {
    cargoLock = workspaceRoot + "/${hostRoot}/Cargo.lock";
    overrideVendorGitCheckout = packages: drv:
      if lib.any (dependency: dependency.name == "symbolica-typst-atom-payload") packages
      then
        drv.overrideAttrs (_: {
          src = builtins.fetchTree {
            type = "github";
            owner = "symbolica-dev";
            repo = "symbolica-typst-plugin";
            rev = lib.last (lib.splitString "#" (builtins.head packages).source);
          };
        })
      else drv;
  };
  cargoExtraArgs = "--locked --manifest-path ${hostRoot}/Cargo.toml --target ${target} --no-default-features --features wasm";
  common = {
    pname = "notebook-wasm";
    version = "3.0.1";
    inherit cargoVendorDir cargoExtraArgs;
    src = notebookSource [];
    nativeBuildInputs = [sdk python docsPkgs.maturin pkgs.pkg-config];
    postPatch = workspace.workspaceMissingCargoTargetsScript;
    CARGO_BUILD_TARGET = target;
    CARGO_TARGET_DIR = "target";
    CARGO_NET_OFFLINE = "true";
    CARGO_PROFILE_RELEASE_LTO = "thin";
    RUSTFLAGS = "-C link-arg=-sSIDE_MODULE=2";
    PYO3_CROSS = "1";
    PYO3_CROSS_PYTHON_VERSION = "3.14";
    PYO3_CROSS_LIB_DIR = "${crossPython}/sysconfigdata";
    PYTHON = "${python}/bin/python";
    SOURCE_DATE_EPOCH = "315532800";
    preConfigure = ''
      export EM_CACHE="$TMPDIR/emscripten-cache"
      cp -R ${sdk}/emscripten/cache "$EM_CACHE"
      chmod -R u+w "$EM_CACHE"
    '';
    doCheck = false;
    doNotLinkInheritedArtifacts = true;
    strictDeps = true;
  };
  notebookWasmDeps = crossCrane.buildDepsOnly (common
    // {
      pname = "notebook-wasm-deps";
      extraDummyScript = ''
        install -D -m 0644 ${workspaceRoot + "/${hostRoot}/Cargo.lock"} "$out/${hostRoot}/Cargo.lock"
      '';
    });
  notebookWasmHostWheel = crossCrane.mkCargoDerivation (common
    // {
      pname = "notebook-wasm-wheel";
      cargoArtifacts = notebookWasmDeps;
      # These host shims are symlinks into the owning binding packages. Python
      # source belongs to wheel packaging, not the Cargo dependency producer.
      src = notebookSource [
        (workspaceRoot + "/crates/feynkit-py/python/symbolica/community/hepkit/__init__.py")
        (workspaceRoot + "/crates/linnet-py/python/symbolica/community/graph/__init__.py")
      ];
      MATURIN_PYEMSCRIPTEN_PLATFORM_VERSION = "2026_0";
      buildPhaseCargoCommand = ''
        maturin build --release ${cargoExtraArgs} --out wheels
        ${python}/bin/python ${hostRoot}/compress_wheel.py wheels/*.whl compressed
      '';
      installPhaseCommand = ''
        mkdir -p "$out"
        cp compressed/*.whl "$out/"
      '';
    });
  # Keep the reusable cross-compilation output independent of the host-owned
  # resource bundle. Repacking updates RECORD before the final compression.
  notebookWasmWheel =
    pkgs.runCommand "notebook-wasm-wheel-with-data-3.0.1" {
      nativeBuildInputs = [python];
    } ''
      python -m wheel unpack ${notebookWasmHostWheel}/*.whl
      mkdir -p symbolica-3.0.1/symbolica/_data
      cp -R ${workspace.symbolicaRenderData} symbolica-3.0.1/symbolica/_data/fonts
      chmod -R u+w symbolica-3.0.1
      mkdir -p wheels
      python -m wheel pack symbolica-3.0.1 -d wheels
      python ${notebookSource []}/${hostRoot}/compress_wheel.py wheels/*.whl compressed
      mkdir -p "$out"
      cp compressed/*.whl "$out/"
    '';
in {
  inherit notebookWasmDeps notebookWasmWheel;
}
