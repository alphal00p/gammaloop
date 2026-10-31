{
  pkgs,
  fenix,
  fork,
  vendorExpression,
  pythonPackage ? pkgs.python313,
}: let
  inherit (pkgs) lib;
  forkOutputs = (import (fork + "/flake.nix")).outputs {
    self = {};
    nixpkgs = {inherit lib;};
  };
  forkPackages = pkgs.extend forkOutputs.overlays.default;
  driver = forkPackages.cargo-dyndrv;
  worker = forkPackages.cargo-dyndrv-worker;
  toolchain = fenix.stable.withComponents [
    "cargo"
    "clippy"
    "rustc"
    "rust-std"
  ];
  workerPython = pkgs.python313.withPackages (ps: [ps.tomli-w]);
  # These programs are available to native build scripts and the compiler. Keep
  # session orchestration out of their PATH and immutable input closure.
  compilerPrograms = [
    toolchain
    pythonPackage
    pkgs.bash
    pkgs.coreutils
    pkgs.findutils
    pkgs.gnutar
    pkgs.pkg-config
    pkgs.gcc
    pkgs.gnum4
    pkgs.gnumake
    pkgs.gnused
    pkgs.gnugrep
    pkgs.gawk
    pkgs.diffutils
    pkgs.patch
    pkgs.gzip
    pkgs.bzip2
    pkgs.xz
    pkgs.perl
    pkgs.form
  ];
  runtimePrograms = [
    worker
    driver
    (lib.hiPrio workerPython)
    fenix.stable.rustfmt
    pkgs.nix
    pkgs.git
    pkgs.jujutsu
    pkgs.rsync
    pkgs.zstd
    pkgs.util-linux
    pkgs.cargo-nextest
    pkgs.maturin
    pkgs.uv
    pkgs.patchelf
  ];
  libraries = [
    pkgs.openssl
    pkgs.gmp
    pkgs.mpfr
    pkgs.libmpc
    pythonPackage
  ];
  tools = pkgs.buildEnv {
    name = "gammaloop-dyndrv-tools";
    # The worker needs tomli-w even when PyO3 selects another interpreter.
    paths = runtimePrograms ++ compilerPrograms;
  };
  compilerInputs = lib.unique (compilerPrograms ++ libraries ++ map lib.getDev libraries ++ [pkgs.stdenv]);
  compilerEnvironment =
    pkgs.runCommand "gammaloop-dyndrv-compiler-environment.json"
    {
      nativeBuildInputs = compilerPrograms;
      buildInputs = libraries;
      RUSTC = "${toolchain}/bin/rustc";
      RUSTDOC = "${toolchain}/bin/rustdoc";
      CC = "${pkgs.stdenv.cc}/bin/cc";
      CXX = "${pkgs.stdenv.cc}/bin/c++";
      RUSTFLAGS = "-C linker=${pkgs.stdenv.cc}/bin/cc";
      PYO3_PYTHON = "${pythonPackage}/bin/python3";
      PYTHON = "${pythonPackage}/bin/python3";
      LD_LIBRARY_PATH = lib.makeLibraryPath (libraries ++ [pkgs.stdenv.cc.cc.lib]);
      # This output is metadata, not a native library search directory.
      NIX_NO_SELF_RPATH = 1;
      NIX_OUTPATH_USED_AS_RANDOM_SEED = "gammaloop-dyndrv-v1";
    }
    ''
      export PATH=${lib.makeBinPath compilerPrograms}
      ${pythonPackage}/bin/python3 - "$out" <<'PY'
      import json, os, pathlib, sys
      names = {"PATH", "CC", "CXX", "RUSTFLAGS", "RUSTC", "RUSTDOC", "PYO3_PYTHON", "PYTHON",
               "LD_LIBRARY_PATH", "PKG_CONFIG_PATH", "PKG_CONFIG_FOR_BUILD", "CONFIG_SHELL", "SHELL"}
      environment = {name: value for name, value in os.environ.items()
                     if name in names or (name.startswith("NIX_") and name not in {"NIX_BUILD_TOP", "NIX_LOG_FD"})}
      pathlib.Path(sys.argv[1]).write_text(json.dumps({
        "environment": environment,
        "inputs": ${builtins.toJSON (map toString compilerInputs)},
      }, sort_keys=True) + "\n")
      PY
    '';
  compilerClosure = pkgs.closureInfo {
    rootPaths = compilerInputs ++ [compilerEnvironment];
  };
  # Runtime metadata is assembled separately so session-tool updates cannot
  # change the compiler environment or bring their references into its closure.
  environment = pkgs.runCommand "gammaloop-dyndrv-environment.json" {} ''
    ${workerPython}/bin/python3 - "$out" <<'PY'
    import json, pathlib, sys
    spec = json.loads(pathlib.Path("${compilerEnvironment}").read_text())
    spec.update({
      "compiler_closure": "${compilerClosure}",
      "runtime_path": "${tools}/bin",
      "nix": "${pkgs.nix}/bin/nix",
      "cargo": "${toolchain}/bin/cargo",
      "rustdoc": "${toolchain}/bin/rustdoc",
      "driver": "${driver}/bin/cargo-dyndrv",
      "worker": "${worker}/bin/cargo-dyndrv-worker",
      "worker_script": "${worker}/share/cargo-dyndrv/worker.py",
      "vendor_expression": ${builtins.toJSON vendorExpression},
      "vendor_lockfile_only": True,
      "runtime_environment_names": ["SYMBOLICA_LICENSE", "SYMBOLICA_LICENSE_PATH", "SYMBOLICA_LICENSE_SIGNED", "SYMBOLICA_LICENSE_FILE", "FEYNKIT_ORACLE_DUMP"],
      "runtime_environment_prefixes": ["GL_", "GAMMALOOP_"],
      "private_environment_names": ["SYMBOLICA_LICENSE", "SYMBOLICA_LICENSE_PATH", "SYMBOLICA_LICENSE_SIGNED", "SYMBOLICA_LICENSE_FILE"],
      "forward_environment_names": ["NO_SYMBOLICA_OEM_LICENSE"],
      "workspace_environment_names": ["GAMMALOOP_WORKSPACE_ROOT", "NEXTEST_WORKSPACE_ROOT", "INSTA_WORKSPACE_ROOT"],
      "action_executables": {"python": "${pkgs.maturin}/bin/maturin"},
      "post_commands": {"python": [[
        "${workerPython}/bin/python3", "-c",
        "import json,sys; from pathlib import Path; sys.path.insert(0,sys.argv[1]); from python_wheel import wheel_render_data; args=sys.argv[3:]; wheel_render_data(json.loads(Path(sys.argv[2]).read_text()),args[args.index('--manifest-path')+1],Path(args[args.index('--out')+1])/'render-data')",
        "{workspace}/scripts", "{bootstrap}/environment.json",
      ]]},
      "maturin": "${pkgs.maturin}/bin/maturin",
      "nextest": "${pkgs.cargo-nextest}/bin/cargo-nextest",
      "python": "${workerPython}/bin/python3",
    })
    pathlib.Path(sys.argv[1]).write_text(json.dumps(spec, sort_keys=True) + "\n")
    PY
  '';
  closure = pkgs.closureInfo {
    rootPaths =
      compilerInputs
      ++ runtimePrograms
      ++ [
        compilerClosure
        compilerEnvironment
        tools
        environment
      ];
  };
in
  pkgs.linkFarm "gammaloop-dyndrv-bootstrap" [
    {
      name = "worker-tests";
      path = worker.testsSource + "/tests/worker";
    }
    {
      name = "tools";
      path = tools;
    }
    {
      name = "environment.json";
      path = environment;
    }
    {
      name = "compiler-environment.json";
      path = compilerEnvironment;
    }
    {
      name = "compiler-closure";
      path = compilerClosure;
    }
    {
      name = "closure";
      path = closure;
    }
  ]
