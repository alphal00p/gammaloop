{
  description = "Gammaloop";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

    crane = {
      url = "github:ipetkov/crane";
      # inputs.nixpkgs.follows = "nixpkgs";
    };

    fenix = {
      url = "github:nix-community/fenix";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.rust-analyzer-src.follows = "";
    };

    flake-utils.url = "github:numtide/flake-utils";

    advisory-db = {
      url = "github:rustsec/advisory-db";
      flake = false;
    };
  };

  outputs = {
    self,
    nixpkgs,
    crane,
    fenix,
    flake-utils,
    advisory-db,
    ...
  }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = nixpkgs.legacyPackages.${system};
      inherit (pkgs) lib;

      craneLib =
        (crane.mkLib nixpkgs.legacyPackages.${system}).overrideToolchain
        fenix.packages.${system}.stable.toolchain;

      src = craneLib.cleanCargoSource (craneLib.path ./.);

      # Host Rust target triple, e.g. x86_64-unknown-linux-gnu
      rustTarget = pkgs.stdenv.hostPlatform.rust.rustcTarget;

      # Env var name Cargo uses to pick the linker for this target
      cargoLinkerVar = "CARGO_TARGET_${lib.toUpper (lib.replaceStrings ["-"] ["_"] rustTarget)}_LINKER";

      # Force the "Nix cc wrapper" as both C compiler and Rust linker.
      nixCc = "${pkgs.stdenv.cc}/bin/cc";
      nixCxx = "${pkgs.stdenv.cc}/bin/c++";

      # Runtime library search path for locally-built binaries and for maturin/auditwheel
      # (so it can locate libpython and libs like libmpfr at repair time).
      runtimeLibPath = lib.makeLibraryPath [
        pkgs.python313
        pkgs.gmp
        pkgs.mpfr
        pkgs.libmpc
        pkgs.openssl
        pkgs.stdenv.cc.cc.lib
      ];

      # Common arguments can be set here to avoid repeating them later
      commonArgs = {
        inherit src;
        strictDeps = true;

        nativeBuildInputs =
          [
            pkgs.pkg-config
            pkgs.gcc
            pkgs.git
            pkgs.python313
          ]
          ++ lib.optionals pkgs.stdenv.isDarwin [
            pkgs.darwin.cctools
          ];

        buildInputs =
          [
            pkgs.openssl

            # System GMP/MPFR/MPC (for gmp-mpfr-sys with feature use-system-libs)
            pkgs.gmp
            pkgs.gmp.dev
            pkgs.mpfr
            pkgs.mpfr.dev
            pkgs.libmpc

            pkgs.python313
          ]
          ++ lib.optionals pkgs.stdenv.isDarwin [
            pkgs.libiconv
          ];

        # Hard override: use Nix's compiler wrapper everywhere.
        CC = nixCc;
        CXX = nixCxx;

        # Hard override: use Nix's cc wrapper as Rust linker for the host target.
        "${cargoLinkerVar}" = nixCc;

        # Final override in case something injects `-C linker=clang`.
        RUSTFLAGS = "-C linker=${nixCc}";

        # Make runtime libs discoverable inside Nix builds too (useful for build scripts/tests
        # that execute freshly-built binaries).
        LD_LIBRARY_PATH = runtimeLibPath;
        DYLD_LIBRARY_PATH = runtimeLibPath;
      };

      ciArgs =
        commonArgs
        // {
          buildType = "release";

          PYO3_PYTHON = "${pkgs.python313}/bin/python3";
          PYTHONPATH = "${pkgs.python313}/lib/python3.13/site-packages";
        };

      craneLibLLvmTools =
        craneLib.overrideToolchain
        (fenix.packages.${system}.stable.withComponents [
          "cargo"
          "llvm-tools"
          "rustc"
        ]);

      cargoArtifacts = craneLib.buildDepsOnly ciArgs;

      gammaloop = craneLib.buildPackage (commonArgs
        // {
          inherit cargoArtifacts;
        });
    in {
      checks = {
        inherit gammaloop;

        gammaloop-clippy = craneLib.cargoClippy (ciArgs
          // {
            inherit cargoArtifacts;
            cargoClippyExtraArgs = "--all-targets -- --deny warnings";
          });

        gammaloop-doc = craneLib.cargoDoc (ciArgs
          // {
            inherit cargoArtifacts;
          });

        gammaloop-fmt = craneLib.cargoFmt {
          inherit src;
        };

        gammaloop-audit = craneLib.cargoAudit {
          inherit src advisory-db;
        };

        gammaloop-deny = craneLib.cargoDeny {
          inherit src;
        };

        gammaloop-nextest = craneLib.cargoNextest (ciArgs
          // {
            inherit cargoArtifacts;
            cargoNextestExtraArgs = "--profile ci --test-threads 0 --no-fail-fast --final-status-level fail";
          });
      };

      packages =
        {
          default = gammaloop;
          inherit cargoArtifacts;
        }
        // lib.optionalAttrs (!pkgs.stdenv.isDarwin) {
          gammaloop-llvm-coverage = craneLibLLvmTools.cargoLlvmCov (commonArgs
            // {
              inherit cargoArtifacts;
            });
        };

      apps.default = flake-utils.lib.mkApp {
        drv = gammaloop;
      };

      devShells.default = craneLib.devShell {
        checks = self.checks.${system};

        RUST_SRC_PATH = "${pkgs.rustPlatform.rustLibSrc}";
        GLIBC_TUNABLES = "glibc.rtld.optional_static_tls=10000";

        # Mirror the same hard overrides in the interactive shell.
        CC = nixCc;
        CXX = nixCxx;
        "${cargoLinkerVar}" = nixCc;
        RUSTFLAGS = "-C linker=${nixCc}";

        # Make libpython + libgmp/libmpfr/libmpc visible to:
        # - target/debug/stub_gen (runtime)
        # - maturin/auditwheel (wheel repair step)
        LD_LIBRARY_PATH = runtimeLibPath;
        DYLD_LIBRARY_PATH = runtimeLibPath;

        packages = with pkgs; [
          tdf
          cargo-flamegraph
          yaml-language-server
          zellij
          vim
          helix
          jujutsu
          just
          dot-language-server
          cargo-insta
          cargo-udeps
          cargo-machete
          openssl
          pyright
          gmp
          mpfr
          libmpc
          form
          gnum4
          nickel
          nls
          typst
          cargo-nextest
          pkg-config
          cargo-deny
          cargo-edit
          cargo-watch
          bacon
          gfortran
          gcc
          rust-script
          uv
          ghostscript
          graphviz
          mupdf
          tinymist
          typstyle
          poppler-utils
          rust-analyzer
          maturin
          virtualenv
          (pkgs.rustPlatform.buildRustPackage rec {
            pname = "clinnet";
            version = "0.1.8";
            src = pkgs.fetchCrate {
              inherit pname version;
              sha256 = "sha256-CbZBHbf+8bIkdiSI5LMFO2Qc3zDr9UEBEry+fZOuep8=";
            };
            cargoHash = "sha256-GTixU2ZJZVMrEWLOfWjEnXMVLG2+cpkPbJuNnkTuFfo=";
          })
        ];
      };
    });
}
