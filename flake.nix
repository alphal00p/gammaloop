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

      # Common arguments can be set here to avoid repeating them later
      commonArgs = {
        inherit src;
        strictDeps = true;

        buildInputs =
          [
            # Add additional build inputs here
            pkgs.openssl
            pkgs.pkg-config
            pkgs.gmp
            pkgs.mpfr
            pkgs.libmpc
            pkgs.python313
          ]
          ++ lib.optionals pkgs.stdenv.isDarwin [
            # Additional darwin specific inputs can be set here
            pkgs.libiconv
          ];

        nativeBuildInputs =
          [
            pkgs.pkg-config
            pkgs.clang
            pkgs.gmp.dev
            pkgs.mpfr.dev
            pkgs.libmpc
            # pkgs.gcc
            pkgs.git
            pkgs.python313
            pkgs.gnum4
            pkgs.gmp
            pkgs.mpfr
          ]
          ++ lib.optionals pkgs.stdenv.isDarwin [
            pkgs.darwin.cctools
          ];

        # Additional environment variables can be set directly
        # Help gmp-mpfr-sys find system libraries
        CPPFLAGS = "-I${pkgs.gmp.dev}/include -I${pkgs.mpfr.dev}/include -I${pkgs.libmpc}/include";
        LDFLAGS = "-L${pkgs.gmp}/lib -L${pkgs.mpfr}/lib -L${pkgs.libmpc}/lib";
        PKG_CONFIG_PATH = "${pkgs.gmp.dev}/lib/pkgconfig:${pkgs.mpfr.dev}/lib/pkgconfig:${pkgs.libmpc}/lib/pkgconfig";

        # Force gmp-mpfr-sys to use system libraries
        MPFR_LIB_DIR = "${pkgs.mpfr}/lib";
        MPFR_INCLUDE_DIR = "${pkgs.mpfr.dev}/include";
        GMP_LIB_DIR = "${pkgs.gmp}/lib";
        GMP_INCLUDE_DIR = "${pkgs.gmp.dev}/include";
        MPC_LIB_DIR = "${pkgs.libmpc}/lib";
        MPC_INCLUDE_DIR = "${pkgs.libmpc}/include";

        # Additional variables to force system lib usage
        GMP_MPFR_SYS_LIBRARY = "1";
        LIBGMP_LIB_DIR = "${pkgs.gmp}/lib";
        LIBGMP_INCLUDE_DIR = "${pkgs.gmp.dev}/include";
        LIBMPFR_LIB_DIR = "${pkgs.mpfr}/lib";
        LIBMPFR_INCLUDE_DIR = "${pkgs.mpfr.dev}/include";
        LIBMPC_LIB_DIR = "${pkgs.libmpc}/lib";
        LIBMPC_INCLUDE_DIR = "${pkgs.libmpc}/include";

        # Disable source builds
        GMP_NO_SYS = "0";
        MPFR_NO_SYS = "0";
        MPC_NO_SYS = "0";

        # Rug-specific environment variables for system GMP
        CARGO_FEATURE_USE_SYSTEM_LIBS = "1";
        RUG_GMP_DIR = "${pkgs.gmp}";
        RUG_MPFR_DIR = "${pkgs.mpfr}";
        RUG_MPC_DIR = "${pkgs.libmpc}";
      };

      # CI-specific args that disable Python API to avoid pyo3 build issues
      ciArgs =
        commonArgs
        // {
          buildType = "release";
          # Set PyO3 environment variables to help it find Python
          PYO3_PYTHON = "${pkgs.python313}/bin/python3";
          PYTHONPATH = "${pkgs.python313}/lib/python3.13/site-packages";
          # Help gmp-mpfr-sys find system libraries
          CPPFLAGS = "-I${pkgs.gmp.dev}/include -I${pkgs.mpfr.dev}/include -I${pkgs.libmpc}/include";
          LDFLAGS = "-L${pkgs.gmp}/lib -L${pkgs.mpfr}/lib -L${pkgs.libmpc}/lib";
          PKG_CONFIG_PATH = "${pkgs.gmp.dev}/lib/pkgconfig:${pkgs.mpfr.dev}/lib/pkgconfig:${pkgs.libmpc}/lib/pkgconfig";

          # Force gmp-mpfr-sys to use system libraries
          MPFR_LIB_DIR = "${pkgs.mpfr}/lib";
          MPFR_INCLUDE_DIR = "${pkgs.mpfr.dev}/include";
          GMP_LIB_DIR = "${pkgs.gmp}/lib";
          GMP_INCLUDE_DIR = "${pkgs.gmp.dev}/include";
          MPC_LIB_DIR = "${pkgs.libmpc}/lib";
          MPC_INCLUDE_DIR = "${pkgs.libmpc}/include";

          # Additional variables to force system lib usage
          GMP_MPFR_SYS_LIBRARY = "1";
          LIBGMP_LIB_DIR = "${pkgs.gmp}/lib";
          LIBGMP_INCLUDE_DIR = "${pkgs.gmp.dev}/include";
          LIBMPFR_LIB_DIR = "${pkgs.mpfr}/lib";
          LIBMPFR_INCLUDE_DIR = "${pkgs.mpfr.dev}/include";
          LIBMPC_LIB_DIR = "${pkgs.libmpc}/lib";
          LIBMPC_INCLUDE_DIR = "${pkgs.libmpc}/include";

          # Disable source builds
          GMP_NO_SYS = "0";
          MPFR_NO_SYS = "0";
          MPC_NO_SYS = "0";

          # Rug-specific environment variables for system GMP
          CARGO_FEATURE_USE_SYSTEM_LIBS = "1";
          RUG_GMP_DIR = "${pkgs.gmp}";
          RUG_MPFR_DIR = "${pkgs.mpfr}";
          RUG_MPC_DIR = "${pkgs.libmpc}";
        };

      craneLibLLvmTools =
        craneLib.overrideToolchain
        (fenix.packages.${system}.stable.withComponents [
          "cargo"
          "llvm-tools"
          "rustc"
        ]);

      # Build *just* the cargo dependencies, so we can reuse
      # all of that work (e.g. via cachix) when running in CI
      cargoArtifacts = craneLib.buildDepsOnly ciArgs;

      # Build the actual crate itself, reusing the dependency
      # artifacts from above.
      gammaloop = craneLib.buildPackage (commonArgs
        // {
          inherit cargoArtifacts;
        });
    in {
      checks = {
        # Build the crate as part of `nix flake check` for convenience
        inherit gammaloop;

        # Run clippy (and deny all warnings) on the crate source,
        # again, reusing the dependency artifacts from above.
        #
        # Note that this is done as a separate derivation so that
        # we can block the CI if there are issues here, but not
        # prevent downstream consumers from building our crate by itself.
        gammaloop-clippy = craneLib.cargoClippy (ciArgs
          // {
            inherit cargoArtifacts;
            cargoClippyExtraArgs = "--all-targets -- --deny warnings";
          });

        gammaloop-doc = craneLib.cargoDoc (ciArgs
          // {
            inherit cargoArtifacts;
          });

        # Check formatting
        gammaloop-fmt = craneLib.cargoFmt {
          inherit src;
        };

        # Audit dependencies
        gammaloop-audit = craneLib.cargoAudit {
          inherit src advisory-db;
        };

        # Audit licenses
        gammaloop-deny = craneLib.cargoDeny {
          inherit src;
        };

        # Run tests with cargo-nextest
        # Consider setting `doCheck = false` on `gammaloop` if you do not want
        # the tests to run twice
        gammaloop-nextest = craneLib.cargoNextest (ciArgs
          // {
            inherit cargoArtifacts;
            cargoNextestExtraArgs = "--profile ci --test-threads 0 --no-fail-fast --final-status-level fail";
          });

        # Individual partitioned test checks for CI
        gammaloop-nextest-partition-1 = craneLib.cargoNextest (ciArgs
          // {
            inherit cargoArtifacts;
            cargoNextestExtraArgs = "--profile ci --test-threads 0 --no-fail-fast --final-status-level fail --partition hash:1/6";
          });

        gammaloop-nextest-partition-2 = craneLib.cargoNextest (ciArgs
          // {
            inherit cargoArtifacts;
            cargoNextestExtraArgs = "--profile ci --test-threads 0 --no-fail-fast --final-status-level fail --partition hash:2/6";
          });

        gammaloop-nextest-partition-3 = craneLib.cargoNextest (ciArgs
          // {
            inherit cargoArtifacts;
            cargoNextestExtraArgs = "--profile ci --test-threads 0 --no-fail-fast --final-status-level fail --partition hash:3/6";
          });

        gammaloop-nextest-partition-4 = craneLib.cargoNextest (ciArgs
          // {
            inherit cargoArtifacts;
            cargoNextestExtraArgs = "--profile ci --test-threads 0 --no-fail-fast --final-status-level fail --partition hash:4/6";
          });

        gammaloop-nextest-partition-5 = craneLib.cargoNextest (ciArgs
          // {
            inherit cargoArtifacts;
            cargoNextestExtraArgs = "--profile ci --test-threads 0 --no-fail-fast --final-status-level fail --partition hash:5/6";
          });

        gammaloop-nextest-partition-6 = craneLib.cargoNextest (ciArgs
          // {
            inherit cargoArtifacts;
            cargoNextestExtraArgs = "--profile ci --test-threads 0 --no-fail-fast --final-status-level fail --partition hash:6/6";
          });
      };

      packages =
        {
          default = gammaloop;
          # Expose cargoArtifacts for CI caching
          inherit cargoArtifacts;
        }
        // lib.optionalAttrs (!pkgs.stdenv.isDarwin) {
          gammaloop-llvm-coverage = craneLibLLvmTools.cargoLlvmCov (commonArgs
            // {
              inherit cargoArtifacts;
            });
        };

      apps = {
        default = flake-utils.lib.mkApp {
          drv = gammaloop;
        };
      };

      devShells.default = craneLib.devShell {
        # Inherit inputs from checks.
        checks = self.checks.${system};

        # Additional dev-shell environment variables can be set directly
        # MY_CUSTOM_DEVELOPMENT_VAR = "something else";
        RUST_SRC_PATH = "${pkgs.rustPlatform.rustLibSrc}";
        GLIBC_TUNABLES = "glibc.rtld.optional_static_tls=10000";

        # LD_LIBRARY_PATH = "${pkgs.stdenv.cc.cc.lib}/lib";

        # Extra inputs can be added here; cargo and rustc are provided by default.
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
          # pkgs.ripgrep
          cargo-insta
          cargo-udeps
          cargo-insta
          openssl
          pyright
          gmp
          mpfr
          form
          gnum4
          nickel
          nls
          typst
          cargo-nextest
          # gcc_debug.out
          # stdenv.cc.cc.lib
          pkg-config
          cargo-deny
          cargo-edit
          cargo-watch
          bacon
          # python311
          texlive.combined.scheme-medium
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
