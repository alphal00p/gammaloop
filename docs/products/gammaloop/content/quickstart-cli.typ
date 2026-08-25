#import "../../shared.typ": callout

#let quickstart-cli = [
= Using GammaLoop from the command line

This guide builds one scalar one-loop topology, saves it as a GammaLoop state, and evaluates
one point. It uses a built-in model and does not require a run card, a UFO model, or FORM. Once the
executable is available, the example itself takes only a few seconds.

== Get the CLI

#callout("Nix is recommended, not required", [
  If you already use Nix, or want the most reproducible installation, follow the first route.
  Otherwise skip directly to the Cargo route below. Both routes provide the same `gammaloop`
  command, and neither requires a GammaLoop source checkout.
])

=== Recommended: the Nix package

The repository flake exposes GammaLoop as a packaged application. Install that package directly;
a source checkout and development shell are not required:

// docs-example: syntax
```sh
nix profile add github:alphal00p/gammaloop#gammaloop
gammaloop --help
```

Nix downloads the package from a configured binary cache when it is available and otherwise
builds the same derivation locally. The flake currently supports x86-64 and ARM64 Linux, plus
Apple-silicon macOS.

To try the same packaged CLI without adding it to your profile, run:

// docs-example: syntax
```sh
nix run github:alphal00p/gammaloop#gammaloop -- --help
```

Use `nix run .#gammaloop -- --help` instead when testing the package from an existing source
checkout. `nix develop` followed by `just build-cli` remains useful for contributors changing
GammaLoop, but it is not the normal installation path.

=== No Nix: install with Cargo

Install a stable Rust toolchain with #link("https://rustup.rs/")[rustup], then install the small
set of native build tools used by GammaLoop. On Debian or Ubuntu:

// docs-example: syntax
```sh
sudo apt-get update
sudo apt-get install -y build-essential m4 diffutils
```

On macOS, install Apple's command-line build tools if they are absent:

// docs-example: syntax
```sh
xcode-select --install
```

Then install the `gammaloop-api` crate. Its executable is named `gammaloop`:

// docs-example: syntax
```sh
SYMBOLICA_OEM_LICENSE=SYMBOLICA_OEM_GAMMALOOP \
  cargo install --locked gammaloop-api
gammaloop --help
```

This registry command becomes available with the first automated crates.io release. Until then,
install the release candidate directly from its tested Git revision, without making a checkout:

// docs-example: syntax
```sh
SYMBOLICA_OEM_LICENSE=SYMBOLICA_OEM_GAMMALOOP \
  cargo install --locked \
  --git https://github.com/alphal00p/gammaloop.git \
  --rev eb1dab23c24345def88226ef10874f0dde6c4aa5 \
  gammaloop-api
```

Cargo compiles the CLI locally, so its first installation takes longer than downloading a cached
Nix build. Its default build vendors the required GMP, MPFR, and MPC sources; OpenSSL, Python,
FORM, and system copies of those libraries are not prerequisites for this example. On Windows,
use WSL 2 and follow the Debian or Ubuntu instructions above. The `SYMBOLICA_OEM_LICENSE` value
is GammaLoop's public build selector, not a private license key. Add the `ufo_support` feature and
Python development headers only when external UFO model loading is needed. Symbolica's license
terms still apply to both installation paths.

#callout("Cargo install is Nix-free, but it is still a source build", [
  The current release workflow publishes crate sources and Python distributions, but not
  target-specific GammaLoop CLI archives or Cargo Binstall metadata. `cargo binstall` therefore
  cannot yet provide a supported prebuilt binary and may only fall back to a source build.
])

== Generate a one-loop graph

Users who installed the Nix package or Cargo crate can run this command as written. Contributors
using `just build-cli` in a source checkout should replace `gammaloop` with `./gammaloop`:

#callout("Choose a disposable state path", [
  `--clean-state` removes and replaces the state directory passed with `-s`. Use the fresh
  `./gammaloop_state/quickstart` path shown below, or substitute another path whose existing
  contents you are prepared to delete.
])

// docs-example: syntax gammaloop-quickstart-scalar-bubble
```sh
gammaloop --clean-state -s ./gammaloop_state/quickstart run -c \
  'import model scalars-default.json; set global kv global.generation.uv.generate_integrated=false; generate amp scalar_1 > scalar_1 [{1}] --allowed-vertex-interactions V_3_SCALAR_122 -p bubble -i one_loop; display processes; quit -o'
```

The command imports the built-in scalar model, generates one one-loop bubble amplitude, prints
the process named `bubble` and integrand named `one_loop`, and writes a reusable state under
`gammaloop_state/quickstart/`. Integrated UV generation is disabled deliberately so that this
first topology does not invoke FORM.

Check the result by looking for `bubble` and `one_loop` in the process listing and for
`gammaloop_state/quickstart/state_manifest.toml`. This graph is a compact software exercise,
not a normalized collider prediction.

== Evaluate one point

Load the state read-only and inspect the integrand at one loop-momentum point:

// docs-example: syntax
```sh
gammaloop --read-only-state -s ./gammaloop_state/quickstart run -c \
  'inspect -p bubble -i one_loop -x 0.1 0.2 0.3; quit'
```

A successful run prints the parameterization Jacobian and a finite complex value without
changing the saved state. Continue with the #link("tutorial/")[first-state tutorial] to use a
maintained run card, inspect the complete persistence lifecycle, and prepare a calculation whose
physics scope and normalization are explicit.
]
