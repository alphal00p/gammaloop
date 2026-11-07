# Gammaloop build and development commands

ci_cargo_profile := "ci-optim"

mod linnet 'crates/linnet/Justfile'
mod tydenso 'tydenso/Justfile'

# Imported recipes keep the same top-level commands and repository working directory.
import 'just/ci.just'
import 'just/tests.just'
import 'just/drawing.just'

alias sync-draw-assets := sync-drawing-assets

drawing_asset_sync := "none"

_sync-drawing-assets:
    #!/usr/bin/env bash
    set -euo pipefail
    case "{{ drawing_asset_sync }}" in
      nix)
        just sync-drawing-assets
        ;;
      cargo)
        just sync-drawing-assets-cargo
        ;;
      none|off|false)
        ;;
      *)
        echo "unknown drawing_asset_sync={{ drawing_asset_sync }}; expected nix, cargo, or none" >&2
        exit 2
        ;;
    esac

# Build gammaloop Python CLI with UFO support and dev-optim profile
build-cli: _sync-drawing-assets
    cargo build -p gammaloop-api --bin gammaloop --features ufo-import --profile dev-optim

# Build gammaloop CLI without Python (no pyo3)
build-cli-no-pyo3: _sync-drawing-assets
    cargo build -p gammaloop-api --bin gammaloop --no-default-features --features cli --profile dev-optim

# Build gammaloop Python CLI with UFO support and stable ABI (dev-optim profile)
build-cli-abi: _sync-drawing-assets
    cargo build -p gammaloop-api --bin gammaloop --features ufo-import,abi3 --profile dev-optim

# Build gammaloop Python CLI in release mode
build-cli-release: _sync-drawing-assets
    cargo build -p gammaloop-api --bin gammaloop --features ufo-import --release

# Build gammaloop Python CLI in release mode with stable ABI
build-cli-release-abi: _sync-drawing-assets
    cargo build -p gammaloop-api --bin gammaloop --features ufo-import,abi3 --release

# Build gammaloop Python API with UFO support and dev-optim profile
build-api: _sync-drawing-assets
    maturin develop -m crates/gammaloop-api/Cargo.toml --no-default-features --features=extension-module --profile=dev-optim

# Build gammaloop Python API with UFO support and stable ABI (dev-optim profile)
build-api-abi: _sync-drawing-assets
    maturin develop -m crates/gammaloop-api/Cargo.toml --no-default-features --features=extension-module,abi3 --profile=dev-optim

# Build gammaloop Python API with UFO support and dev-optim profile
build-api-wheel: _sync-drawing-assets
    maturin build -m crates/gammaloop-api/Cargo.toml --no-default-features --features=extension-module --profile=dev-optim

# Build gammaloop Python API with UFO support and stable ABI (dev-optim profile)
build-api-abi-wheel: _sync-drawing-assets
    maturin build -m crates/gammaloop-api/Cargo.toml --no-default-features --features=extension-module,abi3 --profile=dev-optim

# Build gammaloop Python API with UFO support and release profile
build-api-release: _sync-drawing-assets
    maturin develop -m crates/gammaloop-api/Cargo.toml --no-default-features --features=extension-module --profile=release

# Build gammaloop Python API with UFO support and stable ABI (release profile)
build-api-abi-release: _sync-drawing-assets
    maturin develop -m crates/gammaloop-api/Cargo.toml --no-default-features --features=extension-module,abi3 --profile=release

# Build gammaloop Python API wheel with UFO support and release profile
build-api-release-wheel: _sync-drawing-assets
    maturin build -m crates/gammaloop-api/Cargo.toml --no-default-features --features=extension-module --profile=release

# Build gammaloop Python API wheel with UFO support and stable ABI (release profile)
build-api-abi-release-wheel: _sync-drawing-assets
    maturin build -m crates/gammaloop-api/Cargo.toml --no-default-features --features=extension-module,abi3 --profile=release

# Build all workspace packages
build-all:
    cargo build --workspace

# Clean build artifacts
clean:
    cargo clean

# Guard mutually exclusive Symbolica backends in Hakari's native-only graph.
check-symbolica-feature-isolation:
    #!/usr/bin/env bash
    set -euo pipefail
    root="{{ justfile_directory() }}"
    hack="$root/crates/gammaloop-workspace-hack/Cargo.toml"
    if grep -Eq '"(wasm|integer-malachite|float-astro)"' "$hack"; then
      echo "the native GammaLoop workspace hack contains a Wasm Symbolica backend" >&2
      exit 1
    fi
    for feature in integer-gmp float-mpfr native_code_generation; do
      if ! grep -q "\"$feature\"" "$hack"; then
        echo "the native GammaLoop workspace hack is missing Symbolica feature $feature" >&2
        exit 1
      fi
    done
    symbolica_source="$(sed -n '/^name = "symbolica"$/,/^\[\[package\]\]/s/^source = "\(.*\)"$/\1/p' "$root/Cargo.lock")"
    if ! [[ "$symbolica_source" =~ ^git\+https://github\.com/symbolica-dev/symbolica\?branch=community#[0-9a-f]{40}$ ]]; then
      echo "the root lock must select one Symbolica revision from the official community branch" >&2
      exit 1
    fi
    for lock in "$root/Cargo.lock" "$root/tydenso/Cargo.lock"; do
      for package in symbolica numerica graphica; do
        source="$symbolica_source"
        version="3.0.1"
        expected="$(printf 'version = "%s"\nsource = "%s"\n' "$version" "$source")"
        resolved="$(sed -n "/^name = \"$package\"$/,/^\[\[package\]\]/p" "$lock" | grep -E '^(version|source) = ')"
        if [ "$resolved" != "$expected" ]; then
          echo "$lock must resolve exactly one $package at $version from $source; found ${resolved:-none}" >&2
          exit 1
        fi
      done
    done

# Keep Pyodide's portable arithmetic isolated from native code generation and Hakari.
check-spynso-wasm:
    #!/usr/bin/env bash
    set -euo pipefail
    dependencies=$(cargo tree --locked -p spynso3 --no-default-features --features wasm \
      --target wasm32-unknown-emscripten -e normal,build --prefix none -f '{p} {f}')
    if grep -E '^(gammaloop-workspace-hack |symbolica .*\b(integer-gmp|float-mpfr|native_code_generation)\b)' <<< "$dependencies"; then
      echo "Spynso's Pyodide dependencies enable native-only features" >&2
      exit 1
    fi
    PYO3_CROSS=1 PYO3_CROSS_PYTHON_VERSION=3.14 \
      cargo check --locked -p spynso3 --no-default-features --features wasm \
      --target wasm32-unknown-emscripten

# Check code without building
check: check-symbolica-feature-isolation
    cargo check --workspace --all-targets --locked

doc:
    cargo doc --workspace --no-deps --locked --profile {{ ci_cargo_profile }}

# Render the website illustrations and Linnest example into target/docs-svg-assets.
docs-svg-assets:
    scripts/render-docs-svg-assets.sh

# Compare rendered diagrams and marks with their reviewed snapshots.
docs-svg-assets-check:
    #!/usr/bin/env bash
    set -euo pipefail
    output=$(mktemp -d)
    trap 'rm -rf -- "$output"' EXIT
    bash scripts/render-docs-svg-assets.sh "$output" --check

# Build one product documentation site, or all five sites.
docs-site PRODUCT="all" CHANNEL="latest" SNAPSHOT_TAG="" OUTPUT="target/alphal00p-docs":
    #!/usr/bin/env bash
    set -euo pipefail

    channel={{ quote(CHANNEL) }}
    snapshot_tag={{ quote(SNAPSHOT_TAG) }}
    args=(
        build
        --product {{ quote(PRODUCT) }}
        --channel "$channel"
        --output {{ quote(OUTPUT) }}
    )

    case "$channel" in
        latest)
            if [ -n "$snapshot_tag" ]; then
                echo "SNAPSHOT_TAG is only valid when CHANNEL=snapshot" >&2
                exit 2
            fi
            ;;
        snapshot)
            if [ -z "$snapshot_tag" ]; then
                echo "SNAPSHOT_TAG is required when CHANNEL=snapshot" >&2
                exit 2
            fi
            args+=(--snapshot-tag "$snapshot_tag")
            ;;
        *)
            echo "CHANNEL must be latest or snapshot" >&2
            exit 2
            ;;
    esac

    unset GITHUB_REF_NAME
    just docs-check
    cargo run --locked -p alphal00p-docs-builder -- "${args[@]}"

# Add executable cells to a product in a previously built documentation site.
docs-notebooks WHEEL PRODUCT="linnet" OUTPUT="target/alphal00p-docs":
    uv run --no-project --with marimo==0.24.0 python \
        crates/linnet-py/examples/export_wasm.py --docs {{ quote(PRODUCT) }} \
        --wheel {{ quote(WHEEL) }} \
        --output {{ quote(OUTPUT + "/products/" + PRODUCT + "/latest/assets/notebooks") }}

# Validate the five-product documentation registry and generated inputs.
docs-check:
    cargo run --locked -p alphal00p-docs-catalogs --features gammaloop-reference --bin alphal00p-docs-gammaloop-reference -- --check
    cargo run --locked -p alphal00p-docs-catalogs --features vakint-reference --bin alphal00p-docs-vakint-reference -- --check
    cargo run --locked -p alphal00p-docs-python-exporter --features gammaloop -- gammaloop-python docs/api/python/gammaloop-python.pyi --check
    cargo run --locked -p alphal00p-docs-python-exporter --features linnet -- linnet-graph docs/api/python/linnet-graph.pyi --check
    cargo run --locked -p alphal00p-docs-python-exporter --features spenso -- spynso3 docs/api/python/spynso3.pyi --check
    cargo run --locked -p alphal00p-docs-python-exporter --features vakint -- vakint-community docs/api/python/vakint-community.pyi --check
    cargo test --locked -p alphal00p-docs-python-exporter
    cargo test --locked -p alphal00p-docs-python-exporter --features gammaloop gammaloop_runtime_surface_and_signatures_match_the_docs_stub
    cargo test --locked -p alphal00p-docs-python-exporter --features linnet linnet_package_and_docs_share_the_typed_stub_info_surface
    cargo test --locked -p alphal00p-docs-examples
    cargo run --locked -p alphal00p-docs-builder -- check
    just docs-svg-assets-check
    just docs-linnet-graph-check

# Render the front-page tour to target/showcase/ (needs Node, Playwright's Chromium, and ffmpeg).
showcase-video *ARGS:
    node scripts/render-showcase-video.mjs {{ ARGS }}

# Install the compiled Community graph module and run its native Python suites.
docs-linnet-graph-check:
    #!/usr/bin/env bash
    set -euo pipefail
    root="{{ justfile_directory() }}"
    scratch=$(mktemp -d)
    trap 'rm -rf -- "$scratch"' EXIT
    cd "$scratch"
    mkdir -p src python/symbolica/community
    cp -R "$root/crates/linnet-py/python/symbolica/community/graph" python/symbolica/community/
    touch python/symbolica/community/__init__.py
    echo 'from .core import *' >python/symbolica/__init__.py
    cat >pyproject.toml <<'TOML'
    [build-system]
    requires = ["maturin>=1,<2"]
    build-backend = "maturin"
    [project]
    name = "symbolica"
    version = "3.0.1"
    [tool.maturin]
    module-name = "symbolica.core"
    python-source = "python"
    TOML
    cp "$root/Cargo.lock" Cargo.lock
    python3 - "$root" <<'PYTHON'
    import json, pathlib, sys, tomllib
    root = pathlib.Path(sys.argv[1])
    workspace = tomllib.loads((root / "Cargo.toml").read_text())
    with open("Cargo.toml", "w") as fp:
        fp.write('[workspace]\n[package]\nname="linnet-graph-host"\nversion="0.1.0"\nedition="2024"\n[lib]\nname="core"\ncrate-type=["cdylib"]\n[dependencies]\n')
        fp.write('linnet-py = {path=' + json.dumps(str(root / "crates/linnet-py")) + '}\n')
        fp.write('pyo3 = {version="0.28",features=["abi3-py310","extension-module"]}\nsymbolica = {version="=3.0.1",default-features=false,features=["python_export","integer-gmp","float-mpfr"]}\n')
        fp.write('[patch.crates-io]\n')
        for name, props in workspace.get('patch', {}).get('crates-io', {}).items():
            props = dict(props)
            if 'path' in props:
                props['path'] = str(root / props['path'])
            fp.write(name + ' = {' + ', '.join(k + '=' + json.dumps(v) for k, v in props.items()) + '}\n')
    PYTHON
    cat >src/lib.rs <<'RUST'
    use pyo3::{prelude::*, wrap_pyfunction};
    use pyo3::types::PyModule;
    use symbolica::api::python::{create_symbolica_module, SymbolicaCommunityModule};
    #[pyfunction]
    fn initialize_module(py: Python) -> PyResult<()> { linnet_py::GraphModule::initialize(py) }
    #[pymodule]
    fn core(module: &Bound<'_, PyModule>) -> PyResult<()> {
        create_symbolica_module(module)?;
        let graph = PyModule::new(module.py(), "graph_native")?;
        graph.add_function(wrap_pyfunction!(initialize_module, &graph)?)?;
        linnet_py::GraphModule::register_module(&graph)?;
        module.add_submodule(&graph)?;
        module.py().import("sys")?.getattr("modules")?.set_item("symbolica.community.graph_native", graph)?;
        Ok(())
    }
    RUST
    cargo metadata --format-version 1 > /dev/null
    python3 - "$root/Cargo.lock" <<'PYTHON'
    import pathlib, sys, tomllib
    def packages(file):
        return {(p["name"], p["version"], p["source"]) for p in tomllib.loads(pathlib.Path(file).read_text())["package"] if "source" in p}
    extra = packages("Cargo.lock") - packages(sys.argv[1])
    if extra:
        raise SystemExit(f"test host resolved dependencies outside the workspace lock: {sorted(extra)}")
    PYTHON
    uv venv venv --python "${PYTHON_BIN_PATH:-python3}"
    maturin build --locked --out wheels
    uv pip install --python venv/bin/python --no-deps wheels/*.whl
    python_check="$scratch/venv/bin/python"
    cd "$root"
    for suite in test_basic.py test_wasm.py test_streaming.py; do
        "$python_check" -m unittest discover -s crates/linnet-py/tests -p "$suite"
    done

# Regenerate source-backed CLI/settings and topology/dependency snapshots.
docs-generated:
    cargo run --locked -p alphal00p-docs-catalogs --features gammaloop-reference --bin alphal00p-docs-gammaloop-reference
    cargo run --locked -p alphal00p-docs-catalogs --features vakint-reference --bin alphal00p-docs-vakint-reference

# Run the provisioned FORM-backed Vakint documentation examples without pySecDec.
docs-vakint-form-check:
    VAKINT_SKIP_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test tensor_reduction_tests
    VAKINT_SKIP_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_evaluation_freeform_tests test_integrate_1l_decorated_indices_matad -- --exact
    VAKINT_SKIP_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_evaluation_freeform_tests test_integrate_1l_decorated_indices_fmft -- --exact
    VAKINT_SKIP_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_alphaloop_vs_matad_tests

# Run one numerical pySecDec smoke comparison. The caller must provision pySecDec.
docs-vakint-pysecdec-smoke-check:
    RUN_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_evaluation_pysecdec_tests test_integrate_1l_simple -- --exact

# Run the complete opt-in pySecDec and analytic-backend comparison tier.
docs-vakint-pysecdec-full-check:
    RUN_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_evaluation_pysecdec_tests
    RUN_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --locked --package vakint --no-default-features --test integral_comparison_vs_pysecdec_tests

# Continuously build and serve one product with native browser reloads.
docs-watch PRODUCT PORT="8000" OUTPUT="target/alphal00p-docs-watch":
    #!/usr/bin/env bash
    set -euo pipefail

    exec cargo run --locked --profile docs-watch -p alphal00p-docs-builder \
        --features watch -- \
        watch \
        --product {{ quote(PRODUCT) }} \
        --port {{ quote(PORT) }} \
        --output {{ quote(OUTPUT) }}

# Backward-compatible name for the native live watcher.
docs-serve PRODUCT PORT="8000" OUTPUT="target/alphal00p-docs-serve":
    just docs-watch {{ quote(PRODUCT) }} {{ quote(PORT) }} {{ quote(OUTPUT) }}

# Format code
fmt *lint_args:
    cargo fmt --all {{ lint_args }}

# Run clippy linter
clippy *lint_args:
    cargo clippy --workspace --all-targets --locked --profile {{ ci_cargo_profile }} {{ lint_args }}

# Build everything in release mode for maximum performance
build-release-all:
    cargo build --workspace --release

# Build in release mode using the local Cargo artifacts.
dev-release: build-release-all

# Run gammaloop
run *ARGS:
    cargo run -p gammaloop-api --bin gammaloop --features ufo-import --profile dev-optim -- {{ ARGS }}
