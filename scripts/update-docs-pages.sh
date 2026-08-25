#!/usr/bin/env bash
set -euo pipefail

readonly products=(gammaloop linnet spenso idenso vakint)

usage() {
    echo "usage: $0 <latest|snapshot> <build-root> <gh-pages-worktree> [snapshot-tag]" >&2
    exit 2
}

fail() {
    echo "update-docs-pages: $*" >&2
    exit 1
}

canonical_directory() {
    (cd "$1" && pwd -P)
}

require_generated_reference_surface() {
    local bundle=$1
    local product=$2
    python3 - "$bundle" "$product" <<'PY'
import json
import sys
from pathlib import Path

bundle = Path(sys.argv[1])
product = sys.argv[2]


def reject(message):
    print(f"update-docs-pages: incomplete product bundle: {message}", file=sys.stderr)
    raise SystemExit(1)


def load_json(path):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        reject(f"cannot read {path}: {error}")


def pages_below(route):
    root = bundle / route
    return {
        f"{path.parent.relative_to(bundle).as_posix()}/"
        for path in root.rglob("index.html")
    } if root.is_dir() else set()


def require_search_coverage(label, pages, expected_count, kind):
    if len(pages) != expected_count:
        reject(f"{label} has {len(pages)} pages, expected {expected_count}")
    indexed = {entry.get("href") for entry in search if entry.get("kind") == kind}
    if missing := pages - indexed:
        reject(f"{label} pages are absent from {kind} search: {sorted(missing)}")


def supported_symbols(scope):
    return sum(item.get("supported") is True for item in scope.get("items", {}).values()) + sum(
        supported_symbols(child) for child in scope.get("scopes", {}).values()
    )


search = load_json(bundle / "search-index.json")
if not isinstance(search, list):
    reject(f"{bundle / 'search-index.json'} is not an entry list")

catalogs = [load_json(path) for path in sorted((bundle / "catalogs").glob("*.json"))]

rust_catalogs = [
    catalog for catalog in catalogs if catalog.get("component", {}).get("language") == "rust"
]
if not rust_catalogs:
    reject(f"{bundle / 'catalogs'} has no Rust catalog")
for catalog in rust_catalogs:
    package = catalog["component"]["package"]
    crate_roots = [
        entry.get("href")
        for entry in search
        if entry.get("kind") == "rust-api"
        and entry.get("title") == f"{package} Rustdoc"
    ]
    if len(crate_roots) != 1:
        reject(f"{package} has {len(crate_roots)} crate-root search entries, expected 1")
    crate_root = crate_roots[0]
    route_parts = Path(crate_root).parts if isinstance(crate_root, str) else ()
    if (
        not isinstance(crate_root, str)
        or not crate_root.startswith("reference/rust/")
        or not crate_root.endswith("/")
        or len(route_parts) != 3
        or ".." in route_parts
    ):
        reject(f"{package} has invalid crate-root route: {crate_root!r}")
    crate_root_page = bundle / crate_root / "index.html"
    if not crate_root_page.is_file():
        reject(f"{crate_root_page} is missing")

python_catalogs = [
    catalog for catalog in catalogs if catalog.get("component", {}).get("language") == "python"
]
if not python_catalogs:
    reject(f"{bundle / 'catalogs'} has no Python catalog")

supported_python = 0
for catalog in python_catalogs:
    component = catalog["component"]["id"]
    module_route = f"reference/python/{component}/"
    pages = pages_below(module_route)
    if module_route not in pages:
        reject(f"{bundle / module_route / 'index.html'} is missing")
    symbol_count = supported_symbols(catalog["root"])
    supported_python += symbol_count
    require_search_coverage(
        f"{component} Python symbols", pages - {module_route}, symbol_count, "python-api"
    )
if not supported_python:
    reject(f"{bundle / 'reference/python'} has no Python symbol page")

if product != "gammaloop":
    raise SystemExit(0)

reference = load_json(bundle / "reference/generated/gammaloop-reference.json")
public_commands = sum(
    not command.get("hidden") and not command.get("generated_help")
    for command in reference["commands"]
)
require_search_coverage(
    "public CLI commands", pages_below("reference/cli/commands"), public_commands, "command"
)

setting_paths = [setting["path"] for setting in reference["settings"]]
setting_titles = [entry.get("title") for entry in search if entry.get("kind") == "setting"]
if len(setting_titles) != len(setting_paths) or set(setting_titles) != set(setting_paths):
    reject("setting search titles do not match the generated settings schema")

namespaces = set()
for path in setting_paths:
    parts = path.split(".")
    namespaces.update(".".join(parts[:end]) for end in range(1, len(parts)))
namespace_pages = pages_below("reference/cli/settings") - {"reference/cli/settings/"}
require_search_coverage(
    "settings namespaces", namespace_pages, len(namespaces), "settings namespace"
)
PY
}

require_bundle() {
    local bundle=$1
    local product=$2
    local required
    for required in \
        .note \
        index.html \
        manual.pdf \
        snapshot.json \
        search-index.json \
        tutorial/index.html \
        reference/interfaces/index.html \
        version-history/index.html \
        manual/interfaces/index.html \
        manual/releases/index.html \
        assets/site.css \
        assets/site.js \
        assets/local-unitarity-light.svg \
        assets/local-unitarity-dark.svg \
        assets/spensologo.svg \
        assets/gammalooplogo-light.svg \
        assets/gammalooplogo-dark.svg \
        reference/index.html \
        reference/rust/index.html \
        reference/python/index.html
    do
        [ -f "$bundle/$required" ] || fail "incomplete product bundle: $bundle/$required is missing"
    done
    case "$product" in
        gammaloop)
            [ -f "$bundle/reference/cli/index.html" ] ||
                fail "incomplete product bundle: $bundle has no generated CLI reference"
            [ -f "$bundle/reference/cli/settings/index.html" ] ||
                fail "incomplete product bundle: $bundle has no settings namespace index"
            [ -f "$bundle/reference/generated/gammaloop-reference.json" ] ||
                fail "incomplete product bundle: $bundle has no machine-readable CLI catalogue"
            ;;
        vakint)
            [ -f "$bundle/reference/topologies/index.html" ] ||
                fail "incomplete product bundle: $bundle has no generated topology reference"
            [ -f "$bundle/reference/generated/vakint-reference.json" ] ||
                fail "incomplete product bundle: $bundle has no machine-readable topology catalogue"
            ;;
    esac
    require_generated_reference_surface "$bundle" "$product"
}

[ "$#" -ge 3 ] && [ "$#" -le 4 ] || usage

mode=$1
build_root_input=$2
pages_root_input=$3
snapshot_tag=${4-}

[ -d "$build_root_input" ] || fail "build root is not a directory: $build_root_input"
[ -d "$pages_root_input" ] || fail "Pages worktree is not a directory: $pages_root_input"
[ -e "$pages_root_input/.git" ] || fail "Pages destination is not a Git worktree: $pages_root_input"

build_root=$(canonical_directory "$build_root_input")
pages_root=$(canonical_directory "$pages_root_input")
[ "$build_root" != "$pages_root" ] || fail "build root and Pages worktree must be different directories"

case "$mode" in
    latest)
        [ "$#" -eq 3 ] || usage
        [ -f "$build_root/index.html" ] || fail "latest build has no repository portal"
        [ -f "$build_root/search-index.json" ] ||
            fail "latest build has no federated search index"
        [ -f "$build_root/.nojekyll" ] || fail "latest build has no .nojekyll marker"
        [ -f "$build_root/assets/site.css" ] || fail "latest build has no portal stylesheet"
        [ -f "$build_root/assets/site.js" ] || fail "latest build has no portal script"
        [ -f "$build_root/assets/local-unitarity-light.svg" ] ||
            fail "latest build has no light collaboration mark"
        [ -f "$build_root/assets/local-unitarity-dark.svg" ] ||
            fail "latest build has no dark collaboration mark"
        for asset in \
            about-double-triangle-light.svg \
            about-double-triangle-dark.svg \
            about-local-unitarity-equation-light.svg \
            about-local-unitarity-equation-dark.svg; do
            [ -f "$build_root/assets/$asset" ] ||
                fail "latest build has no About-page Typst asset: $asset"
        done
        for graph_id in \
            aa-2l-gl00 \
            aa-2l-gl08 \
            aa-3l-gl000 \
            aa-3l-gl150 \
            aa-3l-gl300 \
            gg-hhh-3l \
            gg-hhh-1l \
            qq-aaa-pentabox \
            ad-ad-gluon \
            epem-bbx \
            epem-ttbar-cut; do
            for theme in light dark; do
                graph="portal-graph-$graph_id-$theme.svg"
                [ -f "$build_root/assets/graphs/$graph" ] ||
                    fail "latest build has no Typst graph asset: $graph"
            done
        done
        [ -f "$build_root/assets/gammalooplogo-light.svg" ] ||
            fail "latest build has no light project wordmark"
        [ -f "$build_root/assets/gammalooplogo-dark.svg" ] ||
            fail "latest build has no dark project wordmark"
        [ -f "$build_root/assets/spensologo.svg" ] ||
            fail "latest build has no Spenso project wordmark"
        [ -f "$build_root/assets/publications.json" ] ||
            fail "latest build has no publication cache"
        [ -f "$build_root/assets/people/valentin.webp" ] ||
            fail "latest build has no people portraits"
        [ -f "$build_root/people/index.html" ] ||
            fail "latest build has no people page"
        [ -f "$build_root/about/index.html" ] ||
            fail "latest build has no about page"
        [ -f "$build_root/talks/index.html" ] ||
            fail "latest build has no talks page"
        [ -f "$build_root/publications/index.html" ] ||
            fail "latest build has no publications page"
        [ -f "$build_root/citations/index.html" ] ||
            fail "latest build has no citations page"
        [ -f "$build_root/developers/.note" ] ||
            fail "latest build has no developer publication note"
        [ -f "$build_root/developers/index.html" ] ||
            fail "latest build has no developer architecture landing"
        [ -f "$build_root/developers/search-index.json" ] ||
            fail "latest build has no developer search index"
        [ -f "$build_root/developers/assets/site.css" ] ||
            fail "latest build has no developer stylesheet"
        [ -f "$build_root/developers/assets/site.js" ] ||
            fail "latest build has no developer script"
        [ -f "$build_root/developers/architecture/gammaloop-architecture/index.html" ] ||
            fail "latest build has no implemented GammaLoop architecture note"
        for product in linnet spenso idenso vakint; do
            [ -f "$build_root/developers/architecture/$product-architecture/index.html" ] ||
                fail "latest build has no implemented $product architecture note"
        done
        [ -f "$build_root/developers/architecture/spenso-parsing-flow/index.html" ] ||
            fail "latest build has no native Spenso parsing page"
        [ ! -e "$build_root/developers/architecture/spenso-parsing-flow/diagram.html" ] ||
            fail "latest build still contains the legacy Spenso parsing diagram"

        for product in "${products[@]}"; do
            require_bundle "$build_root/products/$product/latest" "$product"
            [ -f "$build_root/products/$product/index.html" ] ||
                fail "latest build has no redirect for $product"
        done

        install -m 0644 "$build_root/index.html" "$pages_root/index.html"
        install -m 0644 "$build_root/search-index.json" "$pages_root/search-index.json"
        install -m 0644 "$build_root/.nojekyll" "$pages_root/.nojekyll"
        rm -rf -- "$pages_root/assets"
        cp -a "$build_root/assets" "$pages_root/assets"
        for portal_page in about people talks publications citations; do
            rm -rf -- "$pages_root/$portal_page"
            cp -a "$build_root/$portal_page" "$pages_root/$portal_page"
        done
        rm -rf -- "$pages_root/developers"
        cp -a "$build_root/developers" "$pages_root/developers"
        mkdir -p "$pages_root/products"

        for product in "${products[@]}"; do
            product_root="$pages_root/products/$product"
            mkdir -p "$product_root"
            rm -rf -- "$product_root/latest"
            cp -a "$build_root/products/$product/latest" "$product_root/latest"
            install -m 0644 \
                "$build_root/products/$product/index.html" \
                "$product_root/index.html"
        done
        ;;
    snapshot)
        [ "$#" -eq 4 ] || usage
        [ -n "$snapshot_tag" ] || usage
        [[ "$snapshot_tag" =~ ^[A-Za-z0-9][A-Za-z0-9._+-]*$ ]] ||
            fail "unsafe snapshot tag: $snapshot_tag"
        [[ "$snapshot_tag" != *..* ]] || fail "unsafe snapshot tag: $snapshot_tag"

        for product in "${products[@]}"; do
            source_bundle="$build_root/products/$product/snapshots/$snapshot_tag"
            destination_bundle="$pages_root/products/$product/snapshots/$snapshot_tag"
            require_bundle "$source_bundle" "$product"
            if [ -e "$destination_bundle" ] &&
                ! diff -qr "$source_bundle" "$destination_bundle" >/dev/null
            then
                diff -qr "$source_bundle" "$destination_bundle" >&2 || true
                fail "immutable snapshot differs: products/$product/snapshots/$snapshot_tag"
            fi
        done

        # A snapshot publication never changes the portal, redirects, or latest
        # documentation. Check every existing snapshot before copying any new
        # bundle so a conflict leaves the worktree untouched.
        for product in "${products[@]}"; do
            source_bundle="$build_root/products/$product/snapshots/$snapshot_tag"
            destination_bundle="$pages_root/products/$product/snapshots/$snapshot_tag"
            if [ -e "$destination_bundle" ]; then
                echo "Snapshot already exists and is byte-identical: $product/$snapshot_tag"
            else
                mkdir -p "$(dirname "$destination_bundle")"
                cp -a "$source_bundle" "$destination_bundle"
            fi
        done
        ;;
    *)
        usage
        ;;
esac
