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
    find "$bundle/reference/python" -mindepth 2 -maxdepth 2 -name index.html -print -quit |
        grep -q . || fail "incomplete product bundle: $bundle has no structured Python component page"
    find "$bundle/reference/rust/supported" -mindepth 2 -maxdepth 2 -name index.html -print -quit |
        grep -q . || fail "incomplete product bundle: $bundle has no structured Rust component page"

    case "$product" in
        gammaloop)
            [ -f "$bundle/reference/cli/index.html" ] ||
                fail "incomplete product bundle: $bundle has no generated CLI reference"
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
        [ -f "$build_root/.nojekyll" ] || fail "latest build has no .nojekyll marker"
        [ -f "$build_root/assets/site.css" ] || fail "latest build has no portal stylesheet"
        [ -f "$build_root/assets/site.js" ] || fail "latest build has no portal script"
        [ -f "$build_root/assets/local-unitarity-light.svg" ] ||
            fail "latest build has no light collaboration mark"
        [ -f "$build_root/assets/local-unitarity-dark.svg" ] ||
            fail "latest build has no dark collaboration mark"
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
        [ -f "$build_root/developers/architecture/spenso-parsing-flow/diagram.html" ] ||
            fail "latest build has no Spenso parsing diagram"

        for product in "${products[@]}"; do
            require_bundle "$build_root/products/$product/latest" "$product"
            [ -f "$build_root/products/$product/index.html" ] ||
                fail "latest build has no redirect for $product"
        done

        install -m 0644 "$build_root/index.html" "$pages_root/index.html"
        install -m 0644 "$build_root/.nojekyll" "$pages_root/.nojekyll"
        rm -rf -- "$pages_root/assets"
        cp -a "$build_root/assets" "$pages_root/assets"
        for portal_page in people publications citations; do
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
