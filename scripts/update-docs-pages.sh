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
        reference/rust/index.html \
        reference/python/index.html
    do
        [ -f "$bundle/$required" ] || fail "incomplete product bundle: $bundle/$required is missing"
    done
    find "$bundle/reference/python" -mindepth 2 -maxdepth 2 -name index.html -print -quit |
        grep -q . || fail "incomplete product bundle: $bundle has no structured Python component page"
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

        for product in "${products[@]}"; do
            require_bundle "$build_root/products/$product/latest"
            [ -f "$build_root/products/$product/index.html" ] ||
                fail "latest build has no redirect for $product"
        done

        install -m 0644 "$build_root/index.html" "$pages_root/index.html"
        install -m 0644 "$build_root/.nojekyll" "$pages_root/.nojekyll"
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
            require_bundle "$source_bundle"
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
