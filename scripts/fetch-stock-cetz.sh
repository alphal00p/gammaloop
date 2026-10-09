#!/usr/bin/env bash
# Install stock CeTZ 0.5.2 from Typst Universe into a Typst package directory,
# for Typst-route tests outside the development shell (the wheels bundle no CeTZ).
# usage: scripts/fetch-stock-cetz.sh [PACKAGE_DIR]   (default: target/typst-packages)
set -euo pipefail

package_dir=${1:-target/typst-packages}
destination="$package_dir/preview/cetz/0.5.2"
archive=$(mktemp)
trap 'rm -f "$archive"' EXIT

curl -fsSL https://packages.typst.org/preview/cetz-0.5.2.tar.gz -o "$archive"
echo "77cf8490114ae04c6e665a11efa691d284a0cadb9719771b5708c1197292f23f  $archive" | sha256sum -c --quiet -
mkdir -p "$destination"
tar -xzf "$archive" -C "$destination"
