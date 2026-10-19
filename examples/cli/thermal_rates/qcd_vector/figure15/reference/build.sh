#!/usr/bin/env bash
set -euo pipefail
reference_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build_directory="${1:?usage: build.sh BUILD_DIRECTORY}"
mkdir -p "$build_directory"
build_directory="$(cd -- "$build_directory" && pwd)"
archive="$build_directory/spectral-v0.9-spf.zip"
archive_sha256=207a2cb3d39933a57a5c1a8592d664c7dc9d5750eb8cd10801e98ac9560446fb
curl --fail --location https://zenodo.org/api/records/3478144/files/gw3g/spectral-v0.9-spf.zip/content --output "$archive"
printf '%s  %s\n' "$archive_sha256" "$archive" | sha256sum --check -
unzip -qo "$archive" -d "$build_directory"
spectral_source="$build_directory/gw3g-spectral-dca510d"
gsl_config="${GSL_CONFIG:-gsl-config}"
gsl_cflags_text="$("$gsl_config" --cflags)"
read -r -a gsl_cflags <<< "$gsl_cflags_text"
gsl_libraries_text="$("$gsl_config" --libs)"
read -r -a gsl_libraries <<< "$gsl_libraries_text"
gsl_libdir="$("$gsl_config" --prefix)/lib"
compile=("${CXX:-c++}" -std=c++11 -O2 "${gsl_cflags[@]}" -I "$spectral_source" -I "$spectral_source/inc" "$spectral_source"/src/*.cpp "$reference_directory/driver.cpp" "${gsl_libraries[@]}" "-Wl,-rpath,$gsl_libdir" -o "$build_directory/reference")
printf '%q ' "${compile[@]}" > "$build_directory/build-command.txt"
printf '\n' >> "$build_directory/build-command.txt"
"${compile[@]}"
"${CXX:-c++}" --version > "$build_directory/compiler-version.txt"
"$gsl_config" --version > "$build_directory/gsl-version.txt"
sha256sum "$archive" "$reference_directory/driver.cpp" "$build_directory/reference" > "$build_directory/SHA256SUMS"
printf 'Built %s\n' "$build_directory/reference"
