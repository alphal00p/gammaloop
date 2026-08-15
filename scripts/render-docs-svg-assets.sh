#!/usr/bin/env bash
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
invocation_dir=$(pwd -P)
if [[ -n ${1-} ]]; then
    case $1 in
        /*) output_root=$1 ;;
        *) output_root="$invocation_dir/$1" ;;
    esac
else
    output_root=$repo_root
fi

mkdir -p "$output_root"
output_root=$(CDPATH= cd -- "$output_root" && pwd -P)
cd "$repo_root"

case $(typst --version) in
    "typst 0.15.0"*) ;;
    *)
        printf 'website SVG rendering requires Typst 0.15.0\n' >&2
        exit 1
        ;;
esac

render() {
    local source=$1
    local destination=$2
    local theme=${3-light}
    mkdir -p "$(dirname -- "$destination")"
    typst compile \
        --root "$repo_root" \
        --format svg \
        --creation-timestamp 0 \
        --input "theme=$theme" \
        "$repo_root/$source" \
        "$destination"
    grep -Fq '<svg' "$destination"
}

wait_batch() {
    local status=0
    local pid
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            status=1
        fi
    done
    (( status == 0 )) || exit "$status"
    pids=()
}

queue_render() {
    render "$@" &
    pids+=("$!")
    if (( ${#pids[@]} == 4 )); then
        wait_batch
    fi
}

pids=()
graph_sources=(docs/assets/typst/portal-graphs/graphs/portal-graph-*.typ)
test "${#graph_sources[@]}" -eq 11
for source in "${graph_sources[@]}"; do
    name=$(basename -- "$source" .typ)
    for theme in light dark; do
        queue_render \
            "$source" \
            "$output_root/docs/assets/graphs/$name-$theme.svg" \
            "$theme"
    done
done

queue_render \
    docs/assets/typst/marks/local-unitarity.typ \
    "$output_root/docs/assets/local-unitarity-light.svg" \
    light
queue_render \
    docs/assets/typst/marks/local-unitarity.typ \
    "$output_root/docs/assets/local-unitarity-dark.svg" \
    dark
queue_render \
    docs/assets/typst/marks/gammaloop.typ \
    "$output_root/assets/gammalooplogo-light.svg" \
    light
queue_render \
    docs/assets/typst/marks/gammaloop.typ \
    "$output_root/assets/gammalooplogo-dark.svg" \
    dark
queue_render \
    docs/assets/typst/marks/gammaloop.typ \
    "$output_root/assets/gammalooplogo.svg" \
    light
queue_render \
    docs/assets/typst/marks/spenso.typ \
    "$output_root/docs/assets/spensologo.svg" \
    light
wait_batch

generated=(
    "$output_root"/docs/assets/graphs/portal-graph-*.svg
    "$output_root"/docs/assets/local-unitarity-*.svg
    "$output_root"/docs/assets/spensologo.svg
    "$output_root"/assets/gammalooplogo*.svg
)
test "${#generated[@]}" -eq 28
