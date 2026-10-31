#!/usr/bin/env bash
# The host path deliberately needs only Bash and Cargo.
set -euo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$root"
export GAMMALOOP_WORKSPACE_ROOT="$root"
fail() { echo "error: $*" >&2; exit 2; }
backend=auto
case "${1-}" in
    --backend) [[ $# -ge 2 && -n $2 ]] || fail 'missing backend'; backend=$2; shift 2 ;;
    --backend=*) backend=${1#*=}; shift ;;
esac
case "${1-}" in --backend|--backend=*) fail 'backend may be given only once before the action' ;; esac
case "$backend" in auto|nix|cargo|dyndrv) ;; *) fail "invalid backend: $backend" ;; esac
if [[ $backend == auto ]]; then
    if command -v nix >/dev/null 2>&1; then backend=nix; else backend=cargo; fi
fi
[[ $backend != nix && $backend != dyndrv ]] || command -v nix >/dev/null 2>&1 || fail 'Nix was explicitly selected but is not installed'
action=${1-}; [[ -n $action ]] || fail 'missing action'; shift
case "$action" in test|list|build|check|clippy|doc|run|fmt|clean|python|upload-preview|stop|gc) ;; *) fail "unknown action: $action" ;; esac
case "$action" in stop|gc) [[ $backend == dyndrv ]] || fail "$action requires backend=dyndrv" ;; esac
if [[ -z ${GAMMALOOP_BUILD_PROFILE:-} ]]; then
    GAMMALOOP_BUILD_PROFILE=$(< "$root/.config/build-profile")
fi
export GAMMALOOP_BUILD_PROFILE
echo "Backend: $backend" >&2
load_license() {
    if [[ -z ${SYMBOLICA_LICENSE:-} ]]; then
        if [[ -n ${SYMBOLICA_LICENSE_SIGNED:-} ]]; then
            export SYMBOLICA_LICENSE=$SYMBOLICA_LICENSE_SIGNED
        else
            local license_file=${XDG_CONFIG_HOME:-${HOME:-}/.config}/symbolica/license
            if [[ -f $license_file ]]; then SYMBOLICA_LICENSE=$(cat "$license_file"); export SYMBOLICA_LICENSE; fi
        fi
    fi
}
if [[ $backend == nix ]]; then
    source "$root/scripts/develop_nix.sh"
    nix_main "$@"
    exit
fi
if [[ $backend == dyndrv ]]; then
    source "$root/scripts/develop_dyndrv.sh"
    # Leave native command construction below shared with ordinary Cargo.
    cargo_execute() { dyndrv_main "$@"; }
else
    cargo_execute() { exec cargo "$@"; }
fi
case "$action" in
    python) load_license; if [[ $backend == dyndrv ]]; then dyndrv_main "$@"; exit; fi; exec python3 "$root/scripts/python_wheel.py" --backend cargo "$@" ;;
    upload-preview) [[ $backend != dyndrv ]] || dyndrv_main "$@"; fail "$action requires Nix" ;;
    stop|gc) dyndrv_main "$@"; exit ;;
    fmt) cargo_execute fmt --all "$@"; exit ;;
    clean) cargo_execute clean "$@"; exit ;;
esac
load_license
args=(); profile_set=false; selection=false; runtime_profile=ci_gammaloop; no_tests=fail
if [[ $action == test || $action == list ]]; then
    profile=$GAMMALOOP_BUILD_PROFILE; release=false; workspace=false; packages=false; filters=(); forwarded=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --) forwarded+=("$@"); break ;;
            --cargo-profile) [[ $# -ge 2 ]] || fail 'missing Cargo profile'; profile=$2; shift 2 ;;
            --cargo-profile=*) profile=${1#*=}; shift ;;
            --release) release=true; shift ;;
            -P|--profile) [[ $# -ge 2 ]] || fail 'missing Nextest profile'; runtime_profile=$2; shift 2 ;;
            --profile=*|-P=*) runtime_profile=${1#*=}; shift ;;
            -P?*) runtime_profile=${1#-P}; shift ;;
            -E|--filterset|--filter-expr) [[ $# -ge 2 ]] || fail 'missing filter'; filters+=("$2"); shift 2 ;;
            --filterset=*|--filter-expr=*|-E=*) filters+=("${1#*=}"); shift ;;
            -E?*) filters+=("${1#-E}"); shift ;;
            --no-tests) [[ $# -ge 2 ]] || fail 'missing no-tests value'; no_tests=$2; shift 2 ;;
            --no-tests=*) no_tests=${1#*=}; shift ;;
            -p|--package) [[ $# -ge 2 ]] || fail 'missing package'; packages=true; selection=true; args+=("$1" "$2"); shift 2 ;;
            --workspace) workspace=true; selection=true; args+=("$1"); shift ;;
            --package=*|-p?*) packages=true; selection=true; args+=("$1"); shift ;;
            *) forwarded+=("$1"); shift ;;
        esac
    done
    "$selection" || fail 'select a crate with --package, or use --workspace'
    [[ $workspace == false || $packages == false ]] || fail 'select packages or the workspace, not both'
    "$release" && profile=release
    expression='all()'; union=''
    for item in ${filters[@]+"${filters[@]}"}; do union="${union:+$union | }($item)"; done
    [[ -z $union ]] || expression="(all()) & ($union)"
    export NO_SYMBOLICA_OEM_LICENSE=1
    if [[ $runtime_profile == ci_slow ]]; then
        export RUST_MIN_STACK=134217728 RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-${NEXTEST_TEST_THREADS:-16}} GAMMALOOP_TESTS_NO_CLEAN_STATE=1
    fi
    command=(cargo nextest run --locked --cargo-profile "$profile" ${args[@]+"${args[@]}"} --profile "$runtime_profile")
    if [[ $action == list ]]; then command[2]=list; command+=(--message-format json); else command+=(--no-tests "$no_tests"); fi
    if [[ $backend == dyndrv ]]; then
        dyndrv_main "${command[@]:1}" -E "$expression" ${forwarded[@]+"${forwarded[@]}"}
    else
        exec "${command[@]}" -E "$expression" ${forwarded[@]+"${forwarded[@]}"}
    fi
    exit
fi
# Cargo validates native flags. Inspect only the build prefix to supply defaults.
for arg in "$@"; do
    [[ $arg != -- ]] || break
    case "$arg" in --profile|--profile=*|--release|-r) profile_set=true ;; -p|--package|--package=*|-p?*|--workspace) selection=true ;; esac
done
command=(cargo "$action" --locked)
"$profile_set" || command+=(--profile "$GAMMALOOP_BUILD_PROFILE")
"$selection" || command+=(-p gammaloop-api)
if [[ $backend == dyndrv ]]; then
    dyndrv_main "${command[@]:1}" "$@"
else
    exec "${command[@]}" "$@"
fi
