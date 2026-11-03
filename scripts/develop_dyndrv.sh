# shellcheck shell=bash
# shellcheck disable=SC2154  # Shared parser state is owned by develop.sh.
# The public Cargo argv is constructed once by develop.sh. Only worker ownership
# and immutable compiler delivery belong to this backend.
dyndrv_main() {
    case "$action" in
        upload-preview)
            fail "dyndrv does not support $action; select backend=nix or backend=cargo explicitly" ;;
        fmt|clean) ;;
    esac
    local bootstrap status=0
    mkdir -p -- "$root/target/dyndrv"
    dyndrv_bootstrap_root=$(mktemp -d "$root/target/dyndrv/.bootstrap.XXXXXXXX")
    trap 'rm -rf -- "$dyndrv_bootstrap_root"' EXIT
    bootstrap=$(
        unset SYMBOLICA_LICENSE SYMBOLICA_LICENSE_PATH SYMBOLICA_LICENSE_SIGNED SYMBOLICA_LICENSE_FILE
        nix build --impure --out-link "$dyndrv_bootstrap_root/root" --print-out-paths --expr '
          (import ((builtins.getEnv "GAMMALOOP_WORKSPACE_ROOT") + "/nix/dyndrv-bootstrap.nix")).tools
        '
    )
    if [[ $action == python ]]; then
        PATH="$bootstrap/tools/bin:$PATH" GAMMALOOP_DYNDRV_BOOTSTRAP="$bootstrap" "$bootstrap/tools/bin/python3" \
            "$root/scripts/python_wheel.py" --backend dyndrv "$@" || status=$?
    else
        "$bootstrap/tools/bin/cargo-dyndrv-worker" \
            --bootstrap "$bootstrap" --workspace "$root" --action "$action" -- "$@" || status=$?
    fi
    rm -rf -- "$dyndrv_bootstrap_root"
    trap - EXIT
    return "$status"
}
