# shellcheck shell=bash
# Sourced by develop.sh, which supplies root, action and load_license.
# shellcheck disable=SC2154
nix_safe() ( unset SYMBOLICA_LICENSE SYMBOLICA_LICENSE_PATH SYMBOLICA_LICENSE_SIGNED SYMBOLICA_LICENSE_FILE GAMMALOOP_NIX_TEST_ARGS; "$@"; )
nix_cleanup() {
    if [[ -n ${pending_root:-} ]]; then
        if "$retain_pending"; then echo "Retained binary build root: $pending_root" >&2; else rm -f -- "$pending_root"; fi
    fi
    [[ -z ${install_temporary:-} ]] || rm -f -- "$install_temporary"
    rm -rf -- "$temporary"
}
nix_build() {
    nix_safe nix build --impure --print-out-paths --no-write-lock-file --accept-flake-config --out-link "$2" "$flake#$1"
}
nix_retain() {
    mkdir -p -- "$(dirname -- "$2")"
    nix_safe nix-store --realise "$1" --add-root "$2" --indirect >/dev/null
}
nix_snapshot() {
    local name parent cursor component
    mkdir "$temporary/source"
    if [[ -d $root/.jj ]]; then
        jj file list --no-pager --color never -T 'path ++ "\0"' > "$temporary/files"
    elif [[ -e $root/.git ]]; then
        git ls-files --cached -z > "$temporary/files"
    else fail 'cannot list tracked source without JJ or Git'; fi
    while IFS= read -r -d '' name; do
        case "$name" in ''|target|target/*|.venv|.venv/*) continue ;; /*|..|../*|*/../*|*/..) fail "tracked path escapes workspace: $name" ;; esac
        parent=${name%/*}
        if [[ $parent != "$name" ]]; then
            cursor=$root
            while [[ -n $parent ]]; do
                component=${parent%%/*}; cursor=$cursor/$component
                [[ ! -L $cursor ]] || fail "tracked path traverses directory symlink: $name"
                if [[ $parent == */* ]]; then parent=${parent#*/}; else parent=''; fi
            done
        fi
        [[ -e $root/$name || -L $root/$name ]] || continue
        [[ ! -d $root/$name || -L $root/$name ]] || fail "tracked path is a directory: $name"
        printf '%s\0' "$name"
    done < "$temporary/files" > "$temporary/selected-files"
    tar -C "$root" --null --verbatim-files-from --no-recursion \
        -T "$temporary/selected-files" -cf - | tar -xf - -C "$temporary/source"
    export GAMMALOOP_FLAKE_SOURCE=$temporary/source
}
nix_main() {
    temporary=$(mktemp -d "${TMPDIR:-/tmp}/gammaloop-development.XXXXXXXX")
    pending_root=''; retain_pending=false; install_temporary=''
    trap nix_cleanup EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM
    # Bootstrap from the lock and the small tool definition before copying source.
    local tools
    # Interpolation in the expression belongs to Nix.
    # shellcheck disable=SC2016
    tools=$(nix_safe nix build --impure --out-link "$temporary/tools" --print-out-paths --expr '
        let
          root = builtins.getEnv "GAMMALOOP_WORKSPACE_ROOT";
          lock = builtins.fromJSON (builtins.readFile (root + "/flake.lock"));
          locked = lock.nodes.${lock.nodes.${lock.root}.inputs.nixpkgs}.locked;
          pkgs = import (builtins.fetchTree locked).outPath {
            system = builtins.currentSystem; config = {}; overlays = [];
          };
        in import (root + "/nix/development-tools.nix") { inherit pkgs; }
    ')
    export PATH="$tools/bin:$PATH"
    nix_snapshot
    local flake
    flake="path:$(jq -nr --arg path "$GAMMALOOP_FLAKE_SOURCE" '$path | split("/") | map(@uri) | join("/")')"
    case "$action" in
        fmt|clean)
            local cmd=(cargo "$action"); [[ $action != fmt ]] || cmd+=(--all)
            nix_safe nix develop --no-write-lock-file "$flake" -c "${cmd[@]}" "$@"; return ;;
        python|upload-preview)
            local python_tools
            python_tools=$(nix_build development-python-tools "$temporary/python")
            export PATH="$python_tools/bin:$PATH"
            load_license
            if [[ $action == python ]]; then
                "$python_tools/bin/python3" "$root/scripts/python_wheel.py" --backend nix "$@"
            else
                "$python_tools/bin/python3" "$root/scripts/cache_publication.py" "$@"
            fi
            return ;;
        test|list) nix_test "$@"; return ;;
    esac
    nix_rust "$@"
}
nix_test() {
    local packages=() features=() forwarded=() workspace=false profile=$GAMMALOOP_BUILD_PROFILE release=false output reports
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --) forwarded+=("$@"); break ;;
            -p|--package) [[ $# -ge 2 ]] || fail 'missing package'; packages+=("$2"); shift 2 ;;
            --package=*) packages+=("${1#*=}"); shift ;;
            -p?*) packages+=("${1#-p}"); shift ;;
            --workspace) workspace=true; shift ;;
            --cargo-profile) [[ $# -ge 2 ]] || fail 'missing Cargo profile'; profile=$2; shift 2 ;;
            --cargo-profile=*) profile=${1#*=}; shift ;;
            --release) release=true; shift ;;
            -F|--features) [[ $# -ge 2 ]] || fail 'missing features'; features+=("$2"); shift 2 ;;
            --features=*) features+=("${1#*=}"); shift ;;
            -F?*) features+=("${1#-F}"); shift ;;
            *) forwarded+=("$1"); shift ;;
        esac
    done
    if "$workspace"; then
        [[ -z ${packages[*]-} ]] || fail 'select packages or the workspace, not both'
        unset GAMMALOOP_NIX_TEST_PACKAGES
    else
        [[ -n ${packages[*]-} ]] || fail 'select a crate with --package, or use --workspace'
        GAMMALOOP_NIX_TEST_PACKAGES=$(jq -cn --args '$ARGS.positional | unique' -- ${packages[@]+"${packages[@]}"})
        export GAMMALOOP_NIX_TEST_PACKAGES
    fi
    "$release" && profile=release
    export GAMMALOOP_NIX_TEST_PROFILE=$profile
    GAMMALOOP_NIX_TEST_FEATURES=$(jq -cn --args '$ARGS.positional' -- ${features[@]+"${features[@]}"})
    export GAMMALOOP_NIX_TEST_FEATURES
    export GAMMALOOP_NIX_TEST_LIST=0
    [[ $action != list ]] || GAMMALOOP_NIX_TEST_LIST=1
    mkdir -p "$root/target/nix" "$root/target/test-reports"
    pending_root=$root/target/nix/.pending-tests-${temporary##*.}
    output=$(nix_build test-selected "$pending_root")
    reports=$(mktemp -d "$root/target/test-reports/nix-XXXXXXXX")
    echo "Test reports: $reports" >&2
    load_license
    "$output/bin/nix-test-selected" "$reports" ${forwarded[@]+"${forwarded[@]}"}
}
nix_rust() {
    local package='' workspace=false profile=$GAMMALOOP_BUILD_PROFILE release=false features=() args=() run_args=() defaults=true open=false
    local output attribute permanent build_root directory binary destination executable='' count=0
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --) if [[ $action == run ]]; then shift; run_args=("$@"); else args+=("$@"); fi; break ;;
            -p|--package) [[ $# -ge 2 ]] || fail 'missing package'; [[ -z $package ]] || fail 'select one package'; package=$2; shift 2 ;;
            --package=*|-p?*) [[ -z $package ]] || fail 'select one package'; case "$1" in --*) package=${1#*=} ;; *) package=${1#-p} ;; esac; shift ;;
            --workspace) workspace=true; shift ;;
            --profile) [[ $# -ge 2 ]] || fail 'missing profile'; profile=$2; shift 2 ;;
            --profile=*) profile=${1#*=}; shift ;;
            --release|-r) release=true; shift ;;
            -F|--features) [[ $# -ge 2 ]] || fail 'missing features'; features+=("$2"); shift 2 ;;
            --features=*) features+=("${1#*=}"); shift ;;
            -F?*) features+=("${1#-F}"); shift ;;
            --no-default-features) defaults=false; shift ;;
            --open) open=true; shift ;;
            *) args+=("$1"); shift ;;
        esac
    done
    "$release" && profile=release
    if "$workspace"; then [[ -z $package ]] || fail 'select package or workspace'; else package=${package:-gammaloop-api}; fi
    [[ $open == false || $action == doc ]] || fail '--open requires doc'
    export GAMMALOOP_DEVELOP_SPEC
    GAMMALOOP_DEVELOP_SPEC=$(jq -cn --arg package "$package" --arg profile "$profile" --argjson defaults "$defaults" --argjson features "$(jq -cn --args '$ARGS.positional | map(split(",")) | add // [] | map(select(length>0)) | unique' -- ${features[@]+"${features[@]}"})" --args '{package:(if $package=="" then null else $package end),profile:$profile,features:$features,default_features:$defaults,args:$ARGS.positional}' -- ${args[@]+"${args[@]}"})
    mkdir -p "$root/target/nix"
    permanent=$root/target/nix/$action-${package:-workspace}-$profile
    attribute=$action; [[ $action != run ]] || attribute=build
    build_root=$permanent
    if [[ $action == build || $action == run ]]; then
        pending_root=$root/target/nix/.pending-$action-${temporary##*.}
        build_root=$pending_root
    fi
    output=$(nix_build "development-$attribute" "$build_root")
    if [[ $action == build || $action == run ]]; then
        directory=$root/target/$profile; [[ $profile != dev ]] || directory=$root/target/debug
        mkdir -p "$directory"
        exec 9> "$root/target/nix/.install-$profile.lock"
        flock -x 9
        for binary in "$output"/bin/*; do
            [[ -f $binary ]] || continue
            destination=$directory/${binary##*/}
            install_temporary=$(mktemp "$directory/.${binary##*/}.XXXXXXXX")
            cp -- "$binary" "$install_temporary"; chmod --reference="$binary" "$install_temporary"; chmod u+w "$install_temporary"
            mv -f -- "$install_temporary" "$destination"; install_temporary=''
            retain_pending=true
            nix_retain "$output" "$root/target/nix/bin/$profile/${binary##*/}"
            retain_pending=false
        done
        nix_retain "$output" "$permanent"
        flock -u 9; exec 9>&-
    fi
    local label
    case "$action" in build) label=Build ;; check) label=Check ;; clippy) label=Clippy ;; doc) label=Doc ;; run) label=Run ;; esac
    echo "$label output: $output"
    if "$open"; then
        local doc_package=${package:-gammaloop-api}
        local index=$output/share/doc/${doc_package//-/_}/index.html
        if command -v xdg-open >/dev/null 2>&1; then xdg-open "$index"; else open "$index"; fi
    fi
    if [[ $action == run ]]; then
        for binary in "$output"/bin/*; do
            if [[ -f $binary ]]; then count=$((count+1)); executable=$binary; fi
        done
        [[ $count -eq 1 ]] || fail 'run requires exactly one binary; select it with --bin'
        load_license
        "$executable" ${run_args[@]+"${run_args[@]}"}
    fi
}
