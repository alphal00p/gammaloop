# Shared runtime lifecycle, appended to a pure generated runner. No credential
# is evaluated by Nix. This file runs only in the invoking user's process tree.
set +x
set -euo pipefail
umask 077
runtime_package= runtime_filter= runtime_output= runtime_channel= runtime_snapshot_tag=
runtime_reports="$PWD/target/ci-results/$runtime_name"
while (($#)); do
  case "$1" in
    --report-dir|--package|--filter|--output|--channel|--snapshot-tag)
      [[ $# -ge 2 && -n "$2" ]] || { echo "Missing argument for $1" >&2; exit 2; }
      case "$1" in
        --report-dir) runtime_reports=$2 ;;
        --package) runtime_package=$2 ;;
        --filter) runtime_filter=$2 ;;
        --output) runtime_output=$2 ;;
        --channel) runtime_channel=$2 ;;
        --snapshot-tag) runtime_snapshot_tag=$2 ;;
      esac
      shift 2 ;;
    --help)
      echo "$runtime_name: --report-dir PATH; Nextest: --package NAME --filter EXPR; doctest: --package NAME; docs build: --channel latest|snapshot --snapshot-tag TAG --output PATH"
      exit 0 ;;
    *) echo "Unknown runtime argument: $1" >&2; exit 2 ;;
  esac
done
runtime_reports=$(realpath -m "$runtime_reports")
# The docs owner must see symlinks and '..' exactly as supplied. Only anchor
# relative paths to the caller's cwd before entering the private source tree.
case "$runtime_output" in
  ""|/*) ;;
  *) runtime_output="$PWD/$runtime_output" ;;
esac
# Resolve the caller's config before replacing HOME. Preserve an exported key.
if [[ -z "${SYMBOLICA_LICENSE:-}" ]]; then
  if [[ -n "${SYMBOLICA_LICENSE_SIGNED:-}" ]]; then
    SYMBOLICA_LICENSE=$SYMBOLICA_LICENSE_SIGNED
  elif [[ -r "${XDG_CONFIG_HOME:-$HOME/.config}/symbolica/license" ]]; then
    SYMBOLICA_LICENSE=$(< "${XDG_CONFIG_HOME:-$HOME/.config}/symbolica/license")
  fi
fi
if [[ -z "${SYMBOLICA_LICENSE:-}" ]]; then
  echo 'Missing SYMBOLICA_LICENSE environment variable' >&2
  exit 1
fi
export SYMBOLICA_LICENSE
unset SYMBOLICA_LICENSE_SIGNED
mkdir -p "$runtime_reports"
runtime_scratch=$(mktemp -d)
trap 'rm -rf "$runtime_scratch"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
export HOME="$runtime_scratch/home" TMPDIR="$runtime_scratch/tmp"
export XDG_CONFIG_HOME="$HOME/.config" XDG_CACHE_HOME="$HOME/.cache"
mkdir -p "$HOME" "$TMPDIR" "$runtime_scratch/source"
runtime_cores=$(nproc)
if [[ "${NIX_BUILD_CORES:-}" =~ ^[1-9][0-9]*$ ]] && ((NIX_BUILD_CORES < runtime_cores)); then
  runtime_cores=$NIX_BUILD_CORES
fi
((runtime_cores > 0)) || runtime_cores=1
((runtime_cores <= 16)) || runtime_cores=16
# The generated body is a pipeline. Keep its failing exit status for the caller;
# reports are written outside scratch and survive the EXIT cleanup.
set +e
