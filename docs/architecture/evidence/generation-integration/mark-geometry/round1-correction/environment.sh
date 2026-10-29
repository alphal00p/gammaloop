if [ -z "${SYMBOLICA_LICENSE:-}" ]; then
  mark_license_file="${XDG_CONFIG_HOME:-$HOME/.config}/symbolica/license"
  if [ -f "$mark_license_file" ]; then export SYMBOLICA_LICENSE="$(cat "$mark_license_file")"; fi
fi
export CARGO_BUILD_JOBS=2 NEXTEST_TEST_THREADS=2 RAYON_NUM_THREADS=2 INSTA_UPDATE=no
unset INSTA_FORCE_PASS
export NIX_CONFIG="${NIX_CONFIG:-}
cores = 2
max-jobs = 1"
