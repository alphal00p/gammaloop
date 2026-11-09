{
  pkgs,
  alphal00pDocsNotebooks,
  notebookWasmWheel,
  notebookUfoWheel,
}:
pkgs.writeShellApplication {
  name = "test-notebook-examples";
  runtimeInputs = [
    pkgs.coreutils
    pkgs.findutils
  ];
  text = ''
    if [[ $# -ne 2 ]]; then
      echo "Usage: test-notebook-examples REPORT_DIRECTORY BUILT_SITE_DIRECTORY" >&2
      exit 2
    fi
    report=$(realpath -m -- "$1")
    site=$(realpath -e -- "$2")
    test -s "$site/index.html"
    mkdir -p "$report"
    mapfile -t wheels < <(find -L ${notebookWasmWheel} -maxdepth 2 -name '*.whl' -type f)
    if [[ ''${#wheels[@]} -ne 1 ]]; then
      echo "Expected one combined browser wheel" >&2
      exit 1
    fi
    for product in linnet gammaloop spenso idenso feynkit; do
      dependencies=()
      cache_checks=()
      if [[ "$product" = feynkit ]]; then
        dependencies=(--dependency-wheel ${notebookUfoWheel}/${notebookUfoWheel.filename})
        cache_checks=(--cache-report "$report/cache-flow.json")
      fi
      ${alphal00pDocsNotebooks}/bin/alphal00p-docs-notebooks \
        --docs "$product" --wheel "''${wheels[0]}" "''${dependencies[@]}" \
        --output "$site/products/$product/latest/assets/notebooks" \
        "''${cache_checks[@]}" --browser-smoke --browser-executable ${pkgs.chromium}/bin/chromium \
        --timeout 300 2>&1 | tee "$report/$product.log"
    done
    printf 'All browser notebook examples passed.\n' | tee "$report/result.txt"
  '';
}
