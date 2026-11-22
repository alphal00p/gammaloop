# Marimo and its plotting stack are fixed by the runner. The extension and
# model loader come from the measured revision, including on the first PR
# introducing these gates (whose base need not contain this expression).
{ runner, subject, system ? builtins.currentSystem }:
let
  harness = builtins.getFlake runner;
  measured = (builtins.getFlake subject).packages.${system};
  pkgs = harness.inputs.nixpkgs.legacyPackages.${system};
  docsPkgs = harness.inputs.nixpkgs-docs.legacyPackages.${system};
  python = harness.packages.${system}.notebook-performance-dependencies;
  browsers = docsPkgs.playwright-driver.browsers.override {
    withChromium = false;
    withFirefox = false;
    withWebkit = false;
    withFfmpeg = false;
    withChromiumHeadlessShell = true;
  };
  spec = builtins.fromJSON (builtins.getEnv "GAMMALOOP_DEVELOP_SPEC");
in
assert pkgs.lib.assertMsg (spec == { package = "alphal00p-notebook-host"; profile = "ci-optim"; })
  "Set GAMMALOOP_DEVELOP_SPEC to the ci-optim alphal00p-notebook-host specification";
pkgs.runCommand "notebook-performance-python" { nativeBuildInputs = [ pkgs.makeWrapper ]; } ''
  mkdir -p "$out/bin" "$out/site-packages"
  for wheel in ${measured.development-wheel}/dist/*.whl ${measured.notebook-ufo-wheel}/*.whl; do
    ${python}/bin/python -m zipfile -e "$wheel" "$out/site-packages"
  done
  makeWrapper ${python}/bin/python "$out/bin/python" \
    --set PYTHONPATH "$out/site-packages" \
    --set PYTHONNOUSERSITE 1 \
    --set MPLBACKEND Agg \
    --set PLAYWRIGHT_BROWSERS_PATH ${browsers} \
    --set FONTCONFIG_FILE ${docsPkgs.makeFontsConf { fontDirectories = [ docsPkgs.roboto docsPkgs.dejavu_fonts ]; }} \
    --prefix PATH : ${pkgs.lib.makeBinPath [ pkgs.stdenv.cc ]}
''
