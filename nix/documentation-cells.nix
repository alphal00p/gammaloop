{
  pkgs,
  workspaceRoot,
  notebookWheel,
  docsTypst,
  docsFontPath,
  linnestWasm,
}:
let
  pins = builtins.fromJSON (builtins.readFile ../docs/assets/documentation-runtime-pins.json);
  fetch = name: pkgs.fetchurl { inherit (pins.${name}) url hash; };
  unpack =
    name:
    pkgs.runCommand "documentation-${name}" { } ''
      mkdir -p "$out"
      tar xf ${fetch name} --strip-components=1 -C "$out"
    '';
  source = pkgs.lib.fileset.toSource {
    root = workspaceRoot;
    fileset = pkgs.lib.fileset.unions [
      (workspaceRoot + "/scripts/prepare-documentation-cells.py")
      (workspaceRoot + "/scripts/test_prepare_documentation_cells.py")
      (workspaceRoot + "/scripts/check-documentation-cells-browser.py")
      (workspaceRoot + "/scripts/check-documentation-editors-browser.py")
      (workspaceRoot + "/scripts/check-notebook-startup-browser.py")
      (workspaceRoot + "/docs/documentation-cell-checks.json")
      (workspaceRoot + "/docs/products")
      (workspaceRoot + "/docs/assets")
      (workspaceRoot + "/crates/linnet-py/examples/notebook-cache-worker.js")
      (workspaceRoot + "/docs/api/python")
      (workspaceRoot + "/docs/documentation-cells.toml")
      (workspaceRoot + "/crates/linnest/typst")
      (workspaceRoot + "/crates/kurvst/typst")
      (workspaceRoot + "/vendor/typst-packages")
      (workspaceRoot + "/crates/linnet-render-py/fonts")
      (workspaceRoot + "/crates/linnet-render-py/src/view.js")
      (workspaceRoot + "/crates/feynkit-model/tests/fixtures")
      (workspaceRoot + "/crates/feynkit-py/tests/fixtures")
    ];
  };
  ty = pkgs.runCommand "documentation-ty-72b1fd3ce" { } ''
    mkdir -p "$out"
    cp ${fetch "tyModule"} "$out/ty_wasm.js"
    cp ${fetch "tyWasm"} "$out/ty_wasm_bg.wasm"
    cp ${fetch "tyLicense"} "$out/LICENSE"
  '';
in
pkgs.runCommand "documentation-cell-assets"
  {
    nativeBuildInputs = [ pkgs.python314 ];
    passthru = { inherit source; };
  }
  ''
    mkdir -p "$out/licenses"
    cp ${fetch "typstLicense"} "$out/licenses/typst-ts-Apache-2.0.txt"
    cp ${fetch "robotoLicense"} "$out/licenses/Roboto-OFL.txt"
    cp ${fetch "pyodideLicense"} "$out/licenses/pyodide-MPL-2.0.txt"
    cp ${fetch "symbolicaLicense"} "$out/licenses/Symbolica-License.md"
    cp ${unpack "typstAssets"}/LICENSE "$out/licenses/typst-assets-Apache-2.0.txt"
    cp ${unpack "typstAssets"}/NOTICE "$out/licenses/typst-assets-fonts-NOTICE.txt"
    PYTHONDONTWRITEBYTECODE=1 python ${source}/scripts/test_prepare_documentation_cells.py
    python ${source}/scripts/prepare-documentation-cells.py --root ${source} --output "$out" \
      --compiler ${unpack "compiler"} --renderer ${unpack "renderer"} \
      --pyodide ${unpack "pyodide"} --ty ${ty} --wheel ${notebookWheel}/*.whl \
      --symbolica-stub ${fetch "symbolicaStub"} \
      --typst-packages ${docsTypst}/lib/typst/packages \
      --font-path ${docsFontPath} \
      --default-fonts ${unpack "typstAssets"}/files/fonts --linnest-wasm ${linnestWasm}
  ''
