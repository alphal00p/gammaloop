{
  self,
  pkgs,
  docsPkgs,
  craneLib,
  workspaceRoot,
  cargoSources,
  nonCargoBuildSources,
  commonArgs,
  dummyCargoTarget,
  normalizeWorkspaceHackBuildScriptTimestampScript,
  workspacePackageSrcFor,
  workspaceMissingCargoTargetsScript,
}:
let
  inherit (pkgs) lib;
  docsCargoProfile = "docs";
  documentationRevision = self.dirtyRev or (self.rev or (self.narHash or "local"));
  typst015 =
    assert lib.assertMsg (docsPkgs.typst.version == "0.15.0")
      "the documentation build requires Typst 0.15.0, but the documentation package set provides ${docsPkgs.typst.version}";
    docsPkgs.typst;

  docsTypst = typst015.withPackages (
    typstPackages: with typstPackages; [
      cetz_0_5_1
      mitex_0_2_6
      tidy_0_4_3
    ]
  );

  docsFontPath = "${docsPkgs.roboto}/share/fonts/truetype";

  repositoryDocumentationSources = lib.fileset.fileFilter (
    file:
    let
      name = lib.toLower file.name;
    in
    lib.any (extension: lib.hasSuffix ".${extension}" name) [
      "typ"
      "md"
      "markdown"
      "mdown"
      "mkd"
      "mdx"
      "html"
      "htm"
      "xhtml"
      "shtml"
      "rst"
      "rest"
      "adoc"
      "asciidoc"
      "org"
    ]
    || file.name == "LICENSE"
  ) workspaceRoot;

  documentationDeveloperScopePaths = lib.unique (
    map (scope: scope.path) (
      lib.concatMap (section: lib.concatMap (note: note.scope or [ ]) (section.note or [ ])) (
        (builtins.fromTOML (builtins.readFile (workspaceRoot + "/docs/developers.toml"))).section or [ ]
      )
    )
  );
  # Developer notes verify these source digests during the hermetic site
  # build, including scopes such as flake.nix that are not Cargo or prose.
  documentationDeveloperScopeSources = map (
    path: workspaceRoot + "/${path}"
  ) documentationDeveloperScopePaths;

  documentationSrc = lib.fileset.toSource {
    root = workspaceRoot;
    fileset = lib.fileset.unions (
      documentationDeveloperScopeSources
      ++ [
        cargoSources
        nonCargoBuildSources
        repositoryDocumentationSources
        (workspaceRoot + "/docs")
        (workspaceRoot + "/scripts/check-docs-html.py")
        (workspaceRoot + "/scripts/render-docs-svg-assets.sh")
        (workspaceRoot + "/scripts/update-docs-pages.sh")
        (workspaceRoot + "/examples/api/python")
        (workspaceRoot + "/examples/cli/aa_aa/2L/graphs/GL00.dot")
        (workspaceRoot + "/examples/cli/aa_aa/2L/graphs/GL08.dot")
        (workspaceRoot + "/examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL000.dot")
        (workspaceRoot + "/examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL150.dot")
        (workspaceRoot + "/examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL300.dot")
        (workspaceRoot + "/examples/cli/aa_aa/3L/aa_aa_generation_timings.csv")
        (
          workspaceRoot
          + "/examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_unfiltered_pre_network.sym"
        )
        (workspaceRoot + "/examples/cli/BNL/profiling/bnl_scalar_alias_captures.ansi.txt")
        (workspaceRoot + "/examples/cli/gg_hhh/3L/3L_graph.dot")
        (workspaceRoot + "/tests/resources/graphs/double_triangle.dot")
        (workspaceRoot + "/tests/resources/graphs/gghhh.dot")
        (workspaceRoot + "/tests/resources/graphs/qqx_aaa_pentabox_user_numerator.dot")
        (workspaceRoot + "/tests/resources/graphs/uv_tests/ad_ad_1L_gluon.dot")
        (workspaceRoot + "/tests/resources/graphs/uv_tests/epem_a_bbx.dot")
        (workspaceRoot + "/tests/resources/graphs/epemttbar.dot")
        (workspaceRoot + "/pyproject.toml")
        (workspaceRoot + "/crates/linnet-py/pyproject.toml")
        (workspaceRoot + "/crates/linnet-py/linnet_py.pyi")
        (workspaceRoot + "/crates/linnet-py/tests/test_basic.py")
      ]
    );
  };

  # commonCargoSources includes every TOML file under the workspace. Keep
  # documentation metadata and test fixtures out of this reusable layer.
  documentationWorkspaceBuildSrc = lib.fileset.toSource {
    root = workspaceRoot;
    fileset = lib.fileset.unions [
      (craneLib.fileset.cargoTomlAndLock workspaceRoot)
      (craneLib.fileset.rust workspaceRoot)
      (workspaceRoot + "/.cargo/config.toml")
      (lib.fileset.difference nonCargoBuildSources (
        lib.fileset.unions [
          (workspaceRoot + "/assets/gammalooplogo.svg")
          (workspaceRoot + "/assets/gammalooplogo-dark.svg")
          (workspaceRoot + "/assets/gammalooplogo-light.svg")
        ]
      ))
    ];
  };

  documentationRegistry = builtins.fromTOML (
    builtins.readFile (workspaceRoot + "/docs/products/registry.toml")
  );
  alphal00pDocsCargoTargetRoot = "target/alphal00p-docs-rustdoc";
  alphal00pDocsCargoTarget = "${alphal00pDocsCargoTargetRoot}/cargo-target-v1";
  alphal00pDocsCargoArgs = commonArgs // {
    buildType = docsCargoProfile;
    CARGO_PROFILE = docsCargoProfile;
    CARGO_TARGET_DIR = alphal00pDocsCargoTarget;
    PYO3_PYTHON = "${pkgs.python313}/bin/python3";
    PYTHONPATH = "${pkgs.python313}/lib/python3.13/site-packages";
    # Keep the compile-time Symbolica setting explicit and identical in
    # both the reusable producer and the documentation consumer.
    SYMBOLICA_OEM_LICENSE =
      (builtins.fromTOML (builtins.readFile (workspaceRoot + "/.cargo/config.toml")))
      .env.SYMBOLICA_OEM_LICENSE.value;
  };
  alphal00pDocsRealCargoArgs = alphal00pDocsCargoArgs // {
    postPatch = normalizeWorkspaceHackBuildScriptTimestampScript;
  };
  documentationRustComponents = lib.concatMap (
    product: product.rust_components or [ ]
  ) documentationRegistry.product;
  documentationCatalogFeatures = "gammaloop-reference,vakint-reference";
  documentationPythonExporterBuildCommands = lib.concatMapStringsSep "\n" (product: ''
    cargoWithProfile build --locked -p alphal00p-docs-python-exporter --features ${lib.escapeShellArg product.id}
  '') documentationRegistry.product;
  documentationRustdocBuildCommands = lib.concatMapStringsSep "\n" (
    component:
    let
      features = component.features or [ ];
      featureArgs = lib.optionalString (features != [ ]) (
        " --features ${lib.escapeShellArg (lib.concatStringsSep "," features)}"
      );
    in
    ''
      cargoWithProfile doc --locked --no-deps --no-default-features -p ${lib.escapeShellArg component.package}${featureArgs}
    ''
  ) documentationRustComponents;

  # Keep real workspace outputs keyed only by Rust and build inputs so
  # documentation-only edits reuse the complete Cargo build. End on the
  # first consumer context after the isolated exporter/Rustdoc matrix.
  alphal00pDocsCargoArtifacts = craneLib.mkCargoDerivation (
    alphal00pDocsRealCargoArgs
    // {
      cargoArtifacts = null;
      doInstallCargoArtifacts = true;
      pname = "alphal00p-docs-cargo";
      src = documentationWorkspaceBuildSrc;
      # Prime content-sensitive test dependency contexts without adding
      # authored documentation to this reusable source boundary.
      postPatch = normalizeWorkspaceHackBuildScriptTimestampScript + ''
        install -D -m 0644 ${dummyCargoTarget} crates/alphal00p-docs-examples/build.rs
        install -D -m 0644 ${dummyCargoTarget} crates/alphal00p-docs-examples/src/lib.rs
      '';
      buildPhaseCargoCommand = ''
        cargoWithProfile build --locked -p alphal00p-docs-catalogs --bin alphal00p-docs-catalogs
        ${documentationPythonExporterBuildCommands}
        cargoWithProfile test --locked --no-run -p alphal00p-docs-examples
        cargo clean --profile ${docsCargoProfile} -p alphal00p-docs-examples
        cargoWithProfile build --locked -p alphal00p-docs-builder
        cargoWithProfile build --locked -p linnet-py --features extension-module,abi3-py310
        ${documentationRustdocBuildCommands}
        cargoWithProfile build --locked -p alphal00p-docs-catalogs \
          --features ${lib.escapeShellArg documentationCatalogFeatures} \
          --bin alphal00p-docs-gammaloop-reference \
          --bin alphal00p-docs-vakint-reference
      '';
      checkPhaseCargoCommand = "";
      doCheck = false;
      installPhaseCommand = "";
    }
  );

  # Exercise the embedded Typst watcher behind its opt-in feature without
  # adding Typst's compiler crates to normal workspace or Pages artifacts.
  persistentTypstArgs = commonArgs // {
    pname = "alphal00p-docs-persistent-typst";
    src = workspacePackageSrcFor "alphal00p-docs-builder";
    buildType = "docs-watch";
    CARGO_PROFILE = "docs-watch";
    cargoExtraArgs = "--locked -p alphal00p-docs-builder --features persistent-typst";
    nativeBuildInputs = (commonArgs.nativeBuildInputs or [ ]) ++ [
      docsTypst
      pkgs.roboto
    ];
    TYPST_FONT_PATHS = docsFontPath;
    TYPST_PACKAGE_CACHE_PATH = "${docsTypst}/lib/typst/packages";
  };

  persistentTypstCargoArtifacts = craneLib.buildDepsOnly persistentTypstArgs;

  persistentTypstCheck = craneLib.mkCargoDerivation (
    persistentTypstArgs
    // {
      cargoArtifacts = persistentTypstCargoArtifacts;
      doNotLinkInheritedArtifacts = true;
      doInstallCargoArtifacts = false;
      postPatch = workspaceMissingCargoTargetsScript;
      buildPhaseCargoCommand = ''
        mkdir -p "$out"
        cargoWithProfile check ${persistentTypstArgs.cargoExtraArgs} --all-targets
        cargoWithProfile test ${persistentTypstArgs.cargoExtraArgs} typst_render::persistent::tests
        cargoWithProfile clippy ${persistentTypstArgs.cargoExtraArgs} --all-targets --no-deps -- --deny warnings
      '';
      checkPhaseCargoCommand = "";
      doCheck = false;
      installPhaseCommand = "";
    }
  );

  alphal00pDocsDerivationArgs = alphal00pDocsRealCargoArgs // {
    cargoArtifacts = alphal00pDocsCargoArtifacts;
    src = documentationSrc;
    nativeBuildInputs = (alphal00pDocsCargoArgs.nativeBuildInputs or [ ]) ++ [
      docsTypst
      pkgs.gitMinimal
      pkgs.jujutsu
      pkgs.maturin
      docsPkgs.roboto
      pkgs.uv
    ];
    TYPST_FONT_PATHS = docsFontPath;
    TYPST_PACKAGE_CACHE_PATH = "${docsTypst}/lib/typst/packages";
    ALPHAL00P_DOCS_CARGO_PROFILE = docsCargoProfile;
    doNotLinkInheritedArtifacts = true;
    doInstallCargoArtifacts = false;
    checkPhaseCargoCommand = "";
    doCheck = false;
    installPhaseCommand = "";
  };

  alphal00pDocsPagesChannel =
    let
      requested = builtins.getEnv "ALPHAL00P_DOCS_CHANNEL";
      channel = if requested == "" then "latest" else requested;
    in
    assert lib.assertMsg (builtins.elem channel [
      "latest"
      "snapshot"
    ]) "ALPHAL00P_DOCS_CHANNEL must be latest or snapshot";
    channel;
  alphal00pDocsPagesSnapshotTag =
    let
      tag = builtins.getEnv "ALPHAL00P_DOCS_SNAPSHOT_TAG";
    in
    assert lib.assertMsg (
      (alphal00pDocsPagesChannel == "latest" && tag == "")
      || (alphal00pDocsPagesChannel == "snapshot" && tag != "")
    ) "ALPHAL00P_DOCS_SNAPSHOT_TAG must be set exactly for snapshot builds";
    tag;
  alphal00pDocsPagesGitCommit =
    let
      commit = builtins.getEnv "ALPHAL00P_DOCS_GIT_COMMIT";
    in
    if commit == "" then documentationRevision else commit;
  alphal00pDocsPagesGitTimestamp =
    let
      timestamp = builtins.getEnv "ALPHAL00P_DOCS_GIT_TIMESTAMP";
    in
    if timestamp == "" then toString self.lastModified else timestamp;

  alphal00pDocsValidationCommands = ''
    cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-catalogs --features ${lib.escapeShellArg documentationCatalogFeatures} --bin alphal00p-docs-gammaloop-reference -- --check
    cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-catalogs --features ${lib.escapeShellArg documentationCatalogFeatures} --bin alphal00p-docs-vakint-reference -- --check
    cargo test --locked --profile ${docsCargoProfile} -p alphal00p-docs-examples
    cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features gammaloop -- gammaloop-python docs/api/python/gammaloop-python.pyi --check
    cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features linnet -- linnet-py docs/api/python/linnet-py.pyi --check
    cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features spenso -- spynso3 docs/api/python/spynso3.pyi --check
    cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features idenso -- idenso-community docs/api/python/idenso-community.pyi --check
    cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features vakint -- vakint-community docs/api/python/vakint-community.pyi --check
    cargo test --locked --profile ${docsCargoProfile} -p alphal00p-docs-python-exporter --features gammaloop gammaloop_runtime_surface_and_signatures_match_the_docs_stub
    linnet_python="$TMPDIR/alphal00p-docs-linnet-python"
    export UV_CACHE_DIR="$TMPDIR/alphal00p-docs-uv-cache"
    uv venv "$linnet_python" --python "$PYO3_PYTHON"
    VIRTUAL_ENV="$linnet_python" maturin develop \
      --uv \
      --offline \
      --locked \
      --profile ${docsCargoProfile} \
      --manifest-path crates/linnet-py/Cargo.toml \
      --features extension-module,abi3-py310
    "$linnet_python/bin/python" -m unittest crates/linnet-py/tests/test_basic.py
    cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-builder -- check
    svg_assets="$TMPDIR/alphal00p-svg-assets"
    bash scripts/render-docs-svg-assets.sh "$svg_assets"
    checked_assets=(
      docs/assets/about-*.svg
      docs/assets/graphs/portal-*.svg
      docs/assets/local-unitarity-*.svg
      docs/assets/spensologo.svg
      assets/gammalooplogo*.svg
    )
    generated_assets=(
      "$svg_assets"/docs/assets/about-*.svg
      "$svg_assets"/docs/assets/graphs/portal-*.svg
      "$svg_assets"/docs/assets/local-unitarity-*.svg
      "$svg_assets"/docs/assets/spensologo.svg
      "$svg_assets"/assets/gammalooplogo*.svg
    )
    test "''${#checked_assets[@]}" -eq 32
    test "''${#generated_assets[@]}" -eq 32
    for checked_asset in "''${checked_assets[@]}"; do
      cmp "$checked_asset" "$svg_assets/$checked_asset"
    done
  '';

  alphal00pDocsDeveloperAssertions = output: ''
    test -s "${output}/developers/architecture/documentation-improvement-plan/index.html"
    plan_page="${output}/developers/architecture/documentation-improvement-plan/index.html"
    grep -Fq 'DOC-010' "$plan_page"
    grep -Fq 'id="executive-assessment"' "$plan_page"
    grep -Fq 'href="#executive-assessment"' "$plan_page"
    grep -Fq 'documentation-improvement-plan.typ' "$plan_page"
    test "$(grep -o '<html' "$plan_page" | wc -l)" -eq 1
    test "$(grep -o '<body' "$plan_page" | wc -l)" -eq 1
    grep -Fq '"title": "Executive assessment"' "${output}/developers/search-index.json"
    grep -Fq '"href": "architecture/documentation-improvement-plan/#executive-assessment"' "${output}/developers/search-index.json"
    for product in gammaloop linnet spenso idenso vakint; do
      test -s "${output}/developers/architecture/$product-architecture/index.html"
    done
    test -s "${output}/developers/architecture/spenso-parsing-flow/index.html"
    test ! -e "${output}/developers/architecture/spenso-parsing-flow/diagram.html"
  '';

  # Keep the immutable-route render independently cacheable. The latest
  # package already validates the shared authored/generated inputs, while
  # this fixture exercises the snapshot-specific route layout and HTML.
  # The merge check below can then test publication policy without putting
  # two complete all-product renders in one time-limited NixCI builder.
  alphal00pDocsSnapshotFixture = craneLib.mkCargoDerivation (
    alphal00pDocsDerivationArgs
    // {
      pname = "alphal00p-docs-snapshot-fixture";
      ALPHAL00P_DOCS_GIT_COMMIT = documentationRevision;
      ALPHAL00P_DOCS_GIT_TIMESTAMP = toString self.lastModified;
      buildPhaseCargoCommand = ''
        # Exercise stale global routes and existing product channels in the
        # single cached snapshot render, then remove the fixture-only files.
        mkdir -p "$out/developers" \
          "$out/.staging" \
          "$out/products/gammaloop/latest" \
          "$out/products/gammaloop/snapshots/legacy"
        printf 'removed developer route\n' > "$out/developers/removed-before-rebuild.txt"
        printf 'stale staging\n' > "$out/.staging/incomplete.txt"
        printf 'latest route\n' > "$out/products/gammaloop/latest/.note"
        printf 'product redirect\n' > "$out/products/gammaloop/index.html"
        printf 'historical snapshot\n' > "$out/products/gammaloop/snapshots/legacy/.note"
        cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-builder -- \
          build \
          --product all \
          --channel snapshot \
          --snapshot-tag v0.3.4 \
          --output "$out" \
          --rustdoc-target-root ${lib.escapeShellArg alphal00pDocsCargoTargetRoot}
        test ! -e "$out/developers/removed-before-rebuild.txt"
        test ! -e "$out/.staging"
        grep -Fq 'latest route' "$out/products/gammaloop/latest/.note"
        grep -Fq 'product redirect' "$out/products/gammaloop/index.html"
        grep -Fq 'historical snapshot' "$out/products/gammaloop/snapshots/legacy/.note"
        rm -r "$out/products/gammaloop/latest" \
          "$out/products/gammaloop/snapshots/legacy" \
          "$out/products/gammaloop/index.html"
        python3 scripts/check-docs-html.py "$out"

        for product in gammaloop linnet spenso idenso vakint; do
          test -s "$out/products/$product/snapshots/v0.3.4/.note"
        done
      '';
    }
  );

  alphal00pDocsCheck =
    pkgs.runCommand "alphal00p-docs-check"
      {
        nativeBuildInputs = [
          pkgs.bash
          pkgs.diffutils
          pkgs.findutils
          pkgs.gnugrep
          pkgs.python313
        ];
      }
      ''
        set -x

        docs_first="$TMPDIR/alphal00p-docs-first"
        docs_snapshot=${alphal00pDocsSnapshotFixture}
        mkdir -p "$docs_first"
        cp -R ${alphal00pDocsPages}/. "$docs_first"
        chmod -R u+w "$docs_first"

        docs_pages_test="$TMPDIR/alphal00p-docs-pages-test"
        mkdir -p "$docs_pages_test/products/gammaloop/snapshots/legacy"
        mkdir -p "$docs_pages_test/developers"
        mkdir -p "$docs_pages_test/.git"
        printf 'historical snapshot\n' > "$docs_pages_test/products/gammaloop/snapshots/legacy/.note"
        printf 'removed developer route\n' > "$docs_pages_test/developers/old.txt"
        mv "$docs_first/search-index.json" "$TMPDIR/federated-search-index.json"
        if bash ${
          (workspaceRoot + "/scripts/update-docs-pages.sh")
        } latest "$docs_first" "$docs_pages_test" \
          2> "$TMPDIR/missing-federated-search-index.err"; then
          missing_search_status=0
        else
          missing_search_status=$?
        fi
        mv "$TMPDIR/federated-search-index.json" "$docs_first/search-index.json"
        test "$missing_search_status" -ne 0
        grep -Fq 'latest build has no federated search index' \
          "$TMPDIR/missing-federated-search-index.err"
        bash ${(workspaceRoot + "/scripts/update-docs-pages.sh")} latest "$docs_first" "$docs_pages_test"
        cmp "$docs_pages_test/search-index.json" "$docs_first/search-index.json"
        test -s "$docs_pages_test/products/gammaloop/snapshots/legacy/.note"
        test ! -e "$docs_pages_test/developers/old.txt"
        test -s "$docs_pages_test/developers/architecture/gammaloop-architecture/index.html"
        test -s "$docs_pages_test/developers/architecture/documentation-improvement-plan/index.html"
        cp "$docs_pages_test/index.html" "$TMPDIR/portal-before-snapshot.html"
        cp "$docs_pages_test/developers/.note" "$TMPDIR/developers-before-snapshot.note"
        cp "$docs_pages_test/products/gammaloop/latest/.note" "$TMPDIR/latest-before-snapshot.note"
        bash ${
          (workspaceRoot + "/scripts/update-docs-pages.sh")
        } snapshot "$docs_snapshot" "$docs_pages_test" v0.3.4
        cmp "$docs_pages_test/index.html" "$TMPDIR/portal-before-snapshot.html"
        cmp "$docs_pages_test/developers/.note" "$TMPDIR/developers-before-snapshot.note"
        cmp "$docs_pages_test/products/gammaloop/latest/.note" "$TMPDIR/latest-before-snapshot.note"
        for product in gammaloop linnet spenso idenso vakint; do
          test -s "$docs_pages_test/products/$product/snapshots/v0.3.4/.note"
        done

        mkdir -p "$out"
        cp -R "$docs_first"/. "$out"/

        test -s "$out/index.html"
        test -e "$out/.nojekyll"
        test -s "$out/assets/site.css"
        test -s "$out/assets/site.js"
        test -s "$out/assets/local-unitarity-light.svg"
        test -s "$out/assets/local-unitarity-dark.svg"
        test -s "$out/assets/gammalooplogo-light.svg"
        test -s "$out/assets/gammalooplogo-dark.svg"
        test -s "$out/developers/.note"
        test -s "$out/developers/index.html"
        test -s "$out/developers/search-index.json"
        test -s "$out/developers/assets/site.css"
        test -s "$out/developers/assets/site.js"
        test -s "$out/developers/architecture/gammaloop-architecture/index.html"
        ${alphal00pDocsDeveloperAssertions "$out"}
        for product in gammaloop linnet spenso idenso vakint; do
          product_root="$out/products/$product"
          test -s "$product_root/index.html"
          test -s "$product_root/latest/index.html"
          test -s "$product_root/latest/.note"
          test -s "$product_root/latest/manual.pdf"
          test -s "$product_root/latest/search-index.json"
          test -s "$product_root/latest/snapshot.json"
          test -s "$product_root/latest/tutorial/index.html"
          test -s "$product_root/latest/reference/interfaces/index.html"
          test -s "$product_root/latest/version-history/index.html"
          test -s "$product_root/latest/manual/interfaces/index.html"
          test -s "$product_root/latest/manual/releases/index.html"
          grep -Fq 'url=../../reference/interfaces/' \
            "$product_root/latest/manual/interfaces/index.html"
          grep -Fq 'url=../../version-history/' \
            "$product_root/latest/manual/releases/index.html"
          test -s "$product_root/latest/assets/site.css"
          test -s "$product_root/latest/assets/site.js"
          test -s "$product_root/latest/assets/local-unitarity-light.svg"
          test -s "$product_root/latest/assets/local-unitarity-dark.svg"
          test -s "$product_root/latest/assets/gammalooplogo-light.svg"
          test -s "$product_root/latest/assets/gammalooplogo-dark.svg"
          test -s "$product_root/latest/reference/rust/index.html"
          test -s "$product_root/latest/reference/python/index.html"
          test -n "$(find "$product_root/latest/reference/python" \
            -mindepth 2 -maxdepth 2 -name index.html -print -quit)"
          ! grep -q "Rustdoc generation was skipped" \
            "$product_root/latest/reference/rust/index.html"
        done
        test -s "$out/products/gammaloop/latest/reference/rust/gammalooprs/index.html"
        test -s "$out/products/gammaloop/latest/reference/rust/gammaloop_api/index.html"
        test -s "$out/products/linnet/latest/reference/rust/linnet/index.html"
        test -s "$out/products/spenso/latest/reference/rust/spenso/index.html"
        test -s "$out/products/spenso/latest/reference/rust/spenso_macros/index.html"
        test -s "$out/products/spenso/latest/reference/rust/spenso_hep_lib/index.html"
        test -s "$out/products/idenso/latest/reference/rust/idenso/index.html"
        test -s "$out/products/vakint/latest/reference/rust/vakint/index.html"
        test -s "$out/products/gammaloop/latest/reference/rust/theme.css"
        grep -Fq 'href="../theme.css"' \
          "$out/products/gammaloop/latest/reference/rust/gammalooprs/index.html"
        grep -Fq 'class="alphal00p-rustdoc-bar"' \
          "$out/products/gammaloop/latest/reference/rust/gammalooprs/index.html"
        test -s \
          "$out/products/gammaloop/latest/reference/rust/src/gammalooprs/lib.rs.html"
        grep -Fq 'href="../../theme.css"' \
          "$out/products/gammaloop/latest/reference/rust/src/gammalooprs/lib.rs.html"
        grep -Fq 'class="alphal00p-rustdoc-bar"' \
          "$out/products/gammaloop/latest/reference/rust/src/gammalooprs/lib.rs.html"
        test -s "$out/products/linnet/latest/reference/typst/index.html"
        for typst_page in graph layout drawing physics subgraph; do
          test -s \
            "$out/products/linnet/latest/reference/typst/$typst_page/index.html"
        done
        grep -Fq '"title": "graph.build"' \
          "$out/products/linnet/latest/search-index.json"
        ! grep -Fq '"title": "Parameters"' \
          "$out/products/linnet/latest/search-index.json"
        grep -Fq 'href="../../../../gammaloop/latest/guides/kurvst/"' \
          "$out/products/linnet/latest/reference/typst/index.html"
      '';

  alphal00pDocsPages = craneLib.mkCargoDerivation (
    alphal00pDocsDerivationArgs
    // {
      pname = "alphal00p-docs-pages";
      ALPHAL00P_DOCS_CHANNEL = alphal00pDocsPagesChannel;
      ALPHAL00P_DOCS_GIT_COMMIT = alphal00pDocsPagesGitCommit;
      ALPHAL00P_DOCS_GIT_TIMESTAMP = alphal00pDocsPagesGitTimestamp;
      ALPHAL00P_DOCS_SNAPSHOT_TAG = alphal00pDocsPagesSnapshotTag;
      buildPhaseCargoCommand = ''
        ${alphal00pDocsValidationCommands}

        docs_args=(
          build
          --product all
          --channel ${lib.escapeShellArg alphal00pDocsPagesChannel}
          --output "$out"
          --rustdoc-target-root ${lib.escapeShellArg alphal00pDocsCargoTargetRoot}
        )
        ${lib.optionalString (alphal00pDocsPagesChannel == "snapshot") ''
          docs_args+=(--snapshot-tag ${lib.escapeShellArg alphal00pDocsPagesSnapshotTag})
        ''}
        cargo run --locked --profile ${docsCargoProfile} -p alphal00p-docs-builder -- "''${docs_args[@]}"
        python3 scripts/check-docs-html.py "$out"

        test -s "$out/index.html"
        test -e "$out/.nojekyll"
        ${alphal00pDocsDeveloperAssertions "$out"}
      '';
    }
  );

in
{
  inherit
    docsTypst
    docsFontPath
    documentationDeveloperScopeSources
    alphal00pDocsCargoArtifacts
    alphal00pDocsPages
    alphal00pDocsSnapshotFixture
    alphal00pDocsCheck
    persistentTypstCheck
    ;
}
