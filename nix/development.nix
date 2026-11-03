{
  pkgs,
  craneLib,
  workspace,
  workspaceRoot,
  spec ? let
    rawSpec = builtins.getEnv "GAMMALOOP_DEVELOP_SPEC";
  in
    if rawSpec == ""
    then {}
    else builtins.fromJSON rawSpec,
}: let
  inherit (pkgs) lib;
  package = spec.package or null;
  profile = spec.profile or workspace.ciCargoProfile;
  features = lib.sort builtins.lessThan (lib.unique (spec.features or []));
  defaultFeatures = spec.default_features or true;
  pythonName = spec.python or "python313";
  python = pkgs.${pythonName} or (throw "Unsupported Nix Python interpreter: ${pythonName}");
  args = spec.args or [];
  cargoArgs = lib.take (lib.lists.findFirstIndex (arg: arg == "--") (builtins.length args) args) args;
  rustcArgs = lib.drop (builtins.length cargoArgs) args;
  packageDirs = workspace.workspaceMemberPackageDirs;
  hostRoot = "examples/notebooks/symbolica-host";
  hostManifest = "${hostRoot}/Cargo.toml";
  manifest =
    spec.manifest or (
      if package == "alphal00p-notebook-host"
      then hostManifest
      else if package == null
      then "Cargo.toml"
      else "${packageDirs.${package}}/Cargo.toml"
    );
  isHost = manifest == hostManifest;
  relativePath = path:
    builtins.isString path
    && path != ""
    && !lib.hasPrefix "/" path
    && !(builtins.elem ".." (lib.splitString "/" path));
  selectedPackages =
    if package == null
    then builtins.attrNames packageDirs
    else [package];
  # Include optional dependencies as well as the native CI closure: Python
  # feature selection can activate crates absent from that narrower closure.
  sourcePackages = lib.unique (
    lib.concatMap
    (
      selected:
        lib.concatMap (workspace.workspaceDependencyClosureFor workspace.workspaceDependencyNamesFor) (
          workspace.workspaceTestSourcePackageNamesFor selected
        )
    )
    (
      if isHost
      then [
        "feynkit-py"
        "linnet-py"
        "spynso3"
      ]
      else selectedPackages
    )
  );
  sourceFor = extraFilesets:
    workspace.workspacePackageSrcForSourcePackages {
      inherit sourcePackages extraFilesets;
      packageSourcePackages =
        if isHost
        then []
        else selectedPackages;
      testSourcePackages =
        if isHost
        then []
        else selectedPackages;
    };
  rustSrc = sourceFor (lib.optional isHost (workspaceRoot + "/${hostRoot}"));

  featureArgs =
    lib.optional (!defaultFeatures) "--no-default-features"
    ++ lib.optionals (features != []) [
      "--features"
      (lib.concatStringsSep "," features)
    ];
  selectionArgs =
    if isHost
    then [
      "--manifest-path"
      hostManifest
    ]
    else if package == null
    then ["--workspace"]
    else [
      "-p"
      package
    ];
  cargoExtraArgs = lib.escapeShellArgs (["--locked"] ++ selectionArgs ++ featureArgs ++ cargoArgs);
  # These affect reporting/documentation, not dependency compilation, and are
  # not accepted by the check command used to seed metadata for static checks.
  dependencyCargoArgs = lib.escapeShellArgs (
    ["--locked"]
    ++ selectionArgs
    ++ featureArgs
    ++ lib.filter (
      arg:
        !(builtins.elem arg [
          "--no-deps"
          "--document-private-items"
          "--document-hidden-items"
          "--open"
        ])
    )
    cargoArgs
  );
  replacePython = inputs:
    lib.unique (map (input:
      if input == pkgs.python313
      then python
      else input)
    inputs);
  runtimeLibraries = lib.unique [
    python
    pkgs.gmp
    pkgs.mpfr
    pkgs.libmpc
    pkgs.openssl
    pkgs.stdenv.cc.cc.lib
    pkgs.gcc.cc.lib
  ];
  # This independent workspace has its own lockfile. Keep the same exact-rev
  # archive substitution used by the main workspace's vendor expression.
  hostVendor = craneLib.vendorCargoDeps {
    cargoLock = workspaceRoot + "/${hostRoot}/Cargo.lock";
    overrideVendorGitCheckout = packages: drv:
      if lib.any (dependency: dependency.name == "symbolica-typst-atom-payload") packages
      then
        drv.overrideAttrs (_: {
          src = builtins.fetchTree {
            type = "github";
            owner = "symbolica-dev";
            repo = "symbolica-typst-plugin";
            rev = lib.last (lib.splitString "#" (builtins.head packages).source);
          };
        })
      else drv;
  };
  commonArgs =
    workspace.commonArgs
    // {
      pname = "gammaloop-development-${
        if package == null
        then "workspace"
        else package
      }";
      src = rustSrc;
      inherit cargoExtraArgs;
      cargoVendorDir =
        if isHost
        then hostVendor
        else workspace.commonArgs.cargoVendorDir;
      buildType = profile;
      CARGO_PROFILE = profile;
      # Independent workspaces otherwise write beside their own manifest, outside
      # the target tree Crane retains and restores for the real wheel build.
      CARGO_TARGET_DIR = "target";
      CARGO_BUILD_INCREMENTAL = "true";
      CARGO_NET_OFFLINE = "true";
      PYO3_PYTHON = "${python}/bin/python3";
      PYTHON = "${python}/bin/python3";
      nativeBuildInputs = replacePython workspace.commonArgs.nativeBuildInputs;
      buildInputs = replacePython workspace.commonArgs.buildInputs;
      LD_LIBRARY_PATH = lib.makeLibraryPath runtimeLibraries;
      DYLD_LIBRARY_PATH = lib.makeLibraryPath runtimeLibraries;
      postPatch = workspace.workspaceMissingCargoTargetsScript;
      doNotLinkInheritedArtifacts = true;
      NIX_NO_SELF_RPATH = 1;
      doCheck = false;
    };
  # These are seeds, not substituted results. Cargo still checks the exact
  # requested features/targets and recompiles mismatched units (including the
  # CLI's OEM-license setting). No CI feature flags enter the requested command.
  matchingCrateArtifacts =
    !isHost
    && package != null
    && profile == "ci-optim"
    && pythonName == "python313"
    && defaultFeatures
    && features == lib.sort builtins.lessThan (workspace.craneCiFeaturesFor package);
  hostDummyLockfile = lib.optionalString isHost ''
    install -D -m 0644 ${workspaceRoot + "/${hostRoot}/Cargo.lock"} "$out/${hostRoot}/Cargo.lock"
  '';
  dependenciesFor = command: extra:
    workspace.buildDepsOnlyWithArtifacts (
      commonArgs
      // {
        # Crane's dummy source keeps only manifests, lockfiles and target layout.
        # Changes to implementation or Python files do not rebuild this seed.
        pname = "${commonArgs.pname}-${command}";
        extraDummyScript = hostDummyLockfile;
        buildPhaseCargoCommand = "cargoWithProfile ${command} ${dependencyCargoArgs} ${extra}";
      }
    );
  artifactsFor = command:
    if matchingCrateArtifacts
    then workspace.cranePackageOutputs."crate-${package}"
    else dependenciesFor command "";
  rustArgs = command: commonArgs // {cargoArtifacts = artifactsFor command;};
  wrapUfoCli =
    !isHost
    && (package == null || package == "gammaloop-api")
    && (
      builtins.elem "--all-features" cargoArgs
      || lib.any
      (feature: builtins.elem feature features || builtins.elem "gammaloop-api/${feature}" features)
      [
        "ufo_support"
        "python_api"
        "python_stubgen"
        "pyo3-extension-module"
      ]
    );

  pyprojectRoot =
    spec.pyproject or (
      if isHost
      then hostRoot
      else "."
    );
  pyprojectPath = "${pyprojectRoot}/pyproject.toml";
  pyproject = builtins.fromTOML (builtins.readFile (workspaceRoot + "/${pyprojectPath}"));
  pythonSource = pyproject.tool.maturin.python-source or ".";
  pythonModule = builtins.head (lib.splitString "." pyproject.tool.maturin.module-name);
  compatibility = spec.compatibility or "linux";
  # Starting Maturin in the pyproject directory with manifest-path in its
  # configuration handles a Python project below or above its Rust manifest.
  manifestFromPyproject =
    lib.concatStrings (
      map (_: "../") (lib.filter (part: part != "." && part != "") (lib.splitString "/" pyprojectRoot))
    )
    + manifest;
  wheelPyproject = (pkgs.formats.toml {}).generate "development-pyproject.toml" (
    lib.recursiveUpdate pyproject {
      tool.maturin = {
        inherit features;
        no-default-features = !defaultFeatures;
        all-features = false;
        manifest-path = manifestFromPyproject;
      };
    }
  );
  wheelSrc = sourceFor (
    [
      (workspaceRoot + "/${pyprojectPath}")
      (workspaceRoot + "/${pyprojectRoot}/${pythonSource}")
    ]
    ++ lib.optionals isHost [
      (workspaceRoot + "/${hostRoot}")
      # The host's Python tree links to these packages outside its directory.
      (workspaceRoot + "/crates/feynkit-py/python")
      (workspaceRoot + "/crates/linnet-py/python")
    ]
    ++ lib.optional (!isHost) (workspaceRoot + "/README.md")
  );
  wheelArgs =
    commonArgs
    // {
      src = wheelSrc;
      # Maturin enables extension linking for recent PyO3 versions. Match it in
      # the dependency seed; native Rust builds keep their embedding behavior.
      PYO3_BUILD_EXTENSION_MODULE = "1";
    };
  wheelDependencies = workspace.buildDepsOnlyWithArtifacts (
    wheelArgs
    // {
      pname = "${commonArgs.pname}-wheel";
      extraDummyScript = hostDummyLockfile;
      buildPhaseCargoCommand = "cargoWithProfile build ${cargoExtraArgs} --lib";
      postInstall = ''
        test -d "target/${
          if profile == "dev"
          then "debug"
          else profile
        }/deps"
        test -n "$(find "target/${
          if profile == "dev"
          then "debug"
          else profile
        }/deps" -type f -print -quit)"
      '';
    }
  );
in
  assert lib.assertMsg (builtins.isAttrs spec) "GAMMALOOP_DEVELOP_SPEC must be a JSON object";
  assert lib.assertMsg (
    package == null || isHost || builtins.hasAttr package packageDirs
  ) "Unknown development Cargo package";
  assert lib.assertMsg (
    relativePath manifest && relativePath pyprojectRoot
  ) "manifest and pyproject must be workspace-relative paths";
  assert lib.assertMsg (
    builtins.elem profile [
      "dev"
      "dev-optim"
      "ci-optim"
      "release"
    ]
    && !(isHost && profile == "dev-optim")
  ) "Unsupported development Cargo profile";
  assert lib.assertMsg (
    builtins.match "python3[0-9]+" pythonName != null && python ? interpreter
  ) "python must name a supported Nix CPython interpreter (for example python313)";
  assert lib.assertMsg (
    builtins.isBool defaultFeatures
    && builtins.all builtins.isString features
    && builtins.all builtins.isString args
  ) "Invalid development feature or argument settings"; {
    build = assert lib.assertMsg (
      rustcArgs == []
    ) "Rustc arguments after -- are supported only for clippy and wheel";
      craneLib.buildPackage (
        rustArgs "build"
        // {
          # Install binaries/shared libraries normally, and retain a standard Crane
          # archive as the useful output even for packages containing only an rlib.
          doInstallCargoArtifacts = true;
          nativeBuildInputs = commonArgs.nativeBuildInputs ++ lib.optional wrapUfoCli pkgs.makeWrapper;
          postFixup =
            ''test -s "$out/target.tar.zst"''
            + lib.optionalString wrapUfoCli ''

              if [ -x "$out/bin/gammaloop" ]; then
                ${
                if python == pkgs.python313
                then ''
                  # Keep the only public executable relocatable: the dispatcher
                  # copies it into target/PROFILE, so its target must be absolute.
                  mkdir -p "$out/libexec"
                  mv "$out/bin/gammaloop" "$out/libexec/gammaloop"
                  # Pin stdlib discovery to the linked interpreter, ignoring an
                  # inherited PYTHONHOME, and retain explicit Python paths last.
                  makeWrapper "$out/libexec/gammaloop" "$out/bin/gammaloop" \
                    --set PYTHONHOME ${python} \
                    --prefix PYTHONPATH : ${workspace.nextestUfoPython}/${python.sitePackages} \
                    --set-default SYMBOLICA_DATA_DIRECTORY ${workspace.symbolicaRenderData}
                ''
                else ''
                  echo "The GammaLoop UFO executable requires python313 to match its packaged runtime" >&2
                  exit 1
                ''
              }
              fi
            '';
        }
      );
    check = assert lib.assertMsg (rustcArgs == []) "cargo check does not accept rustc arguments after --";
      craneLib.mkCargoDerivation (
        rustArgs "check"
        // {
          pnameSuffix = "-check";
          buildPhaseCargoCommand = "cargoWithProfile check ${cargoExtraArgs}";
        }
      );
    clippy = craneLib.cargoClippy (
      rustArgs "check"
      // {
        cargoClippyExtraArgs = lib.escapeShellArgs rustcArgs;
      }
    );
    doc = assert lib.assertMsg (rustcArgs == []) "cargo doc does not accept rustc arguments after --";
      craneLib.cargoDoc (
        rustArgs "check"
        // {
          cargoDocExtraArgs = "";
        }
      );
    wheel =
      # A flake package must be a derivation even without a dispatcher request.
      # Only the wheel gets this default; native commands still select the workspace.
      if spec == {}
      then
        (import ./development.nix {
          inherit
            pkgs
            craneLib
            workspace
            workspaceRoot
            ;
          spec = {
            package = "gammaloop-api";
            inherit profile;
            default_features = false;
            features = [
              "pyo3-extension-module"
            ];
          };
        }).wheel
      else
        assert lib.assertMsg (
          isHost || (package == "gammaloop-api" && manifest == "crates/gammaloop-api/Cargo.toml")
        ) "Only gammaloop-api and the Symbolica Community host are Python distributions";
        assert lib.assertMsg (
          compatibility == null || compatibility == "linux"
        ) "Development wheels retain Nix-store libraries; only native linux compatibility is supported";
          craneLib.mkCargoDerivation (
            wheelArgs
            // {
              pname = "${commonArgs.pname}-wheel";
              cargoArtifacts = wheelDependencies;
              doInstallCargoArtifacts = false;
              nativeBuildInputs = wheelArgs.nativeBuildInputs ++ [pkgs.maturin];
              buildPhaseCargoCommand = ''
                project_root="$PWD"
                export CARGO_TARGET_DIR="$project_root/target"
                cp ${wheelPyproject} ${lib.escapeShellArg pyprojectPath}
                # Configure only the packaged copy. An explicit caller setting wins;
                # the checkout's Python sources and metadata remain untouched.
                ${python}/bin/python3 - ${lib.escapeShellArg "${pyprojectRoot}/${pythonSource}/${pythonModule}/__init__.py"} <<'PY'
                import ast, pathlib, sys
                path = pathlib.Path(sys.argv[1])
                source = path.read_text()
                offset = 0
                for statement in ast.parse(source).body:
                    if (isinstance(statement, ast.ImportFrom) and statement.module == "__future__") or (
                        offset == 0 and isinstance(statement, ast.Expr)
                        and isinstance(statement.value, ast.Constant) and isinstance(statement.value.value, str)
                    ):
                        offset = statement.end_lineno
                    else:
                        break
                lines = source.splitlines(keepends=True)
                lines.insert(offset, '\nimport os as _nix_os\n_nix_os.environ.setdefault("SYMBOLICA_DATA_DIRECTORY", '
                             + repr(${builtins.toJSON (toString workspace.symbolicaRenderData)}) + ')\ndel _nix_os\n')
                path.write_text("".join(lines))
                PY
                cd ${lib.escapeShellArg pyprojectRoot}
                maturin build --locked --offline \
                  --profile ${lib.escapeShellArg profile} \
                  --interpreter ${python}/bin/python3 \
                  --out "$project_root/wheels" \
                  ${lib.optionalString pkgs.stdenv.hostPlatform.isLinux "--compatibility linux --auditwheel skip"} \
                  ${lib.escapeShellArgs (featureArgs ++ args)}
                cd "$project_root"
              '';
              installPhaseCommand = ''
                mkdir -p "$out/dist" "$out/lib" "$out/nix-support"
                wheels=(wheels/*.whl)
                test "''${#wheels[@]}" -eq 1
                test -s "''${wheels[0]}"
                cp "''${wheels[0]}" "$out/dist/"
                # ZIP compression hides native references from Nix's reference scanner.
                # Retain the native payload beside the wheel so its RPATH dependencies
                # remain live for as long as the installer keeps this output GC-rooted.
                ${python}/bin/python3 - "$out" <<'PY'
                import pathlib, sys, zipfile
                output = pathlib.Path(sys.argv[1])
                wheel = next((output / "dist").glob("*.whl"))
                libraries = 0
                with zipfile.ZipFile(wheel) as archive:
                    for name in archive.namelist():
                        if name.endswith((".so", ".dylib")):
                            relative = pathlib.PurePosixPath(name)
                            assert not relative.is_absolute() and ".." not in relative.parts
                            destination = output / "lib" / relative
                            destination.parent.mkdir(parents=True, exist_ok=True)
                            destination.write_bytes(archive.read(name))
                            libraries += 1
                assert libraries, "wheel contains no native extension"
                PY
                printf '%s\n' ${workspace.symbolicaRenderData} > "$out/nix-support/symbolica-data-directory"
                printf '%s\n' ${lib.escapeShellArgs (map toString runtimeLibraries)} > "$out/nix-support/runtime-libraries"
              '';
              # Preserve the extracted bytes: they document the exact packaged payload.
              dontStrip = true;
              dontPatchELF = true;
            }
          );
  }
