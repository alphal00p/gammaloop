use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

use eyre::{Context, ContextCompat, Result, ensure};
use tempfile::Builder as TempDirBuilder;

/// A complete Typst bundle compilation request.
///
/// Both the release builder and live watcher consume this representation. The
/// former executes the pinned Typst CLI while the latter retains an in-process
/// world for each stable `key`.
pub(crate) struct TypstRenderJob {
    pub key: String,
    pub label: String,
    pub wrapper: String,
    pub inputs: Vec<(String, String)>,
    pub timestamp: i64,
    pub output: PathBuf,
    pub dependency_file: Option<PathBuf>,
}

pub(crate) trait TypstRenderer {
    fn render(&mut self, job: TypstRenderJob) -> Result<()>;

    /// Advance process-global caches once after a complete site generation.
    fn finish_generation(&mut self) {}
}

pub(crate) struct CliTypstRenderer {
    root: PathBuf,
}

impl CliTypstRenderer {
    pub(crate) fn new(root: &Path) -> Self {
        Self {
            root: root.to_path_buf(),
        }
    }
}

impl TypstRenderer for CliTypstRenderer {
    fn render(&mut self, job: TypstRenderJob) -> Result<()> {
        let target = self.root.join("target");
        fs::create_dir_all(&target)?;
        let work = TempDirBuilder::new()
            .prefix(&format!("alphal00p-typst-cli-{}-", job.key))
            .tempdir_in(target)?;
        let wrapper = work.path().join("main.typ");
        fs::write(&wrapper, job.wrapper)?;

        let mut command = Command::new("typst");
        command.current_dir(&self.root).args([
            "compile",
            "--root",
            self.root.to_str().context("workspace root is not UTF-8")?,
            "--features",
            "html,bundle",
            "--format",
            "bundle",
            "--creation-timestamp",
            &job.timestamp.to_string(),
        ]);
        for (key, value) in job.inputs {
            command.args(["--input", &format!("{key}={value}")]);
        }
        if let Some(path) = &job.dependency_file {
            if let Some(parent) = path.parent() {
                fs::create_dir_all(parent)?;
            }
            command.args(["--deps-format", "zero", "--deps"]).arg(path);
        }
        let status = command
            .arg(wrapper)
            .arg(&job.output)
            .status()
            .wrap_err("failed to launch Typst 0.15; enter the Nix development shell")?;
        ensure!(status.success(), "Typst failed for {}", job.label);
        Ok(())
    }
}

#[cfg(feature = "persistent-typst")]
mod persistent {
    use std::{
        collections::{BTreeSet, HashMap, hash_map::Entry},
        env,
        sync::Arc,
    };

    use chrono::{Datelike, TimeZone, Timelike, Utc};
    use comemo::evict;
    use eyre::{ContextCompat, bail, eyre};
    use typst::{
        Feature, Library, LibraryExt, World,
        diag::{FileResult, SourceDiagnostic, Warned},
        foundations::{Bytes, Datetime, Dict, Duration, IntoValue},
        syntax::{FileId, RootedPath, Source, VirtualPath, VirtualRoot},
        text::{Font, FontBook},
        utils::LazyHash,
    };
    use typst_bundle::{Bundle, BundleOptions, VirtualFs};
    use typst_kit::{
        datetime::Time,
        diagnostics::{DiagnosticFormat, DiagnosticWorld, emit, termcolor},
        downloader::SystemDownloader,
        files::{FileStore, FsRoot, SystemFiles},
        fonts::{self, FontStore},
        packages::{FsPackages, SystemPackages, UniversePackages},
    };
    use typst_pdf::Timestamp;

    use super::*;

    pub(crate) struct PersistentTypstRenderer {
        root: PathBuf,
        work_root: PathBuf,
        fonts: Option<Arc<FontStore>>,
        sessions: HashMap<String, Session>,
    }

    struct Session {
        world: DocsWorld,
        configuration: Configuration,
        compilations: usize,
    }

    #[derive(Clone, Debug, Eq, PartialEq)]
    struct Configuration {
        inputs: Vec<(String, String)>,
        timestamp: i64,
    }

    struct DocsWorld {
        main: FileId,
        library: LazyHash<Library>,
        fonts: Arc<FontStore>,
        files: FileStore<SystemFiles>,
        now: Time,
    }

    impl PersistentTypstRenderer {
        pub(crate) fn new(root: &Path, work_root: &Path) -> Result<Self> {
            let root = root
                .canonicalize()
                .wrap_err_with(|| format!("failed to resolve {}", root.display()))?;
            fs::create_dir_all(work_root)?;
            let work_root = work_root
                .canonicalize()
                .wrap_err_with(|| format!("failed to resolve {}", work_root.display()))?;
            ensure!(
                work_root.starts_with(&root),
                "persistent Typst work directory must be inside the workspace root"
            );
            Ok(Self {
                root,
                work_root,
                fonts: None,
                sessions: HashMap::new(),
            })
        }

        fn fonts(&mut self) -> Result<Arc<FontStore>> {
            if let Some(fonts) = &self.fonts {
                return Ok(fonts.clone());
            }
            let fonts = Arc::new(discover_fonts()?);
            self.fonts = Some(fonts.clone());
            Ok(fonts)
        }

        fn create_session(&mut self, key: &str, configuration: &Configuration) -> Result<Session> {
            ensure!(
                !key.is_empty()
                    && key
                        .bytes()
                        .all(|byte| byte.is_ascii_alphanumeric() || byte == b'-'),
                "invalid persistent Typst session key {key:?}"
            );
            let directory = self.work_root.join(key);
            fs::create_dir_all(&directory)?;
            let wrapper = directory.join("main.typ");
            if !wrapper.exists() {
                fs::write(&wrapper, "")?;
            }
            let vpath = VirtualPath::virtualize(&self.root, &wrapper)
                .map_err(|error| eyre!("invalid Typst wrapper path: {error}"))?;
            let main = RootedPath::new(VirtualRoot::Project, vpath).intern();
            let world = DocsWorld::new(self.root.clone(), main, self.fonts()?, configuration)?;
            Ok(Session {
                world,
                configuration: configuration.clone(),
                compilations: 0,
            })
        }
    }

    impl TypstRenderer for PersistentTypstRenderer {
        fn render(&mut self, job: TypstRenderJob) -> Result<()> {
            let configuration = Configuration {
                inputs: job.inputs,
                timestamp: job.timestamp,
            };
            if !self.sessions.contains_key(&job.key) {
                let session = self.create_session(&job.key, &configuration)?;
                self.sessions.insert(job.key.clone(), session);
            }
            let session = match self.sessions.entry(job.key.clone()) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(_) => unreachable!("session was inserted above"),
            };
            if session.configuration != configuration {
                session.world.configure(&configuration)?;
                session.configuration = configuration;
            }

            let wrapper = self.work_root.join(&job.key).join("main.typ");
            fs::write(&wrapper, job.wrapper)?;

            let pdf_timestamp = pdf_timestamp(job.timestamp)?;
            let Warned { output, warnings } = typst::compile::<Bundle>(&session.world);
            let output = output.and_then(|bundle| {
                let mut options = BundleOptions::default();
                options.pdf.timestamp = Some(pdf_timestamp);
                typst_bundle::export(&bundle, &options)
            });
            emit_diagnostics(&session.world, warnings.iter())?;
            if let Err(errors) = &output {
                emit_diagnostics(&session.world, errors.iter())?;
            }
            let dependencies = session.world.dependencies();
            session.world.reset();
            session.compilations += 1;
            if let Some(path) = &job.dependency_file {
                write_dependencies(path, &dependencies)?;
            }

            let virtual_fs = match output {
                Ok(output) => output,
                Err(_) => bail!("Typst failed for {}", job.label),
            };
            write_virtual_fs(&job.output, &virtual_fs)
        }

        fn finish_generation(&mut self) {
            // Comemo is process-global. Age it once per site generation, not once for
            // every product world in an all-product build.
            evict(10);
        }
    }

    impl DocsWorld {
        fn new(
            root: PathBuf,
            main: FileId,
            fonts: Arc<FontStore>,
            configuration: &Configuration,
        ) -> Result<Self> {
            let packages = SystemPackages::from_parts(
                env::var_os("TYPST_PACKAGE_PATH")
                    .map(FsPackages::new)
                    .or_else(FsPackages::system_data),
                env::var_os("TYPST_PACKAGE_CACHE_PATH")
                    .map(FsPackages::new)
                    .or_else(FsPackages::system_cache),
                UniversePackages::new(SystemDownloader::new(format!(
                    "alphal00p-docs-builder/{}",
                    env!("CARGO_PKG_VERSION")
                ))),
            );
            Ok(Self {
                files: FileStore::new(SystemFiles::new(FsRoot::new(root.clone()), packages)),
                main,
                library: make_library(&configuration.inputs),
                fonts,
                now: fixed_time(configuration.timestamp)?,
            })
        }

        fn configure(&mut self, configuration: &Configuration) -> Result<()> {
            self.library = make_library(&configuration.inputs);
            self.now = fixed_time(configuration.timestamp)?;
            Ok(())
        }

        fn dependencies(&mut self) -> BTreeSet<PathBuf> {
            let (loader, dependencies) = self.files.dependencies();
            dependencies
                .filter_map(|id| loader.resolve(id).ok())
                .collect()
        }

        fn reset(&mut self) {
            self.files.reset();
            self.now.reset();
        }
    }

    impl World for DocsWorld {
        fn library(&self) -> &LazyHash<Library> {
            &self.library
        }

        fn book(&self) -> &LazyHash<FontBook> {
            self.fonts.book()
        }

        fn main(&self) -> FileId {
            self.main
        }

        fn source(&self, id: FileId) -> FileResult<Source> {
            self.files.source(id)
        }

        fn file(&self, id: FileId) -> FileResult<Bytes> {
            self.files.file(id)
        }

        fn font(&self, index: usize) -> Option<Font> {
            self.fonts.font(index)
        }

        fn today(&self, offset: Option<Duration>) -> Option<Datetime> {
            self.now.today(offset)
        }
    }

    impl DiagnosticWorld for DocsWorld {
        fn name(&self, id: FileId) -> String {
            match id.root() {
                VirtualRoot::Project => id.vpath().get_without_slash().into(),
                VirtualRoot::Package(package) => {
                    format!("{package}{}", id.vpath().get_with_slash())
                }
            }
        }
    }

    fn make_library(inputs: &[(String, String)]) -> LazyHash<Library> {
        let inputs: Dict = inputs
            .iter()
            .map(|(key, value)| (key.as_str().into(), value.as_str().into_value()))
            .collect();
        let features = [Feature::Html, Feature::Bundle].into_iter().collect();
        LazyHash::new(
            Library::builder()
                .with_inputs(inputs)
                .with_features(features)
                .build(),
        )
    }

    fn fixed_time(timestamp: i64) -> Result<Time> {
        Time::fixed_timestamp(timestamp)
            .map_err(|error| eyre!("invalid Typst creation timestamp: {error}"))
    }

    fn pdf_timestamp(timestamp: i64) -> Result<Timestamp> {
        let timestamp = Utc
            .timestamp_opt(timestamp, 0)
            .single()
            .context("invalid PDF creation timestamp")?;
        let datetime = Datetime::from_ymd_hms(
            timestamp.year(),
            timestamp.month() as u8,
            timestamp.day() as u8,
            timestamp.hour() as u8,
            timestamp.minute() as u8,
            timestamp.second() as u8,
        )
        .context("invalid PDF creation timestamp")?;
        Ok(Timestamp::new_utc(datetime))
    }

    fn discover_fonts() -> Result<FontStore> {
        let mut store = FontStore::new();
        if !environment_flag("TYPST_IGNORE_SYSTEM_FONTS")? {
            store.extend(fonts::system());
        }
        if !environment_flag("TYPST_IGNORE_EMBEDDED_FONTS")? {
            store.extend(fonts::embedded());
        }
        if let Some(paths) = env::var_os("TYPST_FONT_PATHS") {
            for path in env::split_paths(&paths) {
                store.extend(fonts::scan(&path));
            }
        }
        Ok(store)
    }

    fn environment_flag(name: &str) -> Result<bool> {
        let Some(value) = env::var_os(name) else {
            return Ok(false);
        };
        let value = value
            .to_str()
            .wrap_err_with(|| format!("{name} is not valid UTF-8"))?;
        value
            .parse()
            .wrap_err_with(|| format!("{name} must be true or false"))
    }

    fn emit_diagnostics<'a>(
        world: &DocsWorld,
        diagnostics: impl IntoIterator<Item = &'a SourceDiagnostic>,
    ) -> Result<()> {
        let mut stream = termcolor::StandardStream::stderr(termcolor::ColorChoice::Auto);
        emit(&mut stream, world, diagnostics, DiagnosticFormat::Human)
            .wrap_err("failed to print Typst diagnostics")
    }

    fn write_dependencies(path: &Path, dependencies: &BTreeSet<PathBuf>) -> Result<()> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        let mut encoded = Vec::new();
        for dependency in dependencies {
            encoded.extend_from_slice(dependency.as_os_str().as_encoded_bytes());
            encoded.push(0);
        }
        fs::write(path, encoded)
            .wrap_err_with(|| format!("failed to write Typst dependencies to {}", path.display()))
    }

    fn write_virtual_fs(root: &Path, virtual_fs: &VirtualFs) -> Result<()> {
        fs::create_dir_all(root)?;
        for (path, data) in virtual_fs {
            let path = path
                .realize(root)
                .map_err(|error| eyre!("failed to realize Typst output path: {error}"))?;
            if let Some(parent) = path.parent() {
                fs::create_dir_all(parent)?;
            }
            fs::write(&path, data.as_slice())
                .wrap_err_with(|| format!("failed to write {}", path.display()))?;
        }
        Ok(())
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn persistent_world_reloads_an_import_without_recreation() {
            let root = tempfile::tempdir().unwrap();
            let work = root.path().join("target/watch/typst");
            let output = root.path().join("output");
            let dependency_file = root.path().join("dependencies/test.deps");
            fs::write(root.path().join("body.typ"), "First version").unwrap();
            let mut renderer = PersistentTypstRenderer::new(root.path(), &work).unwrap();

            let job = |output: &Path| TypstRenderJob {
                key: "test".to_owned(),
                label: "test".to_owned(),
                wrapper: "#document(\"index.html\")[#include \"/body.typ\"]".to_owned(),
                inputs: Vec::new(),
                timestamp: 1_700_000_000,
                output: output.to_path_buf(),
                dependency_file: Some(dependency_file.clone()),
            };
            renderer.render(job(&output)).unwrap();
            assert!(
                fs::read_to_string(output.join("index.html"))
                    .unwrap()
                    .contains("First version")
            );

            fs::write(root.path().join("body.typ"), "Second version").unwrap();
            renderer.render(job(&output)).unwrap();
            assert!(
                fs::read_to_string(output.join("index.html"))
                    .unwrap()
                    .contains("Second version")
            );
            let session = renderer.sessions.get("test").unwrap();
            assert_eq!(session.compilations, 2);
            assert_eq!(renderer.sessions.len(), 1);
            assert!(
                fs::read(&dependency_file)
                    .unwrap()
                    .split(|byte| *byte == 0)
                    .any(|path| path == root.path().join("body.typ").as_os_str().as_encoded_bytes())
            );
        }

        #[test]
        fn persistent_world_recovers_after_a_missing_import() {
            let root = tempfile::tempdir().unwrap();
            let work = root.path().join("target/watch/typst");
            let output = root.path().join("output");
            let dependencies = root.path().join("dependencies/test.deps");
            let mut renderer = PersistentTypstRenderer::new(root.path(), &work).unwrap();
            let job = || TypstRenderJob {
                key: "test".to_owned(),
                label: "test".to_owned(),
                wrapper: "#document(\"index.html\")[#include \"/missing.typ\"]".to_owned(),
                inputs: Vec::new(),
                timestamp: 1_700_000_000,
                output: output.clone(),
                dependency_file: Some(dependencies.clone()),
            };

            assert!(renderer.render(job()).is_err());
            assert!(
                fs::read(&dependencies)
                    .unwrap()
                    .windows("missing.typ".len())
                    .any(|window| window == b"missing.typ")
            );
            fs::write(root.path().join("missing.typ"), "Recovered").unwrap();
            renderer.render(job()).unwrap();
            assert!(
                fs::read_to_string(output.join("index.html"))
                    .unwrap()
                    .contains("Recovered")
            );
            assert_eq!(renderer.sessions["test"].compilations, 2);
        }

        #[test]
        fn persistent_bundle_matches_the_pinned_cli() {
            let Ok(version) = Command::new("typst").arg("--version").output() else {
                eprintln!("skipping CLI parity check because Typst is unavailable");
                return;
            };
            assert!(version.status.success());
            assert!(String::from_utf8_lossy(&version.stdout).contains("0.15.0"));

            let root = tempfile::tempdir().unwrap();
            let persistent_output = root.path().join("persistent-output");
            let cli_output = root.path().join("cli-output");
            let wrapper = "#document(\"index.html\")[Parity] \n\
                           #document(\"manual.pdf\")[Parity]";
            let job = |output: &Path| TypstRenderJob {
                key: "parity".to_owned(),
                label: "parity".to_owned(),
                wrapper: wrapper.to_owned(),
                inputs: vec![("test-input".to_owned(), "same".to_owned())],
                timestamp: 1_700_000_000,
                output: output.to_path_buf(),
                dependency_file: None,
            };

            let mut persistent =
                PersistentTypstRenderer::new(root.path(), &root.path().join("target/watch/typst"))
                    .unwrap();
            persistent.render(job(&persistent_output)).unwrap();
            let mut cli = CliTypstRenderer::new(root.path());
            cli.render(job(&cli_output)).unwrap();

            for file in ["index.html", "manual.pdf"] {
                assert_eq!(
                    fs::read(persistent_output.join(file)).unwrap(),
                    fs::read(cli_output.join(file)).unwrap(),
                    "persistent output differs from the Typst CLI for {file}"
                );
            }
        }
    }
}

#[cfg(feature = "persistent-typst")]
pub(crate) use persistent::PersistentTypstRenderer;

#[cfg(not(feature = "persistent-typst"))]
pub(crate) struct PersistentTypstRenderer;

#[cfg(not(feature = "persistent-typst"))]
impl PersistentTypstRenderer {
    pub(crate) fn new(_root: &Path, _work_root: &Path) -> Result<Self> {
        eyre::bail!(
            "the live watcher requires the persistent-typst feature; run `just docs-watch all` or pass `--features persistent-typst`"
        )
    }
}

#[cfg(not(feature = "persistent-typst"))]
impl TypstRenderer for PersistentTypstRenderer {
    fn render(&mut self, _job: TypstRenderJob) -> Result<()> {
        unreachable!("PersistentTypstRenderer::new always fails without its feature")
    }
}
