use std::{
    collections::BTreeSet,
    env,
    ffi::OsStr,
    fs, io,
    path::{Component, Path, PathBuf},
    process::Command,
    sync::mpsc::{self, Receiver, RecvTimeoutError},
    time::{Duration, Instant},
};

use eyre::{Context, ContextCompat, Result, bail, ensure};
use notify::{
    Event, EventKind, RecommendedWatcher, RecursiveMode, Watcher,
    event::{MetadataKind, ModifyKind},
};
use tempfile::{Builder as TempDirBuilder, TempDir};

use super::{BuildChannel, BuildRequest, typst_render::PersistentTypstRenderer};
use super::{SiteBuilder, absolute_from, copy_tree, server::LiveServer};

const QUIET_PERIOD: Duration = Duration::from_millis(100);
const MAX_BATCH: Duration = Duration::from_millis(500);
const LINNET_PYTHON_CRATES: &[&str] = &["linnet-py", "linnet"];
const SPENSO_PYTHON_CRATES: &[&str] = &["spynso3", "spenso", "spenso-macros", "spenso-hep-lib"];

#[derive(Clone, Debug)]
pub struct WatchRequest {
    pub product: String,
    pub output: PathBuf,
    pub port: u16,
    pub open: bool,
    pub serve: bool,
    pub include_rustdoc: bool,
}

struct SourceWatcher {
    watcher: RecommendedWatcher,
    events: Receiver<notify::Result<Event>>,
    root: PathBuf,
    excluded: Vec<PathBuf>,
    typst_dependencies: BTreeSet<PathBuf>,
    dependency_watch_roots: BTreeSet<PathBuf>,
}

impl SourceWatcher {
    fn new(root: &Path, excluded: Vec<PathBuf>) -> Result<Self> {
        let current_dir = env::current_dir().wrap_err("failed to resolve the current directory")?;
        let root = lexical_absolute(&current_dir, root)
            .canonicalize()
            .wrap_err_with(|| format!("failed to resolve documentation root {}", root.display()))?;
        let excluded = excluded
            .into_iter()
            .map(|path| canonicalize_existing_prefix(&root, &path))
            .collect();
        let (sender, events) = mpsc::channel();
        let mut watcher =
            notify::recommended_watcher(sender).wrap_err("failed to create source watcher")?;
        watcher
            .watch(&root, RecursiveMode::NonRecursive)
            .wrap_err_with(|| format!("failed to watch {}", root.display()))?;
        for directory in ["docs", "crates"] {
            let path = root.join(directory);
            if path.is_dir() {
                watcher
                    .watch(&path, RecursiveMode::Recursive)
                    .wrap_err_with(|| format!("failed to watch {}", path.display()))?;
            }
        }
        Ok(Self {
            watcher,
            events,
            root,
            excluded,
            typst_dependencies: BTreeSet::new(),
            dependency_watch_roots: BTreeSet::new(),
        })
    }

    fn replace_typst_dependencies(
        &mut self,
        dependencies: impl IntoIterator<Item = PathBuf>,
    ) -> Result<()> {
        let dependencies = dependencies
            .into_iter()
            .map(|path| lexical_absolute(&self.root, &path))
            .filter(|path| {
                !ignored_workspace_path(&self.root, path)
                    && !self
                        .excluded
                        .iter()
                        .any(|excluded| path.starts_with(excluded))
            })
            .collect::<BTreeSet<_>>();
        let desired_roots = dependencies
            .iter()
            .filter(|path| !self.baseline_covers(path))
            .filter_map(|path| closest_existing_parent(path))
            .collect::<BTreeSet<_>>();

        // Register replacements before removing stale roots so that an editor save cannot
        // fall into a gap while the dependency set changes.
        for root in desired_roots.difference(&self.dependency_watch_roots) {
            self.watcher
                .watch(root, RecursiveMode::NonRecursive)
                .wrap_err_with(|| {
                    format!("failed to watch Typst dependency parent {}", root.display())
                })?;
        }
        for root in self.dependency_watch_roots.difference(&desired_roots) {
            if let Err(error) = self.watcher.unwatch(root) {
                eprintln!(
                    "could not remove stale Typst dependency watch for {}: {error}",
                    root.display()
                );
            }
        }
        self.typst_dependencies = dependencies;
        self.dependency_watch_roots = desired_roots;
        Ok(())
    }

    fn refresh_typst_dependency_watches(&mut self) -> Result<()> {
        self.replace_typst_dependencies(self.typst_dependencies.clone())
    }

    fn extend_typst_dependencies(
        &mut self,
        dependencies: impl IntoIterator<Item = PathBuf>,
    ) -> Result<()> {
        let mut merged = self.typst_dependencies.clone();
        merged.extend(dependencies);
        self.replace_typst_dependencies(merged)
    }

    fn wait(&self) -> Vec<PathBuf> {
        let mut changed: Vec<PathBuf> = Vec::new();
        while changed.is_empty() {
            match self.events.recv() {
                Ok(Ok(event)) => self.collect(event, &mut changed),
                Ok(Err(error)) => eprintln!("documentation watch error: {error}"),
                Err(_) => return changed,
            }
        }

        let started = Instant::now();
        let mut last = Instant::now();
        loop {
            let quiet_remaining = QUIET_PERIOD.saturating_sub(last.elapsed());
            let batch_remaining = MAX_BATCH.saturating_sub(started.elapsed());
            let remaining = quiet_remaining.min(batch_remaining);
            if remaining.is_zero() {
                break;
            }
            match self.events.recv_timeout(remaining) {
                Ok(Ok(event)) => {
                    let previous = changed.len();
                    self.collect(event, &mut changed);
                    if changed.len() != previous {
                        last = Instant::now();
                    }
                }
                Ok(Err(error)) => eprintln!("documentation watch error: {error}"),
                Err(RecvTimeoutError::Timeout | RecvTimeoutError::Disconnected) => break,
            }
        }
        changed.sort();
        changed.dedup();
        changed
    }

    fn collect(&self, event: Event, changed: &mut Vec<PathBuf>) {
        if !rebuild_event(&event.kind) {
            return;
        }
        changed.extend(
            event
                .paths
                .into_iter()
                .filter(|path| self.relevant(path))
                .map(|path| path.strip_prefix(&self.root).unwrap_or(&path).to_path_buf()),
        );
    }

    fn relevant(&self, path: &Path) -> bool {
        let absolute = absolute_from(&self.root, path);
        if self
            .excluded
            .iter()
            .any(|excluded| absolute.starts_with(excluded))
        {
            return false;
        }
        if self
            .typst_dependencies
            .iter()
            .any(|dependency| absolute == *dependency || dependency.starts_with(&absolute))
        {
            return true;
        }
        let Ok(relative) = absolute.strip_prefix(&self.root) else {
            return false;
        };
        if relative.components().any(|component| {
            matches!(component, Component::Normal(name) if matches!(name.to_str(), Some(".git" | ".jj" | "target")))
        }) {
            return false;
        }
        let first = relative.components().next();
        match first {
            Some(Component::Normal(name)) if name == "docs" => !editor_artifact(relative),
            Some(Component::Normal(name)) if name == "crates" => {
                !editor_artifact(relative) && source_like(relative)
            }
            Some(Component::Normal(name)) => matches!(
                name.to_str(),
                Some("Cargo.toml" | "Cargo.lock" | "README.md" | "pyproject.toml")
            ),
            _ => false,
        }
    }

    fn baseline_covers(&self, path: &Path) -> bool {
        path.parent() == Some(self.root.as_path())
            || ["docs", "crates"]
                .into_iter()
                .any(|directory| path.starts_with(self.root.join(directory)))
    }
}

fn closest_existing_parent(path: &Path) -> Option<PathBuf> {
    let mut parent = path.parent();
    while let Some(candidate) = parent {
        if candidate.is_dir() {
            return Some(candidate.to_path_buf());
        }
        parent = candidate.parent();
    }
    None
}

fn lexical_absolute(root: &Path, path: &Path) -> PathBuf {
    let mut normalized = PathBuf::new();
    for component in absolute_from(root, path).components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                normalized.pop();
            }
            Component::Prefix(_) | Component::RootDir | Component::Normal(_) => {
                normalized.push(component.as_os_str());
            }
        }
    }
    normalized
}

fn canonicalize_existing_prefix(root: &Path, path: &Path) -> PathBuf {
    let absolute = lexical_absolute(root, path);
    let mut prefix = absolute.as_path();
    loop {
        if let Ok(canonical) = prefix.canonicalize() {
            let suffix = absolute.strip_prefix(prefix).unwrap_or(Path::new(""));
            return lexical_absolute(&canonical, suffix);
        }
        let Some(parent) = prefix.parent() else {
            return absolute;
        };
        prefix = parent;
    }
}

fn ignored_workspace_path(root: &Path, path: &Path) -> bool {
    path.strip_prefix(root).is_ok_and(|relative| {
        relative.components().any(|component| {
            matches!(component, Component::Normal(name) if matches!(name.to_str(), Some(".git" | ".jj" | "target")))
        })
    })
}

fn rebuild_event(kind: &EventKind) -> bool {
    match kind {
        EventKind::Create(_) | EventKind::Remove(_) => true,
        EventKind::Modify(ModifyKind::Metadata(MetadataKind::AccessTime)) => false,
        EventKind::Modify(ModifyKind::Metadata(_)) | EventKind::Access(_) => false,
        EventKind::Modify(_) | EventKind::Any => true,
        EventKind::Other => false,
    }
}

fn source_like(path: &Path) -> bool {
    matches!(
        path.extension().and_then(|extension| extension.to_str()),
        Some(
            "rs" | "toml"
                | "lock"
                | "md"
                | "typ"
                | "py"
                | "pyi"
                | "json"
                | "html"
                | "css"
                | "js"
                | "yaml"
                | "yml"
        )
    )
}

fn editor_artifact(path: &Path) -> bool {
    let Some(name) = path.file_name().and_then(|name| name.to_str()) else {
        return false;
    };
    name.starts_with(".#")
        || name.ends_with('~')
        || matches!(
            path.extension().and_then(|extension| extension.to_str()),
            Some("swp" | "swx" | "tmp")
        )
}

impl SiteBuilder {
    pub fn watch(&self, request: WatchRequest) -> Result<()> {
        self.watch_product_route(&request.product)?;
        ensure!(
            request.serve || !request.open,
            "--open requires the live server"
        );

        let target = self.root.join("target");
        fs::create_dir_all(&target)?;
        let session = TempDirBuilder::new()
            .prefix("alphal00p-docs-watch-")
            .tempdir_in(&target)?;
        let mut renderer = PersistentTypstRenderer::new(&self.root, &session.path().join("typst"))?;
        let output = absolute_from(&self.root, &request.output);
        let mut excluded = vec![session.path().to_path_buf()];
        if !request.serve {
            ensure_safe_watch_output(&self.root, &output)?;
            excluded.push(output.clone());
        }
        let mut watcher = SourceWatcher::new(&self.root, excluded)?;
        let server = request
            .serve
            .then(|| LiveServer::start(request.port))
            .transpose()?;
        if let Some(server) = &server {
            eprintln!(
                "serving documentation on http://{}{}",
                server.address(),
                self.watch_route(&request.product)
            );
        } else {
            eprintln!("writing documentation to {}", output.display());
        }

        let mut changed: Vec<PathBuf> = Vec::new();
        let mut open = request.open;
        loop {
            let dependency_output = session.path().join("dependencies");
            if let Some(server) = &server {
                server.set_status("Building documentation…");
            }
            eprintln!(
                "building {} documentation{}",
                request.product,
                changed_summary(&changed)
            );
            let build = self.build_generation(
                &request,
                session.path(),
                &dependency_output,
                &changed,
                &mut renderer,
            );
            let build_succeeded = build.is_ok();
            match read_typst_dependencies(&dependency_output) {
                Ok(Some(dependencies)) => {
                    if build_succeeded {
                        watcher.replace_typst_dependencies(dependencies)?;
                    } else {
                        // An all-product build can fail after emitting only a prefix of its
                        // dependency files. Merge that partial attempt into the last complete
                        // set so later products remain watched until a full build succeeds.
                        watcher.extend_typst_dependencies(dependencies)?;
                    }
                }
                Ok(None) => watcher.refresh_typst_dependency_watches()?,
                Err(error) => {
                    eprintln!(
                        "could not refresh Typst dependency watches; retaining the previous set: {error:#}"
                    );
                    watcher.refresh_typst_dependency_watches()?;
                }
            }
            match build {
                Ok((directory, root)) => {
                    if let Some(server) = &server {
                        server.publish(directory, root)?;
                        eprintln!("documentation rebuilt successfully");
                        if open {
                            let route = format!(
                                "http://{}{}",
                                server.address(),
                                self.watch_route(&request.product)
                            );
                            if let Err(error) = open_browser(&route) {
                                eprintln!("could not open {route}: {error:#}");
                            }
                            open = false;
                        }
                    } else {
                        replace_output_tree(&root, &output)?;
                        eprintln!("documentation rebuilt successfully");
                    }
                }
                Err(error) => {
                    eprintln!(
                        "documentation build failed; retaining the last good build\n{error:#}"
                    );
                    if let Some(server) = &server {
                        server.set_status(format!("Build failed: {error:#}"));
                    }
                }
            }
            changed = watcher.wait();
            if changed.is_empty() {
                bail!("documentation source watcher stopped");
            }
            if let Some(path) = changed.iter().find(|path| watcher_restart_required(path)) {
                bail!(
                    "{} changes the running documentation builder; restart the watcher to load it",
                    path.display()
                );
            }
        }
    }

    fn watch_product_route(&self, product: &str) -> Result<()> {
        if product == "all"
            || self
                .registry
                .product
                .iter()
                .any(|entry| entry.id == product)
        {
            Ok(())
        } else {
            bail!("unknown product {product}")
        }
    }

    fn watch_route(&self, product: &str) -> String {
        if product == "all" {
            "/".to_owned()
        } else {
            format!("/products/{product}/latest/")
        }
    }

    fn build_generation(
        &self,
        request: &WatchRequest,
        session: &Path,
        dependency_output: &Path,
        changed: &[PathBuf],
        renderer: &mut PersistentTypstRenderer,
    ) -> Result<(TempDir, PathBuf)> {
        if dependency_output.exists() {
            fs::remove_dir_all(dependency_output)?;
        }
        fs::create_dir_all(dependency_output)?;
        let directory = TempDirBuilder::new()
            .prefix("generation-")
            .tempdir_in(session)?;
        let root = directory.path().join("site");
        let api_root = directory.path().join("api");
        copy_tree(&self.root.join("docs/api"), &api_root)?;
        if changed.is_empty() {
            // Checked-in `.pyi` files are deliberate normal-watch inputs. Refresh only the
            // source-backed Rust references before the first publication so startup remains
            // Python-free; PyO3 inventories are isolated below and run after relevant changes.
            self.refresh_initial_api(&request.product, &api_root)?;
        } else if changed.iter().any(|path| rust_backed_change(path)) {
            self.refresh_api(&request.product, &api_root, changed)?;
        }
        let rustdoc_target = session.join("rustdoc");
        let mut builder = SiteBuilder::load(self.root.clone())?;
        builder.api_root = api_root;
        builder.build_with_renderer(
            BuildRequest {
                product: request.product.clone(),
                channel: BuildChannel::Latest,
                output: root.clone(),
                snapshot_tag: None,
                include_rustdoc: request.include_rustdoc,
                include_typst: true,
                rustdoc_target_root: Some(rustdoc_target),
                dependency_output: Some(dependency_output.to_path_buf()),
            },
            renderer,
        )?;
        Ok((directory, root))
    }

    fn refresh_initial_api(&self, product: &str, api_root: &Path) -> Result<()> {
        let (gammaloop, vakint) = initial_source_references(product);
        if gammaloop {
            self.refresh_gammaloop_reference(api_root)?;
        }
        if vakint {
            self.refresh_vakint_reference(api_root)?;
        }
        Ok(())
    }

    fn refresh_api(&self, product: &str, api_root: &Path, changed: &[PathBuf]) -> Result<()> {
        let global = changed
            .iter()
            .any(|path| matches!(path.to_str(), Some("Cargo.toml" | "Cargo.lock")));
        if matches!(product, "all" | "gammaloop")
            && (global || changed_in_crates(changed, &["gammaloop-api", "gammalooprs"]))
        {
            self.refresh_gammaloop_reference(api_root)?;
        }
        if matches!(product, "all" | "vakint")
            && (global || changed_in_crates(changed, &["vakint"]))
        {
            self.refresh_vakint_reference(api_root)?;
        }
        for (owner, feature, component, crates) in [
            (
                "gammaloop",
                "gammaloop",
                "gammaloop-python",
                &["gammaloop-api", "gammalooprs"][..],
            ),
            ("linnet", "linnet", "linnet-py", LINNET_PYTHON_CRATES),
            ("spenso", "spenso", "spynso3", SPENSO_PYTHON_CRATES),
            ("idenso", "idenso", "idenso-community", &["idenso"][..]),
            ("vakint", "vakint", "vakint-community", &["vakint"][..]),
        ] {
            if exporter_matches_change(product, owner, crates, global, changed) {
                self.run_exporter(
                    &[
                        "run",
                        "--locked",
                        "-p",
                        "alphal00p-docs-python-exporter",
                        "--features",
                        feature,
                        "--",
                        component,
                    ],
                    &api_root.join("python").join(format!("{component}.pyi")),
                )?;
            }
        }
        Ok(())
    }

    fn refresh_gammaloop_reference(&self, api_root: &Path) -> Result<()> {
        self.run_exporter(
            &[
                "run",
                "--locked",
                "-p",
                "alphal00p-docs-catalogs",
                "--features",
                "gammaloop-reference",
                "--bin",
                "alphal00p-docs-gammaloop-reference",
                "--",
                "--output",
            ],
            &api_root.join("generated/gammaloop-reference.json"),
        )
    }

    fn refresh_vakint_reference(&self, api_root: &Path) -> Result<()> {
        self.run_exporter(
            &[
                "run",
                "--locked",
                "-p",
                "alphal00p-docs-catalogs",
                "--features",
                "vakint-reference",
                "--bin",
                "alphal00p-docs-vakint-reference",
                "--",
                "--output",
            ],
            &api_root.join("generated/vakint-reference.json"),
        )
    }

    fn run_exporter(&self, arguments: &[&str], output: &Path) -> Result<()> {
        let status = Command::new(env::var_os("CARGO").unwrap_or_else(|| "cargo".into()))
            .current_dir(&self.root)
            .args(arguments)
            .arg(output)
            .status()
            .wrap_err("failed to launch isolated documentation exporter")?;
        ensure!(
            status.success(),
            "isolated documentation exporter failed for {}",
            output.display()
        );
        Ok(())
    }
}

fn read_typst_dependencies(directory: &Path) -> Result<Option<BTreeSet<PathBuf>>> {
    if !directory.is_dir() {
        return Ok(None);
    }
    let mut files = fs::read_dir(directory)?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<std::io::Result<Vec<_>>>()?;
    files.retain(|path| path.extension() == Some(OsStr::new("deps")));
    files.sort();
    if files.is_empty() {
        return Ok(None);
    }

    let mut dependencies = BTreeSet::new();
    for file in files {
        let encoded = fs::read(&file).wrap_err_with(|| {
            format!("failed to read Typst dependencies from {}", file.display())
        })?;
        dependencies.extend(
            encoded
                .split(|byte| *byte == 0)
                .filter(|path| !path.is_empty())
                .map(encoded_path),
        );
    }
    Ok(Some(dependencies))
}

fn encoded_path(bytes: &[u8]) -> PathBuf {
    // Typst's `zero` format writes the platform-native encoded bytes. The producer and
    // consumer run in the same process environment, satisfying this API's safety contract.
    PathBuf::from(unsafe { OsStr::from_encoded_bytes_unchecked(bytes) })
}

fn rust_backed_change(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| matches!(name, "Cargo.toml" | "Cargo.lock"))
        || path.extension().and_then(|extension| extension.to_str()) == Some("rs")
}

fn watcher_restart_required(path: &Path) -> bool {
    if matches!(path.to_str(), Some("Cargo.toml" | "Cargo.lock")) {
        return true;
    }
    let implementation_crate = [
        Path::new("crates/alphal00p-docs-builder"),
        Path::new("crates/alphal00p-docs-schema"),
    ]
    .into_iter()
    .any(|crate_root| path.starts_with(crate_root));
    implementation_crate
        && (matches!(
            path.file_name().and_then(|name| name.to_str()),
            Some("Cargo.toml" | "build.rs")
        ) || path.extension().and_then(|extension| extension.to_str()) == Some("rs"))
}

fn changed_in_crates(changed: &[PathBuf], crates: &[&str]) -> bool {
    crates.iter().any(|name| {
        let root = Path::new("crates").join(name);
        changed.iter().any(|path| path.starts_with(&root))
    })
}

fn initial_source_references(product: &str) -> (bool, bool) {
    (
        matches!(product, "all" | "gammaloop"),
        matches!(product, "all" | "vakint"),
    )
}

fn exporter_matches_change(
    selected_product: &str,
    owner: &str,
    crates: &[&str],
    global: bool,
    changed: &[PathBuf],
) -> bool {
    (selected_product == "all" || selected_product == owner)
        && (global || changed_in_crates(changed, crates))
}

fn changed_summary(paths: &[PathBuf]) -> String {
    if paths.is_empty() {
        return String::new();
    }
    let displayed = paths
        .iter()
        .take(4)
        .map(|path| path.display().to_string())
        .collect::<Vec<_>>()
        .join(", ");
    if paths.len() > 4 {
        format!(
            " after changes to {displayed}, and {} more",
            paths.len() - 4
        )
    } else {
        format!(" after changes to {displayed}")
    }
}

fn ensure_safe_watch_output(root: &Path, output: &Path) -> Result<()> {
    let target = root.join("target");
    ensure!(
        output.starts_with(&target) && output != target,
        "continuous documentation output must be a child of {}",
        target.display()
    );
    ensure!(
        !output
            .strip_prefix(&target)
            .expect("target descendant")
            .components()
            .any(|component| component == Component::ParentDir),
        "continuous documentation output cannot contain '..'"
    );
    let canonical_target = fs::canonicalize(&target)
        .wrap_err_with(|| format!("failed to resolve {}", target.display()))?;
    let existing_parent = closest_existing_parent(output)
        .context("continuous documentation output has no existing parent")?;
    let canonical_parent = fs::canonicalize(&existing_parent)
        .wrap_err_with(|| format!("failed to resolve {}", existing_parent.display()))?;
    ensure!(
        canonical_parent.starts_with(canonical_target),
        "continuous documentation output cannot escape the workspace target through a symlink"
    );
    Ok(())
}

fn replace_output_tree(source: &Path, output: &Path) -> Result<()> {
    let parent = output
        .parent()
        .context("documentation output has no parent")?;
    fs::create_dir_all(parent)?;
    let staging = TempDirBuilder::new()
        .prefix(".alphal00p-docs-output-")
        .tempdir_in(parent)?;
    let staged = staging.path().join("site");
    copy_tree(source, &staged)?;

    let backup = TempDirBuilder::new()
        .prefix(".alphal00p-docs-backup-")
        .tempdir_in(parent)?;
    let backup_path = backup.path().to_path_buf();
    backup.close()?;
    if path_exists(output) {
        fs::rename(output, &backup_path)
            .wrap_err_with(|| format!("failed to move {}", output.display()))?;
    }
    if let Err(error) = fs::rename(&staged, output) {
        if path_exists(&backup_path) {
            let _ = fs::rename(&backup_path, output);
        }
        return Err(error).wrap_err_with(|| format!("failed to publish {}", output.display()));
    }
    if path_exists(&backup_path) {
        remove_path(&backup_path)?;
    }
    Ok(())
}

fn path_exists(path: &Path) -> bool {
    fs::symlink_metadata(path).is_ok()
}

fn remove_path(path: &Path) -> Result<()> {
    let metadata = match fs::symlink_metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(()),
        Err(error) => return Err(error.into()),
    };
    if metadata.is_dir() && !metadata.file_type().is_symlink() {
        fs::remove_dir_all(path)?;
    } else {
        fs::remove_file(path)?;
    }
    Ok(())
}

fn open_browser(target: &str) -> Result<()> {
    #[cfg(target_os = "macos")]
    let mut command = {
        let mut command = Command::new("open");
        command.arg(target);
        command
    };
    #[cfg(target_os = "windows")]
    let mut command = {
        let mut command = Command::new("cmd");
        command.args(["/C", "start", "", target]);
        command
    };
    #[cfg(all(unix, not(target_os = "macos")))]
    let mut command = {
        let mut command = Command::new("xdg-open");
        command.arg(target);
        command
    };
    command
        .spawn()
        .wrap_err_with(|| format!("failed to open {target}"))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use notify::event::RenameMode;

    fn workspace() -> TempDir {
        let directory = tempfile::tempdir().unwrap();
        fs::create_dir(directory.path().join("docs")).unwrap();
        fs::create_dir(directory.path().join("crates")).unwrap();
        directory
    }

    #[test]
    fn event_filter_ignores_access_and_metadata_only_changes() {
        assert!(!rebuild_event(&EventKind::Access(
            notify::event::AccessKind::Read
        )));
        assert!(!rebuild_event(&EventKind::Modify(ModifyKind::Metadata(
            MetadataKind::Permissions
        ))));
        assert!(rebuild_event(&EventKind::Modify(ModifyKind::Data(
            notify::event::DataChange::Content
        ))));
        assert!(rebuild_event(&EventKind::Create(
            notify::event::CreateKind::File
        )));
        assert!(rebuild_event(&EventKind::Remove(
            notify::event::RemoveKind::File
        )));
        assert!(!rebuild_event(&EventKind::Other));
    }

    #[test]
    fn runtime_exports_refresh_only_for_rust_backed_changes() {
        assert!(rust_backed_change(Path::new("crates/spynso3/src/lib.rs")));
        assert!(rust_backed_change(Path::new("Cargo.toml")));
        assert!(!rust_backed_change(Path::new(
            "docs/products/spenso/content/api.typ"
        )));
        assert!(!rust_backed_change(Path::new(
            "docs/api/python/spynso3.pyi"
        )));
        let changed = vec![PathBuf::from("crates/spenso/src/lib.rs")];
        assert!(changed_in_crates(&changed, &["spenso"]));
        assert!(!changed_in_crates(&changed, &["spynso3"]));
        assert!(exporter_matches_change(
            "spenso",
            "spenso",
            SPENSO_PYTHON_CRATES,
            false,
            &changed,
        ));
        assert!(!exporter_matches_change(
            "idenso",
            "spenso",
            SPENSO_PYTHON_CRATES,
            false,
            &changed,
        ));
        let linnet_change = vec![PathBuf::from("crates/linnet/src/half_edge.rs")];
        assert!(exporter_matches_change(
            "linnet",
            "linnet",
            LINNET_PYTHON_CRATES,
            false,
            &linnet_change,
        ));
        assert!(!exporter_matches_change(
            "spenso",
            "linnet",
            LINNET_PYTHON_CRATES,
            false,
            &linnet_change,
        ));
        let gamma_change = vec![PathBuf::from("crates/gammalooprs/src/processes/mod.rs")];
        assert!(exporter_matches_change(
            "gammaloop",
            "gammaloop",
            &["gammaloop-api", "gammalooprs"],
            false,
            &gamma_change,
        ));
        assert!(!exporter_matches_change(
            "linnet",
            "gammaloop",
            &["gammaloop-api", "gammalooprs"],
            false,
            &gamma_change,
        ));
    }

    #[test]
    fn watcher_requests_a_restart_only_for_its_own_compiled_inputs() {
        for path in [
            "Cargo.toml",
            "Cargo.lock",
            "crates/alphal00p-docs-builder/Cargo.toml",
            "crates/alphal00p-docs-builder/src/typst_render.rs",
            "crates/alphal00p-docs-schema/src/lib.rs",
        ] {
            assert!(watcher_restart_required(Path::new(path)), "{path}");
        }
        for path in [
            "docs/products/gammaloop/content/tutorial.typ",
            "crates/gammalooprs/src/lib.rs",
            "crates/spenso/Cargo.toml",
        ] {
            assert!(!watcher_restart_required(Path::new(path)), "{path}");
        }
    }

    #[test]
    fn initial_refresh_is_limited_to_python_free_source_references() {
        assert_eq!(initial_source_references("gammaloop"), (true, false));
        assert_eq!(initial_source_references("vakint"), (false, true));
        assert_eq!(initial_source_references("all"), (true, true));
        for product in ["linnet", "spenso", "idenso"] {
            assert_eq!(initial_source_references(product), (false, false));
        }
    }

    #[test]
    fn editor_files_are_not_documentation_inputs() {
        assert!(editor_artifact(Path::new("docs/products/.#overview.typ")));
        assert!(editor_artifact(Path::new("docs/products/overview.typ~")));
        assert!(editor_artifact(Path::new("crates/spenso/src/lib.rs.swp")));
        assert!(!editor_artifact(Path::new("crates/spenso/src/lib.rs")));
    }

    #[test]
    fn nested_repository_and_build_directories_are_irrelevant() {
        let root = workspace();
        let watcher = SourceWatcher::new(root.path(), Vec::new()).unwrap();
        for path in [
            "crates/spenso/target/debug/generated.rs",
            "crates/vakint/.git/index",
            "docs/products/idenso/.jj/repo/op_heads",
        ] {
            assert!(!watcher.relevant(&root.path().join(path)), "{path}");
        }
        assert!(watcher.relevant(&root.path().join("crates/spenso/src/lib.rs")));
    }

    #[test]
    fn zero_dependency_files_are_merged_deterministically() {
        let directory = tempfile::tempdir().unwrap();
        fs::write(
            directory.path().join("spenso.deps"),
            b"docs/a.typ\0docs/b.typ\0",
        )
        .unwrap();
        fs::write(
            directory.path().join("idenso.deps"),
            b"docs/b.typ\0docs/c.typ\0",
        )
        .unwrap();
        fs::write(directory.path().join("ignored.json"), b"[]").unwrap();

        let dependencies = read_typst_dependencies(directory.path())
            .unwrap()
            .unwrap()
            .into_iter()
            .collect::<Vec<_>>();
        assert_eq!(
            dependencies,
            ["docs/a.typ", "docs/b.typ", "docs/c.typ"]
                .map(PathBuf::from)
                .to_vec()
        );
    }

    #[test]
    fn dependency_paths_are_normalized_before_matching_events() {
        let root = workspace();
        let dependency = root.path().join("docs/theme.typ");
        let mut watcher = SourceWatcher::new(root.path(), Vec::new()).unwrap();
        watcher
            .replace_typst_dependencies([PathBuf::from("docs/content/../theme.typ")])
            .unwrap();
        assert_eq!(watcher.typst_dependencies, BTreeSet::from([dependency]));
    }

    #[test]
    fn absent_dependency_files_preserve_the_previous_watch_set() {
        let directory = tempfile::tempdir().unwrap();
        assert!(read_typst_dependencies(directory.path()).unwrap().is_none());
    }

    #[test]
    fn failed_build_dependencies_extend_instead_of_replacing_the_last_complete_set() {
        let root = workspace();
        let mut watcher = SourceWatcher::new(root.path(), Vec::new()).unwrap();
        let first = root.path().join("docs/first.typ");
        let second = root.path().join("docs/second.typ");
        let changed = root.path().join("docs/changed.typ");

        watcher
            .replace_typst_dependencies([first.clone(), second.clone()])
            .unwrap();
        watcher
            .extend_typst_dependencies([changed.clone()])
            .unwrap();
        assert_eq!(
            watcher.typst_dependencies,
            BTreeSet::from([first, second, changed.clone()])
        );

        watcher
            .replace_typst_dependencies([changed.clone()])
            .unwrap();
        assert_eq!(watcher.typst_dependencies, BTreeSet::from([changed]));
    }

    #[test]
    fn missing_dependencies_watch_the_nearest_existing_parent() {
        let root = workspace();
        let external = tempfile::tempdir().unwrap();
        let inputs = external.path().join("inputs");
        fs::create_dir(&inputs).unwrap();
        let nested = inputs.join("theme");
        let dependency = nested.join("colors.typ");
        let mut watcher = SourceWatcher::new(root.path(), Vec::new()).unwrap();

        watcher
            .replace_typst_dependencies([dependency.clone()])
            .unwrap();
        assert_eq!(watcher.dependency_watch_roots, BTreeSet::from([inputs]));

        fs::create_dir(&nested).unwrap();
        let mut changed = Vec::new();
        watcher.collect(
            Event::new(EventKind::Create(notify::event::CreateKind::Folder))
                .add_path(nested.clone()),
            &mut changed,
        );
        assert_eq!(changed, vec![nested.clone()]);
        watcher.refresh_typst_dependency_watches().unwrap();
        assert_eq!(watcher.typst_dependencies, BTreeSet::from([dependency]));
        assert_eq!(watcher.dependency_watch_roots, BTreeSet::from([nested]));
    }

    #[test]
    fn replacing_dependencies_removes_stale_watch_roots() {
        let root = workspace();
        let external = tempfile::tempdir().unwrap();
        let first = external.path().join("first");
        let second = external.path().join("second");
        fs::create_dir(&first).unwrap();
        fs::create_dir(&second).unwrap();
        let mut watcher = SourceWatcher::new(root.path(), Vec::new()).unwrap();

        watcher
            .replace_typst_dependencies([first.join("manual.typ")])
            .unwrap();
        assert_eq!(watcher.dependency_watch_roots, BTreeSet::from([first]));
        watcher
            .replace_typst_dependencies([second.join("manual.typ")])
            .unwrap();
        assert_eq!(watcher.dependency_watch_roots, BTreeSet::from([second]));
    }

    #[test]
    fn deleted_dependency_parent_falls_back_without_stopping_the_watcher() {
        let root = workspace();
        let external = tempfile::tempdir().unwrap();
        let nested = external.path().join("theme");
        fs::create_dir(&nested).unwrap();
        let dependency = nested.join("colors.typ");
        let mut watcher = SourceWatcher::new(root.path(), Vec::new()).unwrap();
        watcher
            .replace_typst_dependencies([dependency.clone()])
            .unwrap();

        fs::remove_dir(&nested).unwrap();
        watcher.refresh_typst_dependency_watches().unwrap();
        assert_eq!(
            watcher.dependency_watch_roots,
            BTreeSet::from([external.path().to_path_buf()])
        );
        assert_eq!(watcher.typst_dependencies, BTreeSet::from([dependency]));
    }

    #[test]
    fn atomic_replacement_reports_only_the_dependency_target() {
        let root = workspace();
        let external = tempfile::tempdir().unwrap();
        let dependency = external.path().join("theme.typ");
        fs::write(&dependency, "#let color = blue").unwrap();
        let temporary = external.path().join(".theme.typ.tmp");
        let mut watcher = SourceWatcher::new(root.path(), Vec::new()).unwrap();
        watcher
            .replace_typst_dependencies([dependency.clone()])
            .unwrap();

        let event = Event::new(EventKind::Modify(ModifyKind::Name(RenameMode::Both)))
            .add_path(temporary)
            .add_path(dependency.clone());
        let mut changed = Vec::new();
        watcher.collect(event, &mut changed);
        assert_eq!(changed, vec![dependency]);
    }

    #[test]
    fn generated_outputs_never_become_dynamic_dependencies() {
        let root = workspace();
        let output = root.path().join("target/docs-watch");
        fs::create_dir_all(&output).unwrap();
        let mut watcher = SourceWatcher::new(root.path(), vec![output.clone()]).unwrap();
        watcher
            .replace_typst_dependencies([output.join("index.html")])
            .unwrap();
        assert!(watcher.typst_dependencies.is_empty());
        assert!(watcher.dependency_watch_roots.is_empty());

        let wrapper = root.path().join("target/alphal00p-typst-build/site.typ");
        fs::create_dir_all(wrapper.parent().unwrap()).unwrap();
        fs::write(&wrapper, "#import \"/docs/theme.typ\"").unwrap();
        let mut watcher = SourceWatcher::new(root.path(), Vec::new()).unwrap();
        watcher.replace_typst_dependencies([wrapper]).unwrap();
        assert!(watcher.typst_dependencies.is_empty());
        assert!(watcher.dependency_watch_roots.is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn symlinked_workspace_uses_canonical_watch_and_exclusion_paths() {
        use std::os::unix::fs::symlink;

        let container = tempfile::tempdir().unwrap();
        let root = container.path().join("workspace");
        fs::create_dir(&root).unwrap();
        fs::create_dir(root.join("docs")).unwrap();
        fs::create_dir(root.join("crates")).unwrap();
        fs::create_dir(root.join("target")).unwrap();
        fs::create_dir(root.join("target/session")).unwrap();
        let alias = container.path().join("workspace-link");
        symlink(&root, &alias).unwrap();

        let existing_exclusion = alias.join("target/session");
        let missing_exclusion = alias.join("generated/site");
        let mut watcher = SourceWatcher::new(
            &alias,
            vec![existing_exclusion.clone(), missing_exclusion.clone()],
        )
        .unwrap();

        assert_eq!(watcher.root, root.canonicalize().unwrap());
        assert_eq!(
            watcher.excluded,
            vec![root.join("target/session"), root.join("generated/site")]
        );

        let wrapper = root.join("target/session/main.typ");
        watcher.replace_typst_dependencies([wrapper]).unwrap();
        assert!(watcher.typst_dependencies.is_empty());
    }

    #[test]
    fn continuous_output_is_confined_below_the_workspace_target() {
        let workspace = tempfile::tempdir().unwrap();
        let root = workspace.path();
        fs::create_dir(root.join("target")).unwrap();
        assert!(ensure_safe_watch_output(root, &root.join("target/docs-watch")).is_ok());
        assert!(ensure_safe_watch_output(root, &root.join("target")).is_err());
        assert!(ensure_safe_watch_output(root, &root.join("docs")).is_err());
        assert!(ensure_safe_watch_output(root, Path::new("/tmp/docs-watch")).is_err());
    }

    #[test]
    fn continuous_output_replaces_an_existing_file() {
        let workspace = tempfile::tempdir().unwrap();
        let output = workspace.path().join("target/docs-watch");
        fs::create_dir_all(output.parent().unwrap()).unwrap();
        fs::write(&output, "stale file").unwrap();
        let source = tempfile::tempdir().unwrap();
        fs::write(source.path().join("index.html"), "new site").unwrap();

        replace_output_tree(source.path(), &output).unwrap();

        assert!(output.is_dir());
        assert_eq!(
            fs::read_to_string(output.join("index.html")).unwrap(),
            "new site"
        );
    }
}
