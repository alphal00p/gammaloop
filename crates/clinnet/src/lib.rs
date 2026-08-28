use std::borrow::Cow;
use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, Mutex, OnceLock};

#[cfg(test)]
use std::sync::atomic::{AtomicU64, Ordering};

use eyre::{Context, Result, bail, eyre};
use rust_embed::RustEmbed;

mod typst;

pub use typst::{TypstAngleUnit, TypstConfig, TypstLengthUnit, TypstMathScript, TypstValue};

const DEFAULT_TEMPLATE_SUBDIR: &str = ".clinnet/templates";
const LINNEST_PACKAGE_DIR: &str = "crates/linnest/typst";
const KURVST_PACKAGE_DIR: &str = "crates/kurvst/typst";
const MINIMUM_TYPST_VERSION: (u32, u32, u32) = (0, 15, 0);
#[cfg(test)]
static SCRATCH_ID: AtomicU64 = AtomicU64::new(0);
static ASSET_STAGING_LOCK: Mutex<()> = Mutex::new(());

#[derive(RustEmbed)]
#[folder = "templates"]
struct EmbeddedTemplates;

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/../linnest/typst"]
#[include = "src/*.typ"]
#[include = "src/**/*.typ"]
#[include = "typst.toml"]
#[include = "linnest.wasm"]
struct EmbeddedLinnestPackage;

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/../kurvst/typst"]
#[include = "src/*.typ"]
#[include = "src/**/*.typ"]
#[include = "typst.toml"]
#[include = "kurvst.wasm"]
struct EmbeddedKurvstPackage;

#[derive(Copy, Clone)]
enum TemplateKind {
    Figure,
    Grid,
}

impl TemplateKind {
    /// Return the embedded filename associated with this template kind.
    fn file_name(self) -> &'static str {
        match self {
            TemplateKind::Figure => "figure.typ",
            TemplateKind::Grid => "grid.typ",
        }
    }

    fn embedded_bytes(self) -> Result<Cow<'static, [u8]>> {
        EmbeddedTemplates::get(self.file_name())
            .map(|file| file.data)
            .ok_or_else(|| eyre!("embedded template {} is missing", self.file_name()))
    }
}

struct RenderContext<'a> {
    input: &'a Path,
    template: &'a Path,
    root: &'a Path,
    output: &'a Path,
}

/// A statically imported Typst module referenced by a generated render config.
///
/// Clinnet validates `alias` before serializing module symbols into the
/// entrypoint. The module path also participates in selection of the Typst
/// project root.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TypstModule {
    alias: String,
    source: TypstModuleSource,
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum TypstModuleSource {
    File(PathBuf),
    Package(String),
}

impl TypstModule {
    /// Import a filesystem module.
    pub fn new(alias: impl Into<String>, path: impl Into<PathBuf>) -> Result<Self> {
        let alias = alias.into();
        Self::validate_alias(&alias)?;
        Ok(Self {
            alias,
            source: TypstModuleSource::File(path.into()),
        })
    }

    /// Import a Typst package such as `@preview/cetz:0.5.1`.
    pub fn package(alias: impl Into<String>, package: impl Into<String>) -> Result<Self> {
        let alias = alias.into();
        let package = package.into();
        Self::validate_alias(&alias)?;
        Self::validate_package(&package)?;
        Ok(Self {
            alias,
            source: TypstModuleSource::Package(package),
        })
    }

    fn validate_alias(alias: &str) -> Result<()> {
        let mut chars = alias.chars();
        let valid = chars
            .next()
            .is_some_and(|first| first == '_' || first.is_ascii_alphabetic())
            && chars.all(|character| character == '_' || character.is_ascii_alphanumeric());
        if !valid || matches!(alias, "_clinnet_template" | "_clinnet_config") {
            bail!("invalid or reserved Typst module alias {alias:?}");
        }
        Ok(())
    }

    fn validate_package(package: &str) -> Result<()> {
        let Some((namespace, name_version)) = package
            .strip_prefix('@')
            .and_then(|value| value.split_once('/'))
        else {
            bail!("invalid Typst package source {package:?}");
        };
        let Some((name, version)) = name_version.split_once(':') else {
            bail!("invalid Typst package source {package:?}");
        };
        let valid_component = |value: &str| {
            !value.is_empty()
                && value.chars().all(|character| {
                    character.is_ascii_alphanumeric() || matches!(character, '_' | '-' | '.' | '+')
                })
        };
        if !valid_component(namespace) || !valid_component(name) || !valid_component(version) {
            bail!("invalid Typst package source {package:?}");
        }
        Ok(())
    }
}

/// One V1 Typst render request.
///
/// Clinnet serializes the closed native `config`, owns filesystem paths,
/// statically imports referenced modules, injects `data-path`, and invokes the
/// selected template's mandatory `render(config)` export.
#[derive(Clone, Debug)]
pub struct TypstRenderRequest {
    dot: String,
    config: TypstConfig,
    template: Option<PathBuf>,
    modules: Vec<TypstModule>,
    dependencies: Vec<PathBuf>,
}

impl TypstRenderRequest {
    pub fn new(dot: impl Into<String>, config: TypstConfig) -> Self {
        Self {
            dot: dot.into(),
            config,
            template: None,
            modules: Vec::new(),
            dependencies: Vec::new(),
        }
    }

    /// Select a V1 template exporting `render(config)`.
    pub fn template(mut self, template: impl Into<PathBuf>) -> Self {
        self.template = Some(template.into());
        self
    }

    /// Statically import one module under the alias used by `config`.
    pub fn module(mut self, module: TypstModule) -> Self {
        self.modules.push(module);
        self
    }

    pub fn modules(mut self, modules: impl IntoIterator<Item = TypstModule>) -> Self {
        self.modules.extend(modules);
        self
    }

    /// Include a transitive source or asset when choosing the Typst root.
    pub fn dependency(mut self, dependency: impl Into<PathBuf>) -> Self {
        self.dependencies.push(dependency.into());
        self
    }

    pub fn dependencies(
        mut self,
        dependencies: impl IntoIterator<Item = impl Into<PathBuf>>,
    ) -> Self {
        self.dependencies
            .extend(dependencies.into_iter().map(Into::into));
        self
    }
}

/// A staged Typst entrypoint and DOT input.
///
/// The entrypoint and project root can be passed to an in-process compiler such
/// as typst-py. The caller must also keep the renderer build directory and any
/// external template, module, or dependency paths alive until compilation
/// finishes. Dropping this value removes its ephemeral entrypoint and DOT input.
#[derive(Debug)]
pub struct PreparedTypstRender {
    _scratch_dir: tempfile::TempDir,
    entrypoint: PathBuf,
    root: PathBuf,
    dot_path: PathBuf,
    template: PathBuf,
}

impl PreparedTypstRender {
    /// Return the generated Typst entrypoint.
    pub fn entrypoint(&self) -> &Path {
        &self.entrypoint
    }

    /// Return the project root containing every filesystem compilation input.
    pub fn root(&self) -> &Path {
        &self.root
    }
}

/// Compile Linnet DOT graphs through an external Typst executable.
///
/// The renderer stages the bundled templates and Typst packages beneath its
/// build directory. Every render compiles an ephemeral V1 entrypoint, keeping
/// typed configuration out of process arguments.
#[derive(Clone, Debug)]
pub struct TypstRenderer {
    build_dir: PathBuf,
    typst_executable: PathBuf,
    checked_version: Arc<OnceLock<String>>,
}

impl TypstRenderer {
    pub fn new(build_dir: impl Into<PathBuf>) -> Self {
        Self {
            build_dir: build_dir.into(),
            typst_executable: PathBuf::from("typst"),
            checked_version: Arc::new(OnceLock::new()),
        }
    }

    /// Select a Typst executable other than the `typst` found in `PATH`.
    pub fn typst_executable(mut self, executable: impl Into<PathBuf>) -> Self {
        self.typst_executable = executable.into();
        self.checked_version = Arc::new(OnceLock::new());
        self
    }

    /// Check that the configured executable is Typst 0.15.0 or newer.
    pub fn check_version(&self) -> Result<String> {
        if let Some(version) = self.checked_version.get() {
            return Ok(version.clone());
        }
        let output = Command::new(&self.typst_executable)
            .arg("--version")
            .output()
            .with_context(|| {
                format!(
                    "failed to run '{} --version'. Is Typst installed and in PATH?",
                    self.typst_executable.display()
                )
            })?;

        if !output.status.success() {
            bail!("{} --version failed", self.typst_executable.display());
        }

        let version_output = String::from_utf8_lossy(&output.stdout);
        let version_line = version_output.lines().next().unwrap_or("unknown");
        Self::require_typst_version(version_line)?;
        let version_line = version_line.to_owned();
        let _ = self.checked_version.set(version_line.clone());
        Ok(version_line)
    }

    /// Ensure every embedded default dependency exists in Clinnet's private build namespace.
    pub fn stage_default_assets(&self) -> Result<()> {
        let _staging = ASSET_STAGING_LOCK
            .lock()
            .map_err(|_| eyre!("embedded asset staging lock was poisoned"))?;
        let template_dir = self.build_dir.join(DEFAULT_TEMPLATE_SUBDIR);
        // Bundled paths always represent this renderer version. Custom V1
        // templates are selected explicitly on the render request.
        Self::write_embedded_assets::<EmbeddedTemplates>(
            &template_dir,
            "template dependency",
            true,
        )?;
        Self::write_embedded_assets::<EmbeddedLinnestPackage>(
            &template_dir.join(LINNEST_PACKAGE_DIR),
            "linnest typst package asset",
            true,
        )?;
        Self::write_embedded_assets::<EmbeddedKurvstPackage>(
            &template_dir.join(KURVST_PACKAGE_DIR),
            "kurvst typst package asset",
            true,
        )?;
        Ok(())
    }

    /// Resolve a figure template, creating the bundled default inside the build tree if needed.
    pub fn resolve_figure_template(&self, requested: impl AsRef<Path>) -> Result<PathBuf> {
        self.resolve_template(requested.as_ref(), TemplateKind::Figure)
    }

    /// Resolve a grid template, creating the bundled default inside the build tree if needed.
    pub fn resolve_grid_template(&self, requested: impl AsRef<Path>) -> Result<PathBuf> {
        self.resolve_template(requested.as_ref(), TemplateKind::Grid)
    }

    /// Render a typed V1 request to a PDF, SVG, or PNG selected by `output`.
    pub fn render(&self, request: &TypstRenderRequest, output: impl AsRef<Path>) -> Result<()> {
        let output = output.as_ref();
        Self::validate_output_path(output)?;
        let prepared = self.prepare(request)?;
        self.check_version()?;
        Self::ensure_parent_dir(output)?;
        let mut command = Command::new(&self.typst_executable);
        command
            .arg("compile")
            .arg(prepared.entrypoint())
            .arg(output)
            .arg("--root")
            .arg(prepared.root());
        self.run_typst(
            &mut command,
            &format!("rendering {}", output.display()),
            Some(RenderContext {
                input: &prepared.dot_path,
                template: &prepared.template,
                root: prepared.root(),
                output,
            }),
        )
    }

    /// Stage one render request without selecting a compiler backend.
    pub fn prepare(&self, request: &TypstRenderRequest) -> Result<PreparedTypstRender> {
        Self::validate_modules(&request.modules)?;
        for alias in request.config.module_aliases() {
            if !request.modules.iter().any(|module| module.alias == alias) {
                bail!("render config references missing Typst module alias {alias:?}");
            }
        }
        fs::create_dir_all(&self.build_dir).with_context(|| {
            format!(
                "failed to create renderer build directory {}",
                self.build_dir.display()
            )
        })?;
        let build_dir = Self::canonicalize_existing(&self.build_dir)?;
        let template = if let Some(template) = &request.template {
            Self::canonicalize_existing(template)
                .with_context(|| format!("failed to read template {}", template.display()))?
        } else {
            self.stage_default_assets()?;
            self.resolve_figure_template(
                self.build_dir
                    .join(DEFAULT_TEMPLATE_SUBDIR)
                    .join(TemplateKind::Figure.file_name()),
            )?
        };
        let scratch_dir = tempfile::Builder::new()
            .prefix(".clinnet-render-")
            .tempdir_in(&build_dir)
            .with_context(|| {
                format!(
                    "failed to create render scratch directory in {}",
                    build_dir.display()
                )
            })?;
        let dot_path = scratch_dir.path().join("diagram.dot");
        fs::write(&dot_path, &request.dot)
            .with_context(|| format!("failed to stage DOT input {}", dot_path.display()))?;

        let module_files = request
            .modules
            .iter()
            .map(|module| {
                if let TypstModuleSource::File(path) = &module.source {
                    Self::canonicalize_existing(path)
                        .with_context(|| format!("failed to read Typst module {}", path.display()))
                        .map(Some)
                } else {
                    Ok(None)
                }
            })
            .collect::<Result<Vec<_>>>()?;
        let dependencies = request
            .dependencies
            .iter()
            .map(|dependency| {
                Self::canonicalize_existing(dependency).with_context(|| {
                    format!("failed to read Typst dependency {}", dependency.display())
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let mut root_inputs = Vec::with_capacity(module_files.len() + dependencies.len() + 2);
        root_inputs.push(template.as_path());
        root_inputs.push(dot_path.as_path());
        root_inputs.extend(module_files.iter().flatten().map(PathBuf::as_path));
        root_inputs.extend(dependencies.iter().map(PathBuf::as_path));
        let root = Self::compilation_root(&root_inputs)?;

        let entrypoint = scratch_dir.path().join("render.typ");
        let source = Self::entrypoint_source(request, &template, &module_files, &dot_path, &root)?;
        fs::write(&entrypoint, source).with_context(|| {
            format!("failed to stage Typst entrypoint {}", entrypoint.display())
        })?;

        Ok(PreparedTypstRender {
            _scratch_dir: scratch_dir,
            entrypoint,
            root,
            dot_path,
            template,
        })
    }

    /// Render a typed V1 request and return the resulting SVG document.
    pub fn to_svg(&self, request: &TypstRenderRequest) -> Result<String> {
        fs::create_dir_all(&self.build_dir).with_context(|| {
            format!(
                "failed to create renderer build directory {}",
                self.build_dir.display()
            )
        })?;
        let output_dir = tempfile::Builder::new()
            .prefix(".clinnet-output-")
            .tempdir_in(&self.build_dir)
            .with_context(|| {
                format!(
                    "failed to create render output directory in {}",
                    self.build_dir.display()
                )
            })?;
        let output = output_dir.path().join("diagram.svg");
        self.render(request, &output)?;
        fs::read_to_string(&output)
            .with_context(|| format!("failed to read rendered SVG {}", output.display()))
    }

    /// Compile a Typst template whose dependencies determine the project root.
    pub fn compile_template(
        &self,
        template: impl AsRef<Path>,
        output: impl AsRef<Path>,
        dependencies: &[&Path],
    ) -> Result<()> {
        let template = template.as_ref();
        let output = output.as_ref();
        Self::validate_output_path(output)?;
        self.check_version()?;
        Self::ensure_parent_dir(output)?;
        let mut root_inputs = Vec::with_capacity(dependencies.len() + 1);
        root_inputs.push(template);
        root_inputs.extend_from_slice(dependencies);
        let root = Self::compilation_root(&root_inputs)?;
        let mut command = Command::new(&self.typst_executable);
        command
            .arg("c")
            .arg(template)
            .arg(output)
            .arg("--root")
            .arg(root);
        self.run_typst(
            &mut command,
            &format!("compiling template {}", output.display()),
            None,
        )
    }

    fn validate_modules(modules: &[TypstModule]) -> Result<()> {
        let mut aliases = std::collections::HashSet::with_capacity(modules.len());
        for module in modules {
            TypstModule::validate_alias(&module.alias)?;
            if !aliases.insert(&module.alias) {
                bail!("duplicate Typst module alias {:?}", module.alias);
            }
        }
        Ok(())
    }

    fn entrypoint_source(
        request: &TypstRenderRequest,
        template: &Path,
        module_files: &[Option<PathBuf>],
        dot_path: &Path,
        root: &Path,
    ) -> Result<String> {
        let template = Self::typst_root_path(template, root)?;
        let dot_path = Self::typst_root_path(dot_path, root)?;
        let mut source = format!("#import {template} as _clinnet_template\n");
        for (module, file) in request.modules.iter().zip(module_files) {
            let module_source = match (&module.source, file) {
                (TypstModuleSource::File(_), Some(path)) => Self::typst_root_path(path, root)?,
                (TypstModuleSource::Package(package), None) => Self::typst_string(package),
                _ => bail!("internal Typst module source mismatch"),
            };
            source.push_str(&format!("#import {module_source} as {}\n", module.alias));
        }
        source.push_str("\n#let _clinnet_config = {\n  let value = (");
        source.push_str(&request.config.source()?);
        source.push_str(
            ")\n  if type(value) != dictionary {\n    panic(\"Clinnet render config must be a dictionary\")\n  }\n  value + (data-path: ",
        );
        source.push_str(&dot_path);
        source.push_str(",)\n}\n\n#_clinnet_template.render(_clinnet_config)\n");
        Ok(source)
    }

    fn typst_root_path(path: &Path, root: &Path) -> Result<String> {
        let path = Self::canonicalize_existing(path)?;
        let relative = path.strip_prefix(root).with_context(|| {
            format!(
                "Typst input {} is outside project root {}",
                path.display(),
                root.display()
            )
        })?;
        Ok(Self::typst_string(&format!(
            "/{}",
            relative.to_string_lossy().replace('\\', "/")
        )))
    }

    fn typst_string(value: &str) -> String {
        format!(
            "\"{}\"",
            value
                .replace('\\', "\\\\")
                .replace('"', "\\\"")
                .replace('\n', "\\n")
                .replace('\r', "\\r")
        )
    }

    fn require_typst_version(version_line: &str) -> Result<()> {
        let version = version_line
            .split_whitespace()
            .nth(1)
            .ok_or_else(|| eyre!("could not parse `typst --version` output: {version_line}"))?;
        let (numeric, prerelease) = version
            .split_once('-')
            .map_or((version, None), |(value, suffix)| (value, Some(suffix)));
        let mut parts = numeric.split('.');
        let parsed = (
            parts.next().and_then(|part| part.parse().ok()),
            parts.next().and_then(|part| part.parse().ok()),
            parts.next().and_then(|part| part.parse().ok()),
        );
        let (Some(major), Some(minor), Some(patch)) = parsed else {
            bail!("could not parse Typst version from: {version_line}");
        };
        if (major, minor, patch) < MINIMUM_TYPST_VERSION
            || ((major, minor, patch) == MINIMUM_TYPST_VERSION && prerelease.is_some())
        {
            bail!("linnet requires Typst 0.15.0 or newer, found {version}");
        }
        Ok(())
    }

    /// Find the narrowest existing directory that contains every compilation input.
    fn compilation_root(paths: &[&Path]) -> Result<PathBuf> {
        let mut directories = Vec::with_capacity(paths.len());
        for path in paths {
            let canonical = Self::canonicalize_existing(path)
                .with_context(|| format!("failed to resolve Typst input {}", path.display()))?;
            directories.push(if canonical.is_dir() {
                canonical
            } else {
                canonical
                    .parent()
                    .ok_or_else(|| eyre!("Typst input {} has no parent directory", path.display()))?
                    .to_path_buf()
            });
        }

        let first = directories
            .first()
            .ok_or_else(|| eyre!("cannot choose a Typst project root without inputs"))?;
        first
            .ancestors()
            .find(|ancestor| directories.iter().all(|path| path.starts_with(ancestor)))
            .map(Path::to_path_buf)
            .ok_or_else(|| eyre!("Typst inputs do not share a project root"))
    }

    fn validate_output_path(output: &Path) -> Result<()> {
        let extension = output
            .extension()
            .and_then(OsStr::to_str)
            .map(str::to_ascii_lowercase);
        if matches!(extension.as_deref(), Some("pdf" | "svg" | "png")) {
            Ok(())
        } else {
            bail!(
                "unsupported Typst output {}; expected a .pdf, .svg, or .png suffix",
                output.display()
            )
        }
    }

    /// Canonicalize a path, returning an error if it does not exist.
    fn canonicalize_existing(path: &Path) -> Result<PathBuf> {
        Ok(fs::canonicalize(path)?)
    }

    /// Resolve a template path, writing the embedded default if the requested file does not exist.
    fn resolve_template(&self, requested: &Path, kind: TemplateKind) -> Result<PathBuf> {
        match fs::canonicalize(requested) {
            Ok(path) => return Ok(path),
            Err(err) if err.kind() == std::io::ErrorKind::NotFound => {}
            Err(err) => {
                return Err(err)
                    .with_context(|| format!("failed to read template {}", requested.display()));
            }
        }

        let mut target = requested.to_path_buf();
        if target.file_name().is_none() {
            target = target.join(kind.file_name());
        }

        let templates_dir = self.build_dir.join(DEFAULT_TEMPLATE_SUBDIR);
        if !target.starts_with(&templates_dir) {
            bail!(
                "template {} not found and automatic creation is limited to {}",
                requested.display(),
                templates_dir.display()
            );
        }

        if target.exists() {
            return Ok(target);
        }
        Self::ensure_parent_dir(&target)?;
        let contents = kind.embedded_bytes()?;
        fs::write(&target, contents.as_ref())
            .with_context(|| format!("failed to write default template {}", target.display()))?;
        Ok(target)
    }

    fn write_embedded_assets<E: RustEmbed>(
        root: &Path,
        description: &str,
        overwrite: bool,
    ) -> Result<()> {
        for path in E::iter() {
            let target = root.join(path.as_ref());
            let contents = E::get(path.as_ref())
                .ok_or_else(|| eyre!("embedded {description} {} is missing", path))?;
            if target.exists()
                && (!overwrite
                    || fs::read(&target).is_ok_and(|existing| existing == contents.data.as_ref()))
            {
                continue;
            }
            Self::ensure_parent_dir(&target)?;
            fs::write(&target, contents.data.as_ref()).with_context(|| {
                format!(
                    "failed to write embedded {description} {}",
                    target.display()
                )
            })?;
        }
        Ok(())
    }

    fn ensure_parent_dir(path: &Path) -> Result<()> {
        if let Some(parent) = path.parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create directory {}", parent.display()))?;
        }
        Ok(())
    }

    fn run_typst(
        &self,
        command: &mut Command,
        description: &str,
        context: Option<RenderContext<'_>>,
    ) -> Result<()> {
        let output = command
            .output()
            .with_context(|| format!("failed to run typst while {description}"))?;

        if !output.status.success() {
            let stdout = String::from_utf8_lossy(&output.stdout);
            let stderr = String::from_utf8_lossy(&output.stderr);

            // Build detailed error message
            let mut error_msg = format!("typst failed while {description}");

            // Show the command that was executed
            let cmd_str = format!("{:?}", command);
            error_msg.push_str(&format!("\n\nCommand executed:\n{cmd_str}"));

            // Add context-specific diagnostics
            if let Some(context) = context {
                error_msg.push_str("\n\nDiagnostic information:");
                error_msg.push_str(&format!("\n  - Input file: {}", context.input.display()));
                error_msg.push_str(&format!("\n  - Template: {}", context.template.display()));
                error_msg.push_str(&format!("\n  - Root directory: {}", context.root.display()));
                error_msg.push_str(&format!("\n  - Output path: {}", context.output.display()));

                // Check if files exist
                if !context.input.exists() {
                    error_msg.push_str("\n  - ERROR: Input file does not exist!");
                }
                if !context.template.exists() {
                    error_msg.push_str("\n  - ERROR: Template file does not exist!");
                }
                if !context.root.exists() {
                    error_msg.push_str("\n  - ERROR: Root directory does not exist!");
                }

                // Check if root contains the template
                if let Ok(template_canonical) = context.template.canonicalize()
                    && let Ok(root_canonical) = context.root.canonicalize()
                    && !template_canonical.starts_with(&root_canonical)
                {
                    error_msg.push_str("\n  - WARNING: Template is outside root directory");
                    error_msg.push_str(
                        "\n    This may cause 'source file must be contained in project root' errors",
                    );
                    error_msg.push_str(
                        "\n    Consider using --root with a parent directory of both template and data files",
                    );
                }
            }

            // Add common solutions
            if stderr.contains("source file must be contained in project root") {
                error_msg.push_str("\n\nCommon solutions for 'project root' errors:");
                error_msg.push_str("\n  1. Ensure --root points to a directory that contains both templates and data files");
                error_msg.push_str("\n  2. Use absolute paths or adjust the working directory");
                error_msg.push_str(
                    "\n  3. Check that template and input files are in the same directory tree",
                );
            }

            if stderr.contains("not found") || stderr.contains("No such file") {
                error_msg.push_str(
                    "\n\nFile not found - check that all paths are correct and files exist",
                );
            }

            // Add the actual typst output
            error_msg.push_str(&format!("\n\nTypst stdout:\n{}", stdout.trim()));
            error_msg.push_str(&format!("\n\nTypst stderr:\n{}", stderr.trim()));

            bail!("{error_msg}");
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn value_dict(fields: impl IntoIterator<Item = (impl Into<String>, TypstValue)>) -> TypstValue {
        TypstValue::Dictionary(
            fields
                .into_iter()
                .map(|(key, value)| (key.into(), value))
                .collect(),
        )
    }

    fn empty_config() -> TypstConfig {
        TypstConfig::new(std::collections::BTreeMap::new()).unwrap()
    }

    #[test]
    fn requires_typst_0_15_or_newer() {
        assert!(TypstRenderer::require_typst_version("typst 0.15.0").is_ok());
        assert!(TypstRenderer::require_typst_version("typst 0.15.1").is_ok());
        assert!(TypstRenderer::require_typst_version("typst 0.15.1-rc1").is_ok());
        assert!(TypstRenderer::require_typst_version("typst 1.0.0").is_ok());
        assert!(TypstRenderer::require_typst_version("typst 0.15.0-rc1").is_err());
        assert!(TypstRenderer::require_typst_version("typst 0.14.2").is_err());
    }

    #[test]
    fn prepared_render_owns_its_ephemeral_inputs() {
        let base = tempfile::tempdir().unwrap();
        let template = base.path().join("template.typ");
        fs::write(&template, "#let render(config) = [ok]").unwrap();
        let request = TypstRenderRequest::new("digraph example { a -> b }", empty_config())
            .template(&template);
        let renderer = TypstRenderer::new(base.path());

        let entrypoint = {
            let prepared = renderer.prepare(&request).unwrap();
            assert!(prepared.entrypoint().starts_with(prepared.root()));
            assert!(prepared.root().is_dir());
            let source = fs::read_to_string(prepared.entrypoint()).unwrap();
            assert!(source.contains("data-path:"));
            assert!(source.contains("_clinnet_template.render(_clinnet_config)"));
            prepared.entrypoint().to_owned()
        };
        assert!(!entrypoint.exists());
    }

    #[test]
    fn rejects_unsupported_output_suffix() {
        let base = std::env::temp_dir().join(format!(
            "clinnet-output-suffix-{}-{}",
            std::process::id(),
            SCRATCH_ID.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&base).unwrap();
        let renderer = TypstRenderer::new(&base)
            .typst_executable("/__clinnet_test_missing_typst_executable__");
        let request = TypstRenderRequest::new("digraph example { a -> b }", empty_config());
        let error = renderer
            .render(&request, base.join("diagram.jpg"))
            .unwrap_err();
        assert!(error.to_string().contains(".pdf, .svg, or .png"));
        assert!(TypstRenderer::validate_output_path(Path::new("diagram.SVG")).is_ok());
        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn reports_command_spawn_failure_without_typst() {
        let base = std::env::temp_dir().join(format!(
            "clinnet-command-failure-{}-{}",
            std::process::id(),
            SCRATCH_ID.fetch_add(1, Ordering::Relaxed)
        ));
        let renderer = TypstRenderer::new(&base)
            .typst_executable("/__clinnet_test_missing_typst_executable__");
        let request = TypstRenderRequest::new("digraph example { a -> b }", empty_config());
        let error = renderer
            .render(&request, base.join("diagram.svg"))
            .unwrap_err();
        assert!(error.to_string().contains("failed to run"));
        fs::remove_dir_all(base).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn renders_dot_to_svg_with_configured_typst_executable() {
        use std::os::unix::fs::PermissionsExt;

        let base = std::env::temp_dir().join(format!(
            "clinnet-svg-render-{}-{}",
            std::process::id(),
            SCRATCH_ID.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&base).unwrap();
        let executable = base.join("typst");
        let captured = base.join("entrypoint.typ");
        fs::write(
            &executable,
            format!(
                "#!/bin/sh\nif [ \"$1\" = \"--version\" ]; then\n  echo 'typst 0.15.0'\nelse\n  [ \"$1\" = \"compile\" ] || exit 2\n  for arg in \"$@\"; do\n    [ \"$arg\" = \"--input\" ] && exit 3\n  done\n  cp \"$2\" '{}'\n  printf '<svg>rendered</svg>' > \"$3\"\nfi\n",
                captured.display()
            ),
        )
        .unwrap();
        let mut permissions = fs::metadata(&executable).unwrap().permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&executable, permissions).unwrap();

        let config = TypstConfig::new(std::collections::BTreeMap::from([
            ("title".to_owned(), TypstValue::String("diagram".to_owned())),
            (
                "payload".to_owned(),
                TypstValue::String("x".repeat(100_000)),
            ),
        ]))
        .unwrap();
        let request = TypstRenderRequest::new("digraph example { a -> b }", config)
            .module(TypstModule::package("cetz", "@preview/cetz:0.5.1").unwrap());
        let svg = TypstRenderer::new(&base)
            .typst_executable(&executable)
            .to_svg(&request)
            .unwrap();
        assert_eq!(svg, "<svg>rendered</svg>");
        let entrypoint = fs::read_to_string(captured).unwrap();
        assert!(entrypoint.contains("_clinnet_template.render(_clinnet_config)"));
        assert!(entrypoint.contains("#import \"@preview/cetz:0.5.1\" as cetz"));
        assert!(entrypoint.contains("data-path:"));
        assert!(entrypoint.len() > 100_000);

        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn embedded_packages_include_nested_sources() {
        assert!(EmbeddedLinnestPackage::get("src/impl/draw.typ").is_some());
        assert!(EmbeddedKurvstPackage::get("src/impl.typ").is_some());
    }

    #[test]
    fn validates_static_module_aliases() {
        assert!(TypstModule::new("style_0", "style.typ").is_ok());
        assert!(TypstModule::new("bad-alias", "style.typ").is_err());
        assert!(TypstModule::new("_clinnet_template", "style.typ").is_err());
        assert!(TypstModule::package("cetz", "@preview/cetz:0.5.1").is_ok());
        assert!(TypstModule::package("cetz", "cetz:0.5.1").is_err());
        assert!(TypstModule::package("cetz", "@preview/cetz:\"bad\"").is_err());
    }

    #[test]
    fn rejects_unimported_config_module_symbols() {
        let base = std::env::temp_dir().join(format!(
            "clinnet-config-module-{}-{}",
            std::process::id(),
            SCRATCH_ID.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&base).unwrap();
        let config = TypstConfig::new(std::collections::BTreeMap::from([(
            "style".to_owned(),
            TypstValue::ModuleSymbol {
                alias: "styles".to_owned(),
                path: vec!["edge_style".to_owned()],
            },
        )]))
        .unwrap();
        let request = TypstRenderRequest::new("digraph example { a -> b }", config);
        let error = TypstRenderer::new(&base)
            .render(&request, base.join("diagram.svg"))
            .unwrap_err();
        assert!(error.to_string().contains("missing Typst module alias"));
        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn default_staging_preserves_user_template_bundle() {
        let base = std::env::temp_dir().join(format!(
            "clinnet-private-templates-{}-{}",
            std::process::id(),
            SCRATCH_ID.fetch_add(1, Ordering::Relaxed)
        ));
        let user_template = base.join("templates/figure.typ");
        fs::create_dir_all(user_template.parent().unwrap()).unwrap();
        fs::write(&user_template, "#let render(config) = [Gamma]").unwrap();

        TypstRenderer::new(&base).stage_default_assets().unwrap();

        assert_eq!(
            fs::read_to_string(&user_template).unwrap(),
            "#let render(config) = [Gamma]"
        );
        assert!(
            base.join(DEFAULT_TEMPLATE_SUBDIR)
                .join("figure.typ")
                .is_file()
        );
        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn default_template_renders_generic_config_and_static_subgraphs_with_typst() {
        let typst = std::env::var_os("TYPST_TEST_EXECUTABLE")
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("typst"));
        let Ok(version) = Command::new(&typst).arg("--version").output() else {
            return;
        };
        if !version.status.success()
            || TypstRenderer::require_typst_version(&String::from_utf8_lossy(&version.stdout))
                .is_err()
        {
            return;
        }

        let base = std::env::temp_dir().join(format!(
            "clinnet-real-typst-{}-{}",
            std::process::id(),
            SCRATCH_ID.fetch_add(1, Ordering::Relaxed)
        ));
        let renderer = TypstRenderer::new(&base).typst_executable(typst);
        let dot = r#"digraph generic {
            0 [id=0]; 1 [id=1];
            ext0 [style=invis]; ext0 -> 0:0 [id=0];
            0:1 -> 1:2 [id=1];
            ext1 [style=invis]; 1:3 -> ext1 [id=2];
        }"#;
        let elements = value_dict([
            ("graph", TypstValue::None),
            (
                "nodes",
                TypstValue::Array(vec![
                    value_dict([("label", TypstValue::String("alpha".to_owned()))]),
                    value_dict([("label", TypstValue::String("beta".to_owned()))]),
                ]),
            ),
            (
                "edges",
                TypstValue::Array(vec![
                    value_dict([("label", TypstValue::String("input".to_owned()))]),
                    value_dict([
                        ("label", TypstValue::String("relation".to_owned())),
                        ("minimum-length", TypstValue::Integer(2)),
                    ]),
                    value_dict([("label", TypstValue::String("output".to_owned()))]),
                ]),
            ),
            ("hedges", TypstValue::Array(Vec::new())),
        ]);
        let all_hedges = TypstValue::Array(vec![
            TypstValue::Bool(true),
            TypstValue::Bool(true),
            TypstValue::Bool(true),
            TypstValue::Bool(true),
        ]);
        let layouts = TypstValue::Array(vec![
            value_dict([
                (
                    "layout-algo",
                    TypstValue::String("stable-layered".to_owned()),
                ),
                (
                    "rank-same",
                    TypstValue::Array(vec![TypstValue::Array(vec![
                        TypstValue::Integer(0),
                        TypstValue::Integer(1),
                    ])]),
                ),
                ("subgraph", all_hedges.clone()),
            ]),
            value_dict([
                ("layout-algo", TypstValue::String("force".to_owned())),
                ("layout-nodes", TypstValue::String("fixed".to_owned())),
                ("steps", TypstValue::Integer(2)),
                ("subgraph", all_hedges),
            ]),
        ]);
        let config = TypstConfig::new(std::collections::BTreeMap::from([
            ("layouts".to_owned(), layouts),
            (
                "draw".to_owned(),
                value_dict([(
                    "subgraph",
                    TypstValue::Array(vec![
                        TypstValue::Array(vec![
                            TypstValue::Bool(true),
                            TypstValue::Bool(false),
                            TypstValue::Bool(true),
                            TypstValue::Bool(false),
                        ]),
                        TypstValue::Array(vec![
                            TypstValue::Bool(false),
                            TypstValue::Bool(true),
                            TypstValue::Bool(false),
                            TypstValue::Bool(true),
                        ]),
                    ]),
                )]),
            ),
            ("elements".to_owned(), elements),
            (
                "options".to_owned(),
                value_dict([("application", TypstValue::String("example".to_owned()))]),
            ),
        ]))
        .unwrap();

        let svg = renderer
            .to_svg(&TypstRenderRequest::new(dot, config))
            .unwrap();
        assert!(svg.contains("<svg"));

        let templates = base.join(DEFAULT_TEMPLATE_SUBDIR);
        let assertions = templates.join("generic-assertions.typ");
        fs::write(
            &assertions,
            r#"#import "layout.typ": (
  _attach-elements,
  _apply-element-layout-constraints,
  _apply-label-offsets,
  _edge-rank-same,
  _element-edge-style,
  _layout-pass,
)
#import "crates/linnest/typst/src/lib.typ": graph

#let source-drawing = (statement: "source-statement", compass: "e")
#source-drawing.insert("port-label", "source-port")
#let sink-drawing = (statement: "sink-statement", compass: "w")
#sink-drawing.insert("port-label", "sink-port")
#let attached = _attach-elements(
  graph.build({
    graph.node(<source>)
    graph.node(<sink>)
    graph.edge(graph.source(<source>, id: 0), graph.sink(<sink>, id: 1))
  }),
  (hedges: (source-drawing, sink-drawing)),
)
#let attached-edge = graph.edges(attached).first()
#assert(attached-edge.source.statement == "source-statement")
#assert(attached-edge.source.at("port-label") == "source-port")
#assert(attached-edge.source.compass == "e")
#assert(attached-edge.sink.statement == "sink-statement")
#assert(attached-edge.sink.at("port-label") == "sink-port")
#assert(attached-edge.sink.compass == "w")

#let constrained-base = graph.build({
  graph.node(<left>, pos: graph.pos(x: 0, y: 0, mode: "pin"))
  graph.node(<right>, pos: graph.pos(x: 2, y: 0, mode: "pin"))
  graph.edge(
    graph.source(<left>),
    graph.sink(<right>),
    pos: graph.pos(x: 1, y: 0, mode: "pin"),
    label-pos: graph.pos(x: 1, y: 0, mode: "pin"),
  )
})
#let constrained-base = graph.map(
  constrained-base,
  node: node => if node.node == 0 {
    (data: (rank: 2, minimum-size: 3))
  } else {
    (data: (maximum-size: 4))
  },
  edge: edge => (data: (
    minimum-length: 3,
    same-rank: true,
    label-offset: 0.5,
  )),
)
#let constrained = graph.style(
  constrained-base,
  node-label: none,
  node-style: node => if node.vid == 1 { (radius: 5) } else { (:) },
)
#let constrained = _apply-element-layout-constraints(constrained)
#let constrained-nodes = graph.nodes(constrained)
#let constrained-edge = graph.edges(constrained).first()
#assert(str(constrained-nodes.first().statements.at("layout-rank")) == "2")
#assert(float(str(constrained-nodes.first().statements.at("layout-width"))) == 3)
#assert(float(str(constrained-nodes.first().statements.at("layout-height"))) == 3)
#assert(float(str(constrained-nodes.last().statements.at("layout-width"))) == 4)
#assert(float(str(constrained-nodes.last().statements.at("layout-height"))) == 4)
#assert(str(constrained-edge.statements.at("minlen")) == "3")
#assert(_edge-rank-same(constrained) == ((0, 1),))
#let constrained-pass = _layout-pass(constrained, (layout-algo: "stable-layered"))
#assert(constrained-pass.at("rank-same").len() == 1)
#assert(type(constrained-pass.at("rank-same").first()) == bytes)
#let offset = _apply-label-offsets(constrained)
#let offset-pos = graph.edges(offset).first().label-pos
#assert(calc.abs(offset-pos.x - 1) < 1e-6)
#assert(calc.abs(offset-pos.y - 0.5) < 1e-6)

#let decorated = _element-edge-style(
  (
    (stroke: (paint: black, thickness: 0.5pt), pattern: "coil"),
    (offset: 0.2, stroke: (paint: black, thickness: 0.4pt), mark: (end: "straight")),
  ),
  (data: (decoration: "wave",)),
)
#assert(decorated.first().pattern == "wave")
#assert(decorated.last().offset == 0.2)
#assert(decorated.last().stroke.paint == black)
#assert(decorated.last().mark.end == "straight")

#let explicit-partial = graph.parse(
  "digraph explicit { 0 [id=0]; ext [style=invis]; ext -> 0 [id=0, pos=\"y:37!\"]; }",
).first()
#let explicit-edge = graph.edges(explicit-partial).first()
#assert(explicit-edge.at("pos-x-set") == false)
#assert(explicit-edge.at("pos-y-set") == true)
#let explicit-roundtrip = graph.parse(graph.dot(explicit-partial)).first()
#let explicit-roundtrip-edge = graph.edges(explicit-roundtrip).first()
#assert(explicit-roundtrip-edge.at("pos-x-set") == false)
#assert(explicit-roundtrip-edge.at("pos-y-set") == true)
[ok]
"#,
        )
        .unwrap();
        renderer
            .compile_template(&assertions, base.join("generic-assertions.pdf"), &[])
            .unwrap();

        let gamma_core = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../assets/embedded/drawing/templates/layout-core.typ")
            .canonicalize()
            .unwrap();
        let physics_module = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../assets/embedded/drawing/templates/physics-edge-style.typ")
            .canonicalize()
            .unwrap();
        let graph_module = templates
            .join("crates/linnest/typst/src/graph.typ")
            .canonicalize()
            .unwrap();
        let subgraph_module = templates
            .join("crates/linnest/typst/src/subgraph.typ")
            .canonicalize()
            .unwrap();
        let gamma_assertions = base.join("gamma-mode-assertions.typ");
        let gamma_source = r#"#import "@CORE@": (
  _attach-elements,
  _apply-element-layout-constraints,
  _apply-label-offsets,
  _edge-label,
  _edge-rank-same,
  _element-edge-style,
  _external-label-options,
  _layout-pass,
  _momentum-edge-label,
  _particle-label,
  _resolved-mode,
  autogen-external-edge-fields,
)
#import "@PHYSICS@" as physics
#import "@GRAPH@" as graph
#import "@SUBGRAPH@" as subgraph
#let fake-style = (map: ("fake": (label: [particle])))
#let explicit-map = ("fake": (label: [override]))
#let particle-map = physics.default-map + fake-style.map + explicit-map
#let edge = (fields: (particle: "fake"), data: (:), edge: 2)
#let gamma-source-drawing = (statement: "gamma-source", compass: "e")
#gamma-source-drawing.insert("port-label", "gamma-port")
#let gamma-attached = _attach-elements(
  graph.build({
    graph.node(<gamma-a>)
    graph.node(<gamma-b>)
    graph.edge(graph.source(<gamma-a>, id: 0), graph.sink(<gamma-b>, id: 1))
  }),
  (hedges: (gamma-source-drawing, (:))),
  graph: graph,
)
#let gamma-source = graph.edges(gamma-attached).first().source
#assert(gamma-source.statement == "gamma-source")
#assert(gamma-source.at("port-label") == "gamma-port")
#assert(gamma-source.compass == "e")
#assert(_external-label-options((show-particle: false, show-edge-index: false)) == (
  include-particle: false,
  include-index: false,
))
#assert(_particle-label(
  edge,
  2,
  "plain",
  physics: physics,
  particle-map: particle-map,
  include-particle: false,
  include-index: false,
) == none)
#let index-only = _particle-label(
  edge,
  2,
  "plain",
  physics: physics,
  particle-map: particle-map,
  include-particle: false,
  include-index: true,
)
#assert(type(index-only) == content)
#assert(repr(index-only).contains("p") and repr(index-only).contains("2"))
#let particle-only = _particle-label(
  edge,
  2,
  "plain",
  physics: physics,
  particle-map: particle-map,
  include-particle: true,
  include-index: false,
)
#assert(type(particle-only) == content)
#assert(repr(particle-only).contains("override"))

#let generated-base = graph.build({
  graph.node(<external>)
  graph.edge(graph.sink(<external>), particle: "fake")
})
#let generated = autogen-external-edge-fields(
  generated-base,
  graph: graph,
  physics: physics,
  particle-map: particle-map,
)
#let generated-edge = graph.edges(generated).first()
#assert(generated-edge.data.keys().contains("mode-label"))
#assert(not generated-edge.data.keys().contains("label"))
#assert(repr(_momentum-edge-label(
  generated-edge,
  "plain",
  (:),
)) == repr(generated-edge.data.at("mode-label")))
#assert(repr(_edge-label([configured], generated-edge)) == repr([configured]))
#assert(_edge-label(none, generated-edge) == none)
#let element-label = graph.map(
  generated,
  edge: edge => (data: edge.data + (label: [element],)),
)
#assert(repr(_edge-label(
  [configured],
  graph.edges(element-label).first(),
)) == repr([element]))

// Match the generated GammaLoop wrapper contract: the explicit map is spread
// once, then merged above the generated and model-neutral maps.
#let generated-map = ("fake": (
  source: (stroke: (paint: blue, thickness: 0.55pt)),
  sink: (stroke: (paint: blue, thickness: 0.55pt)),
))
#let generated-style(map: (:), typst-fields: "plain", ..options) = physics.style(
  map: physics.default-map + generated-map + map,
  typst-fields: typst-fields,
  ..options.named(),
)
#let generated-source-style(edge, typst-fields: "plain", ..options) = {
  let callbacks = generated-style(
    typst-fields: typst-fields,
    ..options.named(),
  )
  (callbacks.source-style)(edge)
}
#let styled-edge = (
  particle: "fake",
  source-half-edge: (hedge: 0),
  sink-half-edge: (hedge: 1),
  orientation: "default",
  data: (:),
)
#let explicit-source = generated-source-style(styled-edge, map: ("fake": (
  source: (stroke: (paint: red, thickness: 0.55pt)),
  sink: (stroke: (paint: red, thickness: 0.55pt)),
),))
#assert(explicit-source.stroke.paint == red)

#let constrained-base = graph.build({
  graph.node(<left>, pos: graph.pos(x: 0, y: 0, mode: "pin"))
  graph.node(<right>, pos: graph.pos(x: 2, y: 0, mode: "pin"))
  graph.edge(
    graph.source(<left>),
    graph.sink(<right>),
    pos: graph.pos(x: 1, y: 0, mode: "pin"),
    label-pos: graph.pos(x: 1, y: 0, mode: "pin"),
  )
})
#let constrained-base = graph.map(
  constrained-base,
  node: node => if node.node == 0 {
    (data: (rank: 2, minimum-size: 3))
  } else {
    (data: (maximum-size: 4))
  },
  edge: edge => (data: (
    minimum-length: 3,
    same-rank: true,
    label-offset: 0.5,
  )),
)
#let constrained = graph.style(
  constrained-base,
  node-label: none,
  node-style: node => if node.vid == 1 { (radius: 5) } else { (:) },
)
#let constrained = _apply-element-layout-constraints(constrained, graph: graph)
#let constrained-nodes = graph.nodes(constrained)
#let constrained-edge = graph.edges(constrained).first()
#assert(str(constrained-nodes.first().statements.at("layout-rank")) == "2")
#assert(float(str(constrained-nodes.first().statements.at("layout-width"))) == 3)
#assert(float(str(constrained-nodes.first().statements.at("layout-height"))) == 3)
#assert(float(str(constrained-nodes.last().statements.at("layout-width"))) == 4)
#assert(float(str(constrained-nodes.last().statements.at("layout-height"))) == 4)
#assert(str(constrained-edge.statements.at("minlen")) == "3")
#assert(_edge-rank-same(constrained, graph: graph) == ((0, 1),))
#let constrained-pass = _layout-pass(
  constrained,
  (layout-algo: "stable-layered"),
  graph: graph,
  subgraph: subgraph,
)
#assert(constrained-pass.at("rank-same").len() == 1)
#assert(type(constrained-pass.at("rank-same").first()) == bytes)
#let offset = _apply-label-offsets(constrained, graph: graph)
#let offset-pos = graph.edges(offset).first().label-pos
#assert(calc.abs(offset-pos.x - 1) < 1e-6)
#assert(calc.abs(offset-pos.y - 0.5) < 1e-6)

#let decorated = _element-edge-style(
  (
    (stroke: (paint: black, thickness: 0.5pt), pattern: "coil"),
    (offset: 0.2, stroke: (paint: black, thickness: 0.4pt), mark: (end: "straight")),
  ),
  (data: (
    decoration: "wave",
    momentum-style: (offset: 0.4, stroke: (paint: red, thickness: 0.8pt)),
  )),
  momentum-arrows: true,
)
#assert(decorated.first().pattern == "wave")
#assert(decorated.last().offset == 0.4)
#assert(decorated.last().stroke.paint == red)
#assert(decorated.last().mark.end == "straight")
#let without-momentum = _element-edge-style(
  decorated,
  (data: (momentum-style: none)),
  momentum-arrows: true,
)
#assert(type(without-momentum) == dictionary)
#assert(without-momentum.pattern == "wave")
#assert(without-momentum.stroke.paint == black)

#let generic = graph.build({
  graph.node(<a>)
  graph.node(<b>)
  graph.edge(graph.source(<a>), graph.sink(<b>))
})
#let amplitude = graph.build({
  graph.node(<a>)
  graph.edge(graph.sink(<a>))
})
#let cross-section = graph.build({
  graph.node(<a>)
  graph.node(<b>)
  graph.edge(graph.sink(<a>), cut-id: "pair")
  graph.edge(graph.source(<b>), cut-id: "pair")
})
#assert(_resolved-mode(generic, graph: graph, auto-mode: true) == (
  amplitude: false,
  cross-section: false,
))
#assert(_resolved-mode(amplitude, graph: graph, auto-mode: true) == (
  amplitude: true,
  cross-section: false,
))
#assert(_resolved-mode(cross-section, graph: graph, auto-mode: true) == (
  amplitude: false,
  cross-section: true,
))
#assert(_resolved-mode(generic, graph: graph, amplitude-mode: true) == (
  amplitude: true,
  cross-section: false,
))
[ok]
"#
        .replace("@CORE@", &gamma_core.to_string_lossy().replace('\\', "/"))
        .replace(
            "@PHYSICS@",
            &physics_module.to_string_lossy().replace('\\', "/"),
        )
        .replace(
            "@GRAPH@",
            &graph_module.to_string_lossy().replace('\\', "/"),
        )
        .replace(
            "@SUBGRAPH@",
            &subgraph_module.to_string_lossy().replace('\\', "/"),
        );
        fs::write(&gamma_assertions, gamma_source).unwrap();
        renderer
            .compile_template(
                &gamma_assertions,
                base.join("gamma-mode-assertions.pdf"),
                &[
                    gamma_core.as_path(),
                    physics_module.as_path(),
                    graph_module.as_path(),
                    subgraph_module.as_path(),
                ],
            )
            .unwrap();

        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn compilation_root_contains_sibling_graph_and_template_trees() {
        let base = std::env::temp_dir().join(format!(
            "clinnet-compilation-root-{}-{}",
            std::process::id(),
            SCRATCH_ID.fetch_add(1, Ordering::Relaxed)
        ));
        let graph = base.join("graphs/example.dot");
        let template = base.join("build/.clinnet/templates/figure.typ");
        fs::create_dir_all(graph.parent().unwrap()).unwrap();
        fs::create_dir_all(template.parent().unwrap()).unwrap();
        fs::write(&graph, "digraph example { a -> b }").unwrap();
        fs::write(&template, "[figure]").unwrap();

        assert_eq!(
            TypstRenderer::compilation_root(&[&graph, &template]).unwrap(),
            fs::canonicalize(&base).unwrap()
        );

        fs::remove_dir_all(base).unwrap();
    }
}
