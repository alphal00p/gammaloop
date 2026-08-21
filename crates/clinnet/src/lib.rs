use std::borrow::Cow;
use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, Mutex, OnceLock};

#[cfg(test)]
use std::sync::atomic::{AtomicU64, Ordering};

use eyre::{Context, Result, bail, eyre};
use pathdiff::diff_paths;
use rust_embed::RustEmbed;

const TEMPLATE_SUBDIR: &str = "templates";
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

/// Compile Linnet DOT graphs through an external Typst executable.
///
/// The renderer stages the bundled templates and Typst packages beneath its
/// build directory when no custom figure template is selected. Builder options
/// apply to both file-backed and in-memory DOT renders.
#[derive(Clone, Debug)]
pub struct TypstRenderer {
    build_dir: PathBuf,
    typst_executable: PathBuf,
    template: Option<PathBuf>,
    title: Option<String>,
    inputs: Vec<(String, String)>,
    checked_version: Arc<OnceLock<String>>,
}

impl TypstRenderer {
    pub fn new(build_dir: impl Into<PathBuf>) -> Self {
        Self {
            build_dir: build_dir.into(),
            typst_executable: PathBuf::from("typst"),
            template: None,
            title: None,
            inputs: Vec::new(),
            checked_version: Arc::new(OnceLock::new()),
        }
    }

    /// Select a Typst executable other than the `typst` found in `PATH`.
    pub fn typst_executable(mut self, executable: impl Into<PathBuf>) -> Self {
        self.typst_executable = executable.into();
        self.checked_version = Arc::new(OnceLock::new());
        self
    }

    /// Use a custom figure template instead of the bundled generic template.
    pub fn template(mut self, template: impl Into<PathBuf>) -> Self {
        self.template = Some(template.into());
        self
    }

    /// Set the figure title exposed through `sys.inputs`.
    pub fn title(mut self, title: impl Into<String>) -> Self {
        self.title = Some(title.into());
        self
    }

    /// Add a value exposed through Typst's `sys.inputs`.
    pub fn input(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.inputs.push((key.into(), value.into()));
        self
    }

    /// Add values exposed through Typst's `sys.inputs`.
    pub fn inputs<K, V>(mut self, inputs: impl IntoIterator<Item = (K, V)>) -> Self
    where
        K: Into<String>,
        V: Into<String>,
    {
        self.inputs.extend(
            inputs
                .into_iter()
                .map(|(key, value)| (key.into(), value.into())),
        );
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

    /// Ensure every embedded default template dependency exists under `build/templates`.
    pub fn stage_default_assets(&self) -> Result<()> {
        let _staging = ASSET_STAGING_LOCK
            .lock()
            .map_err(|_| eyre!("embedded asset staging lock was poisoned"))?;
        let template_dir = self.build_dir.join(TEMPLATE_SUBDIR);
        // Ensure the default layout template exists unless the user already created one.
        Self::write_embedded_assets::<EmbeddedTemplates>(
            &template_dir,
            "template dependency",
            false,
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

    /// Render a DOT file through the selected or bundled figure template.
    pub fn render_file(&self, dot_path: impl AsRef<Path>, output: impl AsRef<Path>) -> Result<()> {
        let dot_path = dot_path.as_ref();
        let output = output.as_ref();
        Self::validate_output_path(output)?;
        self.check_version()?;
        Self::ensure_parent_dir(output)?;

        let template = if let Some(template) = &self.template {
            Self::canonicalize_existing(template)
                .with_context(|| format!("failed to read template {}", template.display()))?
        } else {
            self.stage_default_assets()?;
            self.resolve_figure_template(
                self.build_dir
                    .join(TEMPLATE_SUBDIR)
                    .join(TemplateKind::Figure.file_name()),
            )?
        };
        let dot_path = Self::canonicalize_existing(dot_path)
            .with_context(|| format!("failed to read DOT input {}", dot_path.display()))?;
        let template_dir = template
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        let data_rel = diff_paths(&dot_path, &template_dir).unwrap_or_else(|| dot_path.clone());
        let relative_input = data_rel.to_string_lossy().replace('\\', "/");
        let root = Self::compilation_root(&[&template, &dot_path])?;
        let title = self.title.clone().unwrap_or_else(|| {
            dot_path
                .file_stem()
                .map(|name| name.to_string_lossy().into_owned())
                .unwrap_or_else(|| "diagram".to_owned())
        });
        let mut command = Command::new(&self.typst_executable);
        command
            .arg("c")
            .arg(&template)
            .arg(output)
            .arg("--root")
            .arg(&root)
            .arg("--input")
            .arg(format!("data-path={relative_input}"))
            .arg("--input")
            .arg(format!("data_path={relative_input}"))
            .arg("--input")
            .arg(format!("title={title}"));
        self.append_inputs(&mut command);

        self.run_typst(
            &mut command,
            &format!("building {}", dot_path.display()),
            Some(RenderContext {
                input: &dot_path,
                template: &template,
                root: &root,
                output,
            }),
        )
    }

    /// Render one in-memory DOT graph to a PDF, SVG, or PNG selected by `output`.
    pub fn render_dot(&self, dot: &str, output: impl AsRef<Path>) -> Result<()> {
        fs::create_dir_all(&self.build_dir).with_context(|| {
            format!(
                "failed to create renderer build directory {}",
                self.build_dir.display()
            )
        })?;
        let build_dir = Self::canonicalize_existing(&self.build_dir)?;
        let scratch_parent = if let Some(template) = &self.template {
            let template = Self::canonicalize_existing(template)
                .with_context(|| format!("failed to read template {}", template.display()))?;
            if Self::compilation_root(&[&template, &build_dir]).is_ok() {
                build_dir
            } else {
                template
                    .parent()
                    .ok_or_else(|| {
                        eyre!("template {} has no parent directory", template.display())
                    })?
                    .to_path_buf()
            }
        } else {
            build_dir
        };
        let scratch_dir = tempfile::Builder::new()
            .prefix(".clinnet-render-")
            .tempdir_in(&scratch_parent)
            .with_context(|| {
                format!(
                    "failed to create render scratch directory in {}",
                    scratch_parent.display()
                )
            })?;
        let source = scratch_dir.path().join("diagram.dot");
        fs::write(&source, dot)
            .with_context(|| format!("failed to stage DOT input {}", source.display()))?;
        self.render_file(&source, output)
    }

    /// Render one in-memory DOT graph and return the resulting SVG document.
    pub fn render_dot_to_svg(&self, dot: &str) -> Result<String> {
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
        self.render_dot(dot, &output)?;
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
        self.append_inputs(&mut command);
        self.run_typst(
            &mut command,
            &format!("compiling template {}", output.display()),
            None,
        )
    }

    fn append_inputs(&self, command: &mut Command) {
        // Add additional input arguments
        for (key, value) in &self.inputs {
            command.arg("--input").arg(format!("{key}={value}"));
        }
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

        let templates_dir = self.build_dir.join(TEMPLATE_SUBDIR);
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
    fn rejects_unsupported_output_suffix() {
        let base = std::env::temp_dir().join(format!(
            "clinnet-output-suffix-{}-{}",
            std::process::id(),
            SCRATCH_ID.fetch_add(1, Ordering::Relaxed)
        ));
        let renderer = TypstRenderer::new(&base)
            .typst_executable("/__clinnet_test_missing_typst_executable__");
        let error = renderer
            .render_dot("digraph example { a -> b }", base.join("diagram.jpg"))
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
        let error = renderer
            .render_dot("digraph example { a -> b }", base.join("diagram.svg"))
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
        fs::write(
            &executable,
            "#!/bin/sh\nif [ \"$1\" = \"--version\" ]; then\n  echo 'typst 0.15.0'\nelse\n  found=\n  for arg in \"$@\"; do\n    [ \"$arg\" = \"title=diagram\" ] && found=1\n  done\n  [ \"$found\" = 1 ] || exit 2\n  printf '<svg>rendered</svg>' > \"$3\"\nfi\n",
        )
        .unwrap();
        let mut permissions = fs::metadata(&executable).unwrap().permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(&executable, permissions).unwrap();

        let svg = TypstRenderer::new(&base)
            .typst_executable(&executable)
            .render_dot_to_svg("digraph example { a -> b }")
            .unwrap();
        assert_eq!(svg, "<svg>rendered</svg>");

        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn embedded_packages_include_nested_sources() {
        assert!(EmbeddedLinnestPackage::get("src/impl/draw.typ").is_some());
        assert!(EmbeddedKurvstPackage::get("src/impl.typ").is_some());
    }

    #[test]
    fn compilation_root_contains_sibling_graph_and_template_trees() {
        let base = std::env::temp_dir().join(format!(
            "clinnet-compilation-root-{}-{}",
            std::process::id(),
            SCRATCH_ID.fetch_add(1, Ordering::Relaxed)
        ));
        let graph = base.join("graphs/example.dot");
        let template = base.join("build/templates/figure.typ");
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
