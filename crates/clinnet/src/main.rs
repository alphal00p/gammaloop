use std::borrow::Cow;
use std::collections::{BTreeMap, HashMap};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Arc;
use std::time::Duration;

use blake3::Hasher;
use clap::{ArgAction, Parser};
use eyre::{Context, Result, bail, eyre};
use indicatif::{ProgressBar, ProgressStyle};
use parking_lot::Mutex as ParkingMutex;
use pathdiff::diff_paths;
use rayon::prelude::*;
use rust_embed::RustEmbed;
use serde::{Deserialize, Serialize};
use serde_json::{self, Map as JsonMap, Value as JsonValue};
use walkdir::WalkDir;

const TEMPLATE_SUBDIR: &str = "templates";

#[derive(RustEmbed)]
#[folder = "templates"]
struct EmbeddedTemplates;

#[derive(Copy, Clone)]
enum TemplateKind {
    Figure,
    Grid,
    Plugin,
    Layout,
}

impl TemplateKind {
    /// Return the embedded filename associated with this template kind.
    fn file_name(self) -> &'static str {
        match self {
            TemplateKind::Figure => "figure.typ",
            TemplateKind::Grid => "grid.typ",
            TemplateKind::Plugin => "linnest.wasm",
            TemplateKind::Layout => "layout.typ",
        }
    }

    fn embedded_bytes(self) -> Result<Cow<'static, [u8]>> {
        EmbeddedTemplates::get(self.file_name())
            .map(|file| file.data)
            .ok_or_else(|| eyre!("embedded template {} is missing", self.file_name()))
    }
}

fn main() {
    if let Err(err) = run() {
        eprintln!("error: {err}");
        for cause in err.chain().skip(1) {
            eprintln!("  caused by: {cause}");
        }
        std::process::exit(1);
    }
}

fn check_typst_version() -> Result<()> {
    let output = Command::new("typst")
        .arg("--version")
        .output()
        .context("failed to run 'typst --version'. Is typst installed and in PATH?")?;

    if !output.status.success() {
        bail!("typst --version failed");
    }

    let version_output = String::from_utf8_lossy(&output.stdout);
    let version_line = version_output.lines().next().unwrap_or("unknown");

    println!("Using {}", version_line);

    // Check for minimum version (0.11.0 introduced some breaking changes)
    if let Some(version_part) = version_line.split_whitespace().nth(1) {
        let version = version_part.split(' ').next().unwrap_or("0.0.0");
        if version.starts_with("0.") {
            let parts: Vec<&str> = version.split('.').collect();
            if parts.len() >= 2
                && let (Ok(major), Ok(minor)) = (parts[0].parse::<u32>(), parts[1].parse::<u32>())
                && major == 0
                && minor < 11
            {
                eprintln!(
                    "Warning: Typst version {} may not be fully supported. Consider upgrading to 0.11.0 or later.",
                    version
                );
            }
        }
    }

    Ok(())
}

#[derive(Parser, Debug)]
#[command(
    name = "linnet",
    about = "Incrementally render Typst figures for every .dot file and assemble them into a grid PDF."
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,

    /// Base directory that stores build artifacts.
    #[arg(long, value_name = "DIR", default_value = "build")]
    build_dir: PathBuf,
    /// Disable parallel processing (process figures sequentially).
    #[arg(long, action = ArgAction::SetTrue)]
    no_parallel: bool,
    /// Number of parallel jobs (0 = use all available cores).
    #[arg(long, short = 'j', value_name = "N", default_value_t = 0)]
    jobs: usize,
}

#[derive(Parser, Debug)]
enum Commands {
    /// Draw all figures and grid (default behavior)
    Draw(DrawArgs),
    /// Redraw a single figure using cached metadata
    RedrawFigure(RedrawFigureArgs),
    /// Redraw only the grid PDF using cached metadata
    RedrawGrid(RedrawGridArgs),
}

#[derive(Parser, Debug)]
struct DrawArgs {
    /// Directory to scan for .dot files.
    #[arg(value_name = "ROOT")]
    root: PathBuf,
    /// Shared Typst template used for every individual figure.
    #[arg(long, value_name = "FILE")]
    figure_template: Option<PathBuf>,
    /// Typst template used for the grid document.
    #[arg(long, value_name = "FILE")]
    grid_template: Option<PathBuf>,
    /// Directory for the generated figure PDFs (defaults to <build_dir>/figs).
    #[arg(long, value_name = "DIR")]
    figs_dir: Option<PathBuf>,
    /// Cache file storing hashes (defaults to <build_dir>/.cache/figures.json).
    #[arg(long, value_name = "FILE")]
    cache_file: Option<PathBuf>,
    /// Path for the generated fig-index.typ (defaults to the directory of grid_template).
    #[arg(long, value_name = "FILE")]
    fig_index: Option<PathBuf>,
    /// Extra files whose contents influence incremental rebuilds (e.g., style snippets).
    #[arg(long, value_name = "FILE", action = ArgAction::Append)]
    style: Vec<PathBuf>,
    #[command(flatten)]
    output_args: OutputArgs,
    #[command(flatten)]
    input_args: InputArgs,
}

#[derive(Parser, Debug)]
struct OutputArgs {
    /// Destination PDF for the final grid (defaults to <build_dir>/grid.pdf).
    #[arg(long, short = 'o', value_name = "FILE")]
    output_path: Option<PathBuf>,
}

#[derive(Parser, Debug)]
struct InputArgs {
    /// Add a string key-value pair visible through `sys.inputs` in Typst.
    #[arg(long, value_name = "key=value", action = ArgAction::Append, value_parser = parse_input_pair)]
    input: Vec<(String, String)>,
}

#[derive(Parser, Debug)]
struct RedrawFigureArgs {
    /// Path to the .dot file to rebuild
    #[arg(value_name = "DOT")]
    target: PathBuf,
    #[command(flatten)]
    input_args: InputArgs,
}

#[derive(Parser, Debug)]
struct RedrawGridArgs {
    #[command(flatten)]
    output_args: OutputArgs,
    #[command(flatten)]
    input_args: InputArgs,
}

#[derive(Clone, Serialize, Deserialize)]
struct FigurePlan {
    data_path: PathBuf,
    relative: PathBuf,
    output_path: PathBuf,
}

#[derive(Clone, Serialize, Deserialize)]
struct FigureRecord {
    output_path: PathBuf,
    relative: PathBuf,
}

#[derive(Clone, Serialize, Deserialize)]
struct RunMetadata {
    root: PathBuf,
    cwd: PathBuf,
    build_dir: PathBuf,
    figs_dir: PathBuf,
    figure_template: PathBuf,
    grid_template: PathBuf,
    grid_output: PathBuf,
    fig_index_path: PathBuf,
    cache_file: PathBuf,

    style_files: Vec<PathBuf>,
    input: Vec<(String, String)>,
    plans: Vec<FigurePlan>,
    records: Vec<FigureRecord>,
}

#[derive(Default)]
struct FolderNode {
    figures: Vec<FigureEntry>,
    children: BTreeMap<String, FolderNode>,
}

struct FigureEntry {
    path: String,
    relative: String,
    name: String,
}

/// Entry point for the `linnet` CLI: orchestrates full builds and targeted rebuilds.
fn run() -> Result<()> {
    let cli = Cli::parse();

    // Check typst availability and version early
    check_typst_version()?;

    // Configure rayon thread pool
    if let Some(threads) = if cli.jobs > 0 { Some(cli.jobs) } else { None } {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("failed to configure thread pool")?;
    }

    let cwd = std::env::current_dir().context("failed to resolve working directory")?;
    let build_dir = absolutize(&cwd, &cli.build_dir);
    let metadata_path = metadata_path(&build_dir);

    match &cli.command {
        Some(Commands::RedrawFigure(args)) => {
            let metadata = load_run_metadata(&metadata_path)?;
            let target_abs = absolutize(&cwd, &args.target);
            return rebuild_single_figure(&target_abs, &metadata, &args.input_args.input);
        }
        Some(Commands::RedrawGrid(args)) => {
            let metadata = load_run_metadata(&metadata_path)?;
            return rebuild_grid_with_overrides(&metadata, args, &cwd);
        }
        Some(Commands::Draw(_)) => {
            // Continue with draw logic using args
        }
        None => {
            bail!("Use subcommands: 'draw <ROOT>', 'redraw-figure <DOT>', or 'redraw-grid'");
        }
    }

    // Extract draw args for draw command only
    let draw_args = match &cli.command {
        Some(Commands::Draw(args)) => args,
        _ => unreachable!(), // We handle all cases above
    };

    let root = canonicalize_existing(&draw_args.root)
        .with_context(|| format!("failed to read root directory {}", draw_args.root.display()))?;
    let figs_dir = draw_args
        .figs_dir
        .as_ref()
        .map(|path| absolutize(&cwd, path))
        .unwrap_or_else(|| build_dir.join("figs"));
    let cache_file = draw_args
        .cache_file
        .as_ref()
        .map(|path| absolutize(&cwd, path))
        .unwrap_or_else(|| build_dir.join(".cache").join("figures.json"));
    let grid_output = draw_args
        .output_args
        .output_path
        .as_ref()
        .map(|path| absolutize(&cwd, path))
        .unwrap_or_else(|| build_dir.join("grid.pdf"));
    let requested_figure_template = draw_args
        .figure_template
        .as_ref()
        .map(|path| absolutize(&cwd, path))
        .unwrap_or_else(|| build_dir.join("templates").join("figure.typ"));
    let requested_grid_template = draw_args
        .grid_template
        .as_ref()
        .map(|path| absolutize(&cwd, path))
        .unwrap_or_else(|| build_dir.join("templates").join("grid.typ"));

    fs::create_dir_all(&build_dir)
        .with_context(|| format!("failed to create build directory {}", build_dir.display()))?;
    fs::create_dir_all(&figs_dir)
        .with_context(|| format!("failed to create figures directory {}", figs_dir.display()))?;
    if let Some(parent) = cache_file.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create cache directory {}", parent.display()))?;
    }

    let layout_asset = ensure_layout_asset(&build_dir)?;

    let mut style_files = Vec::new();
    for style in &draw_args.style {
        let canonical = canonicalize_existing(style)
            .with_context(|| format!("failed to read style file {}", style.display()))?;
        style_files.push(canonical);
    }
    let layout_style = canonicalize_existing(&layout_asset).with_context(|| {
        format!(
            "failed to resolve layout template {}",
            layout_asset.display()
        )
    })?;
    style_files.push(layout_style);
    style_files.sort();
    style_files.dedup();

    let figure_template =
        resolve_template(&requested_figure_template, TemplateKind::Figure, &build_dir)?;
    let grid_template = resolve_template(&requested_grid_template, TemplateKind::Grid, &build_dir)?;
    let fig_index_path = draw_args
        .fig_index
        .as_ref()
        .map(|path| absolutize(&cwd, path))
        .unwrap_or_else(|| {
            grid_template
                .parent()
                .unwrap_or(Path::new("."))
                .join("fig-index.typ")
        });
    let _plugin_path = ensure_plugin_asset(&build_dir)?;

    let dot_files = collect_dot_files(&root)?;
    if dot_files.is_empty() {
        println!("No .dot files found under {}", root.display());
        return Ok(());
    }

    let mut plans = Vec::new();
    for data_path in dot_files {
        let relative = diff_paths(&data_path, &root).ok_or_else(|| {
            eyre!(
                "failed to compute relative path for {}",
                data_path.display()
            )
        })?;
        let mut output_path = figs_dir.join(&relative);
        output_path.set_extension("pdf");
        if let Some(parent) = output_path.parent() {
            fs::create_dir_all(parent).with_context(|| {
                format!(
                    "failed to create output directory {} for {}",
                    parent.display(),
                    relative.display()
                )
            })?;
        }
        plans.push(FigurePlan {
            data_path,
            relative,
            output_path,
        });
    }

    plans.sort_by(|a, b| a.relative.cmp(&b.relative));

    let metadata = RunMetadata {
        root: root.clone(),
        cwd: cwd.clone(),
        build_dir: build_dir.clone(),
        figs_dir: figs_dir.clone(),
        figure_template: figure_template.clone(),
        grid_template: grid_template.clone(),
        grid_output: grid_output.clone(),
        fig_index_path: fig_index_path.clone(),
        cache_file: cache_file.clone(),

        style_files: style_files.clone(),
        input: draw_args.input_args.input.clone(),
        plans: plans.clone(),
        records: plans
            .iter()
            .map(|plan| FigureRecord {
                output_path: plan.output_path.clone(),
                relative: plan.relative.clone(),
            })
            .collect(),
    };
    save_run_metadata(&metadata_path, &metadata)?;

    let previous_cache = load_cache(&cache_file)?;
    let new_cache = Arc::new(ParkingMutex::new(BTreeMap::new()));
    let rebuilt = Arc::new(ParkingMutex::new(0usize));
    let reused = Arc::new(ParkingMutex::new(0usize));

    let progress = Arc::new(ProgressBar::new(plans.len() as u64));
    progress.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] {wide_bar:.cyan/blue} {pos}/{len} {msg}",
        )
        .unwrap()
        .progress_chars("#>-"),
    );
    progress.enable_steady_tick(Duration::from_millis(120));
    progress.set_message("Preparing figures...");

    // Process figures in parallel or sequential based on CLI flags
    let results: Result<Vec<_>> = if cli.no_parallel {
        // Sequential processing
        plans
            .iter()
            .map(|plan| {
                let key = path_key(&plan.relative);
                let hash = compute_hash(
                    &plan.data_path,
                    &figure_template,
                    &style_files,
                    &draw_args.input_args.input,
                )?;
                let needs_rebuild = !previous_cache
                    .get(&key)
                    .map(|old| old == &hash)
                    .unwrap_or(false);

                if needs_rebuild {
                    progress.set_message(format!("building {}", plan.relative.display()));
                    build_figure(plan, &figure_template, &root, &draw_args.input_args.input)?;
                    *rebuilt.lock() += 1;
                } else {
                    progress.set_message(format!("cached {}", plan.relative.display()));
                    *reused.lock() += 1;
                }

                new_cache.lock().insert(key.clone(), hash.clone());
                progress.inc(1);
                Ok((key, hash))
            })
            .collect()
    } else {
        // Parallel processing
        plans
            .par_iter()
            .map(|plan| {
                let key = path_key(&plan.relative);
                let hash = compute_hash(
                    &plan.data_path,
                    &figure_template,
                    &style_files,
                    &draw_args.input_args.input,
                )?;
                let needs_rebuild = !previous_cache
                    .get(&key)
                    .map(|old| old == &hash)
                    .unwrap_or(false);

                if needs_rebuild {
                    progress.set_message(format!("building {}", plan.relative.display()));
                    build_figure(plan, &figure_template, &root, &draw_args.input_args.input)?;
                    *rebuilt.lock() += 1;
                } else {
                    progress.set_message(format!("cached {}", plan.relative.display()));
                    *reused.lock() += 1;
                }

                new_cache.lock().insert(key.clone(), hash.clone());
                progress.inc(1);
                Ok((key, hash))
            })
            .collect()
    };

    // Check if any parallel processing failed
    results?;

    progress.finish_with_message("figures ready");

    let final_cache = Arc::try_unwrap(new_cache).unwrap().into_inner();
    let rebuilt = Arc::try_unwrap(rebuilt).unwrap().into_inner();
    let reused = Arc::try_unwrap(reused).unwrap().into_inner();

    save_cache(&cache_file, &final_cache)?;
    let removed = remove_stale_outputs(&previous_cache, &final_cache, &figs_dir)?;
    save_run_metadata(&metadata_path, &metadata)?;

    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::with_template("{spinner:.green} {msg}")
            .unwrap()
            .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    spinner.enable_steady_tick(Duration::from_millis(80));
    spinner.set_message("Writing fig-index & compiling grid...");
    run_grid_from_metadata(&metadata)?;
    spinner.finish_with_message(format!("Grid ready -> {}", grid_output.display()));

    println!(
        "figures: {} built, {} reused{}",
        rebuilt,
        reused,
        if removed.is_empty() {
            String::new()
        } else {
            format!(", {} removed", removed.len())
        }
    );
    println!("grid: {}", grid_output.display());

    Ok(())
}

/// Recursively gather every `.dot` file beneath `root`, sorted for deterministic runs.
fn collect_dot_files(root: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in WalkDir::new(root).into_iter() {
        let entry = entry?;
        if entry.file_type().is_file()
            && entry
                .path()
                .extension()
                .and_then(OsStr::to_str)
                .map(|ext| ext.eq_ignore_ascii_case("dot"))
                .unwrap_or(false)
        {
            files.push(entry.into_path());
        }
    }
    files.sort();
    Ok(files)
}

/// Redraw a single figure using cached metadata from the previous full draw.
fn rebuild_single_figure(
    target: &Path,
    metadata: &RunMetadata,
    current_inputs: &[(String, String)],
) -> Result<()> {
    let canonical = canonicalize_existing(target)
        .with_context(|| format!("failed to read figure source {}", target.display()))?;
    let plan = metadata
        .plans
        .iter()
        .find(|plan| plan.data_path == canonical)
        .ok_or_else(|| eyre!("{} is not part of the cached plan", target.display()))?;

    println!("Redrawing figure: {}", plan.relative.display());
    println!("  Input: {}", plan.data_path.display());
    println!("  Output: {}", plan.output_path.display());
    if !current_inputs.is_empty() {
        println!("  Using {} current input argument(s)", current_inputs.len());
    }

    build_figure(
        plan,
        &metadata.figure_template,
        &metadata.root,
        current_inputs,
    )?;
    let existing = load_cache(&metadata.cache_file)?;
    let mut cache: BTreeMap<String, String> = existing.into_iter().collect();
    let hash = compute_hash(
        &plan.data_path,
        &metadata.figure_template,
        &metadata.style_files,
        current_inputs,
    )?;
    cache.insert(path_key(&plan.relative), hash);
    save_cache(&metadata.cache_file, &cache)?;

    println!(
        "Figure redrawn successfully: {}",
        plan.output_path.display()
    );
    Ok(())
}

/// Redraw grid with CLI overrides for output path and inputs.
fn rebuild_grid_with_overrides(
    metadata: &RunMetadata,
    args: &RedrawGridArgs,
    cwd: &Path,
) -> Result<()> {
    // Create modified metadata with CLI overrides
    let mut modified_metadata = metadata.clone();

    // Override output path if specified
    if let Some(output_path) = args.output_args.output_path.as_ref() {
        modified_metadata.grid_output = absolutize(cwd, output_path);
    }

    // Override inputs if specified
    if !args.input_args.input.is_empty() {
        modified_metadata.input = args.input_args.input.clone();
    }

    println!("Redrawing grid with overrides...");
    println!("  Output: {}", modified_metadata.grid_output.display());
    println!(
        "  Using {} cached figure(s)",
        modified_metadata.records.len()
    );
    if !modified_metadata.input.is_empty() {
        println!("  Input arguments:");
        for (key, value) in &modified_metadata.input {
            println!("    {}={}", key, value);
        }
    }

    run_grid_from_metadata(&modified_metadata)?;

    println!(
        "Grid redrawn successfully: {}",
        modified_metadata.grid_output.display()
    );
    Ok(())
}

/// Write `fig-index.typ` and run Typst for the grid using cached metadata.
fn run_grid_from_metadata(metadata: &RunMetadata) -> Result<()> {
    let grid_dir = metadata
        .grid_template
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."));
    write_fig_index(&metadata.records, &metadata.fig_index_path, &grid_dir)?;
    compile_grid(
        &metadata.grid_template,
        &metadata.grid_output,
        &metadata.root,
        &metadata.input,
    )
}

/// Compute the metadata file path under the build directory.
fn metadata_path(build_dir: &Path) -> PathBuf {
    build_dir.join(".cache").join("run-metadata.json")
}

/// Persist the latest run metadata for future incremental rebuilds.
fn save_run_metadata(path: &Path, metadata: &RunMetadata) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create metadata directory {}", parent.display()))?;
    }
    let serialized = serde_json::to_vec_pretty(metadata).context("failed to serialize metadata")?;
    fs::write(path, serialized)
        .with_context(|| format!("failed to write metadata {}", path.display()))?;
    Ok(())
}

/// Read run metadata, guiding incremental rebuilds.
fn load_run_metadata(path: &Path) -> Result<RunMetadata> {
    let data = match fs::read(path) {
        Ok(bytes) => bytes,
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {
            bail!(
                "run metadata missing at {} – run `linnet <root>` first",
                path.display()
            );
        }
        Err(err) => {
            return Err(err).with_context(|| format!("failed to read metadata {}", path.display()));
        }
    };
    let metadata = serde_json::from_slice(&data)
        .with_context(|| format!("failed to parse metadata {}", path.display()))?;
    Ok(metadata)
}

/// Canonicalize a path, returning an error if it does not exist.
fn canonicalize_existing(path: &Path) -> Result<PathBuf> {
    Ok(fs::canonicalize(path)?)
}

/// Join `path` onto `base` unless `path` is already absolute.
fn absolutize(base: &Path, path: &Path) -> PathBuf {
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        base.join(path)
    }
}

/// Derive a human-readable title from a relative DOT path.
fn derive_title(relative: &Path) -> String {
    let mut parts = Vec::new();
    for component in relative.components() {
        let part = component.as_os_str().to_string_lossy();
        parts.push(part);
    }
    let title = parts.join(" / ");
    if let Some(stripped) = title.strip_suffix(".dot") {
        stripped.to_owned()
    } else {
        title
    }
}

/// Hash the DOT data, figure template, and all style/layout files.
fn compute_hash(
    data: &Path,
    template: &Path,
    styles: &[PathBuf],
    inputs: &[(String, String)],
) -> Result<String> {
    let mut hasher = Hasher::new();
    feed_file(&mut hasher, data)?;
    feed_file(&mut hasher, template)?;
    for style in styles {
        feed_file(&mut hasher, style)?;
    }
    // Include input variables in the hash to invalidate cache when they change
    // Sort inputs for consistent hashing regardless of command line order
    let mut sorted_inputs = inputs.to_vec();
    sorted_inputs.sort_by(|a, b| a.0.cmp(&b.0));
    for (key, value) in sorted_inputs {
        hasher.update(key.as_bytes());
        hasher.update(b"=");
        hasher.update(value.as_bytes());
        hasher.update(b"\n");
    }
    Ok(hasher.finalize().to_hex().to_string())
}

/// Feed an entire file into the running hash.
fn feed_file(hasher: &mut Hasher, path: &Path) -> Result<()> {
    let mut file = File::open(path)
        .with_context(|| format!("failed to open {} for hashing", path.display()))?;
    let mut buffer = [0u8; 8192];
    loop {
        let read = file
            .read(&mut buffer)
            .with_context(|| format!("failed to read {} while hashing", path.display()))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(())
}

fn build_figure(
    plan: &FigurePlan,
    template: &Path,
    root: &Path,
    input_args: &[(String, String)],
) -> Result<()> {
    let template_dir = template
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."));
    let data_rel =
        diff_paths(&plan.data_path, &template_dir).unwrap_or_else(|| plan.data_path.clone());
    let relative_input = data_rel.to_string_lossy().replace('\\', "/");
    let mut command = Command::new("typst");
    let title = derive_title(&plan.relative);
    command
        .arg("c")
        .arg(template)
        .arg(&plan.output_path)
        .arg("--root")
        .arg(root)
        .arg("--input")
        .arg(format!("data_path={}", relative_input))
        .arg("--input")
        .arg(format!("title={}", title));

    // Add additional input arguments
    for (key, value) in input_args {
        command.arg("--input").arg(format!("{}={}", key, value));
    }

    run_typst(
        &mut command,
        &format!("building {}", plan.relative.display()),
        Some((plan, template, root)),
    )
}

fn run_typst(
    command: &mut Command,
    description: &str,
    context: Option<(&FigurePlan, &Path, &Path)>,
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
        error_msg.push_str(&format!("\n\nCommand executed:\n{}", cmd_str));

        // Add context-specific diagnostics
        if let Some((plan, template, root)) = context {
            error_msg.push_str("\n\nDiagnostic information:");
            error_msg.push_str(&format!("\n  - Input file: {}", plan.data_path.display()));
            error_msg.push_str(&format!("\n  - Template: {}", template.display()));
            error_msg.push_str(&format!("\n  - Root directory: {}", root.display()));
            error_msg.push_str(&format!(
                "\n  - Output path: {}",
                plan.output_path.display()
            ));

            // Check if files exist
            if !plan.data_path.exists() {
                error_msg.push_str("\n  - ERROR: Input file does not exist!");
            }
            if !template.exists() {
                error_msg.push_str("\n  - ERROR: Template file does not exist!");
            }
            if !root.exists() {
                error_msg.push_str("\n  - ERROR: Root directory does not exist!");
            }

            // Check if root contains the template
            if let Ok(template_canonical) = template.canonicalize()
                && let Ok(root_canonical) = root.canonicalize()
                && !template_canonical.starts_with(&root_canonical)
            {
                error_msg.push_str("\n  - WARNING: Template is outside root directory");
                error_msg.push_str(
                    "\n    This may cause 'source file must be contained in project root' errors",
                );
                error_msg.push_str("\n    Consider using --root with a parent directory of both template and data files");
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
            error_msg
                .push_str("\n\nFile not found - check that all paths are correct and files exist");
        }

        // Add the actual typst output
        error_msg.push_str(&format!("\n\nTypst stdout:\n{}", stdout.trim()));
        error_msg.push_str(&format!("\n\nTypst stderr:\n{}", stderr.trim()));

        bail!("{}", error_msg);
    }
    Ok(())
}

fn escape_typst_string(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
}

fn path_key(path: &Path) -> String {
    path.to_string_lossy().replace('\\', "/")
}

/// Load the per-figure hash cache from disk.
fn load_cache(path: &Path) -> Result<HashMap<String, String>> {
    if !path.exists() {
        return Ok(HashMap::new());
    }
    let data = fs::read_to_string(path)
        .with_context(|| format!("failed to read cache file {}", path.display()))?;
    if data.trim().is_empty() {
        return Ok(HashMap::new());
    }

    let parsed: JsonValue = serde_json::from_str(&data)
        .with_context(|| format!("failed to parse cache file {}", path.display()))?;
    let mut map = HashMap::new();
    if let Some(figs) = parsed.get("figures").and_then(JsonValue::as_object) {
        for (key, value) in figs {
            if let Some(hash) = value.as_str() {
                map.insert(key.clone(), hash.to_owned());
            }
        }
    }
    Ok(map)
}

/// Persist the updated per-figure hash cache.
fn save_cache(path: &Path, cache: &BTreeMap<String, String>) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create cache parent {}", parent.display()))?;
    }
    let mut figures = JsonMap::new();
    for (key, value) in cache {
        figures.insert(key.clone(), JsonValue::String(value.clone()));
    }
    let mut root = JsonMap::new();
    root.insert("figures".into(), JsonValue::Object(figures));
    let json = JsonValue::Object(root);
    let mut file = File::create(path)
        .with_context(|| format!("failed to open cache file {} for writing", path.display()))?;
    let serialized = serde_json::to_string_pretty(&json)?;
    file.write_all(serialized.as_bytes())
        .with_context(|| format!("failed to write cache file {}", path.display()))?;
    Ok(())
}

/// Delete figure PDFs that have no matching source entry anymore.
fn remove_stale_outputs(
    previous: &HashMap<String, String>,
    current: &BTreeMap<String, String>,
    figs_dir: &Path,
) -> Result<Vec<PathBuf>> {
    let mut removed = Vec::new();
    for key in previous.keys() {
        if current.contains_key(key) {
            continue;
        }
        let relative = Path::new(key);
        let mut pdf_path = figs_dir.join(relative);
        pdf_path.set_extension("pdf");
        match fs::remove_file(&pdf_path) {
            Ok(()) => removed.push(pdf_path),
            Err(err) if err.kind() == std::io::ErrorKind::NotFound => {}
            Err(err) => {
                return Err(err).with_context(|| {
                    format!("failed to delete stale figure {}", relative.display())
                });
            }
        }
    }
    Ok(removed)
}

/// Generate `fig-index.typ`, describing the figure tree for the grid template.
fn write_fig_index(records: &[FigureRecord], index_path: &Path, grid_dir: &Path) -> Result<()> {
    ensure_parent_dir(index_path)?;
    let mut file = File::create(index_path)
        .with_context(|| format!("failed to create {}", index_path.display()))?;
    writeln!(file, "// Auto-generated by linnet. Do not edit manually.")?;

    let mut entries = Vec::new();
    for record in records {
        let rel_path =
            diff_paths(&record.output_path, grid_dir).unwrap_or_else(|| record.output_path.clone());
        let rel_display = record.relative.to_string_lossy().replace('\\', "/");
        let file_stem = record
            .relative
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_default();
        let folder_parts = folder_components(&record.relative);
        let entry = FigureEntry {
            path: rel_path.to_string_lossy().replace('\\', "/"),
            relative: rel_display,
            name: file_stem,
        };
        entries.push((folder_parts, entry));
    }

    let mut root = FolderNode::default();
    for (folders, entry) in entries {
        insert_entry(&mut root, &folders, entry);
    }

    write!(file, "#let tree = ")?;
    write_folder_node(&mut file, &root, 0, false)?;
    writeln!(file)?;
    Ok(())
}

/// Resolve a template path, writing the embedded default if the requested file does not exist.
fn resolve_template(requested: &Path, kind: TemplateKind, build_dir: &Path) -> Result<PathBuf> {
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

    let templates_dir = build_dir.join(TEMPLATE_SUBDIR);
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
    ensure_parent_dir(&target)?;
    let contents = kind.embedded_bytes()?;
    fs::write(&target, contents.as_ref())
        .with_context(|| format!("failed to write default template {}", target.display()))?;
    Ok(target)
}

/// Ensure the embedded linnest plugin is available under `build/templates`.
fn ensure_plugin_asset(build_dir: &Path) -> Result<PathBuf> {
    let target = build_dir
        .join(TEMPLATE_SUBDIR)
        .join(TemplateKind::Plugin.file_name());
    if target.exists() {
        return Ok(target);
    }
    ensure_parent_dir(&target)?;
    let contents = TemplateKind::Plugin.embedded_bytes()?;
    fs::write(&target, contents.as_ref())
        .with_context(|| format!("failed to write embedded plugin {}", target.display()))?;
    Ok(target)
}

/// Ensure the default layout template exists unless the user already created one.
fn ensure_layout_asset(build_dir: &Path) -> Result<PathBuf> {
    let target = build_dir
        .join(TEMPLATE_SUBDIR)
        .join(TemplateKind::Layout.file_name());
    if target.exists() {
        return Ok(target);
    }
    ensure_parent_dir(&target)?;
    let contents = TemplateKind::Layout.embedded_bytes()?;
    fs::write(&target, contents.as_ref()).with_context(|| {
        format!(
            "failed to write default layout template {}",
            target.display()
        )
    })?;
    Ok(target)
}

fn folder_components(path: &Path) -> Vec<String> {
    path.parent()
        .map(|parent| {
            parent
                .components()
                .map(|component| component.as_os_str().to_string_lossy().into_owned())
                .collect()
        })
        .unwrap_or_default()
}

fn typst_string_literal(value: &str) -> String {
    format!("\"{}\"", escape_typst_string(value))
}

fn format_typst_array(items: &[String]) -> String {
    if items.is_empty() {
        "()".to_string()
    } else {
        let body = items
            .iter()
            .map(|item| typst_string_literal(item))
            .collect::<Vec<_>>()
            .join(", ");
        format!("({body},)")
    }
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

/// Compile the grid Typst template into the final PDF.
fn compile_grid(
    template: &Path,
    output: &Path,
    root: &Path,
    input_args: &[(String, String)],
) -> Result<()> {
    ensure_parent_dir(output)?;
    let mut cmd = Command::new("typst");
    cmd.arg("c")
        .arg(template)
        .arg(output)
        .arg("--root")
        .arg(root);

    // Add additional input arguments
    for (key, value) in input_args {
        cmd.arg("--input").arg(format!("{}={}", key, value));
    }

    run_typst(
        &mut cmd,
        &format!("compiling grid {}", output.display()),
        None,
    )
}

fn insert_entry(node: &mut FolderNode, folders: &[String], entry: FigureEntry) {
    if folders.is_empty() {
        node.figures.push(entry);
    } else {
        let (head, tail) = folders.split_first().unwrap();
        let child = node.children.entry(head.clone()).or_default();
        insert_entry(child, tail, entry);
    }
}

fn write_folder_node(
    file: &mut File,
    node: &FolderNode,
    indent: usize,
    trailing_comma: bool,
) -> Result<()> {
    let indent_str = indent_spaces(indent);
    writeln!(file, "{indent_str}(")?;
    let field_indent = indent + 2;
    let field_str = indent_spaces(field_indent);
    write_figures_field(file, &node.figures, field_indent)?;
    let child_names: Vec<String> = node.children.keys().cloned().collect();
    writeln!(
        file,
        "{field_str}order: {},",
        format_typst_array(&child_names)
    )?;
    if child_names.is_empty() {
        writeln!(file, "{field_str}folders: (:),")?;
    } else {
        writeln!(file, "{field_str}folders: (")?;
        let key_indent = indent_spaces(field_indent + 2);
        for (child_name, child_node) in &node.children {
            write!(file, "{key_indent}{}: ", typst_string_literal(child_name))?;
            write_folder_node(file, child_node, field_indent + 4, true)?;
        }
        writeln!(file, "{field_str}),")?;
    }
    writeln!(
        file,
        "{indent_str}){}",
        if trailing_comma { "," } else { "" }
    )?;
    Ok(())
}

fn write_figures_field(file: &mut File, figures: &[FigureEntry], indent: usize) -> Result<()> {
    let indent_str = indent_spaces(indent);
    if figures.is_empty() {
        writeln!(file, "{indent_str}figures: (),")?
    } else {
        writeln!(file, "{indent_str}figures: (")?;
        let entry_indent = indent + 2;
        let entry_str = indent_spaces(entry_indent);
        for figure in figures {
            writeln!(
                file,
                "{entry_str}(path: {}, relative: {}, name: {},),",
                typst_string_literal(&figure.path),
                typst_string_literal(&figure.relative),
                typst_string_literal(&figure.name)
            )?;
        }
        writeln!(file, "{indent_str}),")?;
    }
    Ok(())
}

fn indent_spaces(indent: usize) -> String {
    " ".repeat(indent)
}

/// Parses key/value pairs split by the first equal sign.
///
/// This function will return an error if the argument contains no equals sign
/// or the key (before the equals sign) is empty.
fn parse_input_pair(raw: &str) -> Result<(String, String), String> {
    let (key, val) = raw
        .split_once('=')
        .ok_or("input must be a key and a value separated by an equal sign")?;
    let key = key.trim().to_owned();
    if key.is_empty() {
        return Err("the key was missing or empty".to_owned());
    }
    let val = val.trim().to_owned();
    Ok((key, val))
}
