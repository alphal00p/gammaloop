use std::path::{Path, PathBuf};

use clap::Args;
use color_eyre::Result;
use colored::Colorize;
use eyre::{eyre, Context};
use gammalooprs::utils::F;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use gammalooprs::{
    model::{InputParamCard, Model},
    status_info,
};
use include_dir::{include_dir, Dir, File};
use std::{env, fs};

#[cfg(feature = "ufo_support")]
use pyo3::prelude::*;
#[cfg(feature = "ufo_support")]
use pyo3::sync::GILOnceCell;

use crate::state::State;

static BUILTIN_MODELS: Dir = include_dir!("$CARGO_MANIFEST_DIR/../models/json");

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
/// Generate integrands
pub struct ImportModel {
    #[arg(value_hint = clap::ValueHint::AnyPath)]
    pub path: PathBuf,
    #[arg(short = 's', long,
          action = clap::ArgAction::Set,     // take an optional bool value
          default_value_t = true,
          require_equals = true)] // --simplify-model=false
    pub simplify_model: bool,
}

#[cfg(not(feature = "ufo_support"))]
pub(crate) fn load_ufo_model(
    path: &Path,
    restriction_name: Option<String>,
    simplify_model: bool,
) -> Result<(Model, InputParamCard<F<f64>>)> {
    Err(eyre!(
        "UFO support not enabled in gammaloop. Please recompile with the `ufo_support` feature."
    ))
}

#[cfg(feature = "ufo_support")]
#[pyclass]
struct PyLogStream;

#[cfg(feature = "ufo_support")]
#[pymethods]
impl PyLogStream {
    #[new]
    fn new() -> Self {
        Self
    }
    fn write(&self, s: &str) {
        for line in s.lines() {
            status_info!("{}", line);
        }
    }
    fn flush(&self) {}
}

#[cfg(feature = "ufo_support")]
static PY_LOG_BRIDGE_INIT: GILOnceCell<()> = GILOnceCell::new();

#[cfg(feature = "ufo_support")]
fn ensure_py_log_bridge(py: Python<'_>) -> PyResult<()> {
    PY_LOG_BRIDGE_INIT.get_or_try_init(py, || -> PyResult<()> {
        bridge_py_logging(py)?;
        Ok(())
    })?;
    Ok(())
}

#[cfg(feature = "ufo_support")]
fn bridge_py_logging(py: Python<'_>) -> PyResult<()> {
    let logging = py.import("logging")?;
    let stream = Py::new(py, PyLogStream {})?;

    // Handler that sends logging records to Rust
    let handler = logging
        .getattr("StreamHandler")?
        .call1((stream.clone_ref(py),))?;
    handler.call_method1("setLevel", ("DEBUG",))?;
    // handler.call_method1(
    //     "setFormatter",
    //     (logging
    //         .getattr("Formatter")?
    //         .call1(("%(levelname)s %(name)s: %(message)s",))?,),
    // )?;

    let root = logging.call_method0("getLogger")?;
    root.call_method1("addHandler", (handler,))?;
    root.call_method1("setLevel", ("DEBUG",))?;

    // Also capture plain prints
    let sys = py.import("sys")?;
    sys.setattr("stdout", stream.clone_ref(py))?;
    sys.setattr("stderr", stream)?;
    Ok(())
}

#[cfg(feature = "ufo_support")]
pub(crate) fn load_ufo_model(
    path: &std::path::Path,
    restriction_name: Option<String>,
    simplify_model: bool,
) -> eyre::Result<(Model, InputParamCard<F<f64>>)> {
    use pyo3::{prelude::*, types::PyDict};
    pyo3::prepare_freethreaded_python(); // safe if auto-init is also on

    Python::with_gil(|py| {
        // helpful diag if import fails
        let sys = py.import("sys")?;
        let exe: String = sys.getattr("prefix")?.extract()?;
        let ver: String = sys.getattr("version")?.extract()?;

        ensure_py_log_bridge(py)
            .map_err(|e| eyre::eyre!("Failed to bridge python logging to rust. Error: {}", e))?;

        let commands = py.import("ufo_model_loader.commands").map_err(|e| {
            eyre::eyre!("Failed to import ufo_model_loader: {e}\nPython: {exe}\nVersion: {ver}")
        })?;

        let load_model = commands.getattr("load_model")?;

        // kwargs: input_model_path, restriction_name, simplify_model
        let kwargs = PyDict::new(py);
        kwargs.set_item("input_model_path", path.to_string_lossy().as_ref())?;
        match restriction_name.as_deref() {
            Some(name) => kwargs.set_item("restriction_name", name)?,
            None => kwargs.set_item("restriction_name", py.None())?,
        }
        kwargs.set_item("simplify_model", simplify_model)?;
        kwargs.set_item("wrap_indices_in_lorentz_structures", true)?;

        // call loader.load(...)
        let out = load_model
            .call((), Some(&kwargs))
            .map_err(|e| eyre::eyre!("ufo_model_loader.load failed: {e}"))?;

        // expect a 2-tuple: (model, input_param_card)
        let (py_model, py_card): (Py<PyAny>, Py<PyAny>) = out.extract()?;

        // to_json() -> String
        let model_json: String = py_model.call_method0(py, "to_json")?.extract(py)?;
        let card_json: String = py_card.call_method0(py, "to_json")?.extract(py)?;

        // deserialize into your Rust types
        let model: Model = Model::from_str(model_json, "json")
            .map_err(|e| eyre::eyre!("Failed to deserialize JSON Model: {e}"))?;
        let card: InputParamCard<F<f64>> = InputParamCard::<F<f64>>::from_str(card_json, "json")
            .map_err(|e| eyre::eyre!("Failed to deserialize InputParamCard: {e}"))?;

        Ok((model, card))
    })
}

impl ImportModel {
    pub fn run(&self, state: &mut State) -> Result<()> {
        let model_specification = ModelSpecification::parse(
            self.path
                .to_str()
                .ok_or_else(|| eyre!("Invalid model path"))?,
        )?;
        match model_specification {
            ModelSpecification::JSONModelSpecification {
                model_name,
                restriction_name,
                base_path,
                json_model,
                json_restriction,
            } => {
                if let Some(base_path) = base_path {
                    let json_path = PathBuf::from(base_path).join(model_name.clone());
                    let json_path = json_path.to_string_lossy();
                    if let Some(restriction_name) = &restriction_name {
                        status_info!(
                            "Loading {} model {} from {}",
                            "JSON".blue(),
                            format!("{model_name}-{restriction_name}").green(),
                            json_path.yellow()
                        );
                    } else {
                        status_info!(
                            "Loading {} model {} from {}",
                            "JSON".blue(),
                            format!("{model_name}-full").green(),
                            json_path.yellow()
                        );
                    }
                } else if let Some(restriction_name) = &restriction_name {
                    status_info!(
                        "Loading {} {} model {}",
                        "built-in".yellow(),
                        "JSON".blue(),
                        format!("{model_name}-{restriction_name}").green(),
                    );
                } else {
                    status_info!(
                        "Loading {} {} model {}",
                        "built-in".yellow(),
                        "JSON".blue(),
                        format!("{model_name}-full").green(),
                    );
                }
                state.model = Model::from_str(json_model, "json")?;
                state.model_parameters = if let Some(json_restriction) = json_restriction {
                    let mut param_card = InputParamCard::from_str(json_restriction, "json")?;
                    if self.simplify_model {
                        state.model.simplify(&mut param_card)?;
                    } else {
                        state.model.apply_param_card(&param_card)?;
                    }
                    param_card
                } else {
                    InputParamCard::default_from_model(&state.model)
                };
            }
            ModelSpecification::UFOModelSpecification {
                base_path,
                model_name,
                restriction_name,
            } => {
                let ufo_path = base_path.join(&model_name);
                let ufo_path_string = ufo_path.to_string_lossy();
                if let Some(restriction_name) = &restriction_name {
                    status_info!(
                        "Loading {} model {} from {}",
                        "UFO".blue(),
                        format!("{model_name}-{restriction_name}").green(),
                        ufo_path_string.yellow()
                    );
                } else {
                    status_info!(
                        "Loading {} model {} from {}",
                        "UFO".blue(),
                        format!("{model_name}-full").green(),
                        ufo_path_string.yellow()
                    );
                }
                (state.model, state.model_parameters) =
                    load_ufo_model(&ufo_path, restriction_name, self.simplify_model)?;
            }
        }

        Ok(())
    }
}

#[derive(Debug, Clone)]
pub enum ModelSpecification {
    JSONModelSpecification {
        model_name: String,
        restriction_name: Option<String>,
        base_path: Option<String>, // abs dir of the JSON model if from disk, else None
        json_model: String,        // JSON content
        json_restriction: Option<String>, // restriction JSON content if any
    },
    UFOModelSpecification {
        base_path: PathBuf, // abs dir containing the UFO folder
        model_name: String,
        restriction_name: Option<String>,
    },
}

#[derive(Debug, Clone, Copy)]
enum Origin {
    BuiltinJson,
    DiskJson,
    UfoDisk,
}

impl ModelSpecification {
    pub fn parse(input: &str) -> Result<Self> {
        let is_json_path = input.ends_with(".json");
        let has_slash = input.contains('/') || input.contains('\\');
        let cwd_opt = env::current_dir().ok();

        let (origin, base_dir, model_name, restriction_candidate, disk_model_path_opt): (
            Origin,
            Option<PathBuf>,
            String,
            Option<String>,
            Option<PathBuf>,
        ) = if !has_slash && !input.starts_with("./") {
            if is_json_path {
                let (name, r_opt) = split_name_and_restriction(input, false)?;
                if builtin_has_model(&name) {
                    (Origin::BuiltinJson, None, name, r_opt, None)
                } else {
                    let cwd = cwd_opt.ok_or_else(|| eyre!("cannot fetch current working directory; supply an absolute path or use a builtin model name"))?;
                    let (model_file, base) = find_disk_json_model(&cwd, &name, r_opt.as_deref())?;
                    (Origin::DiskJson, Some(base), name, r_opt, Some(model_file))
                }
            } else {
                let (name, r_opt) = split_name_and_restriction(input, true)?;
                if builtin_has_model(&name) {
                    (Origin::BuiltinJson, None, name, r_opt, None)
                } else {
                    let cwd = cwd_opt.ok_or_else(|| {
                        eyre!("cannot fetch current working directory; supply an absolute path")
                    })?;
                    (Origin::UfoDisk, Some(cwd), name, r_opt, None)
                }
            }
        } else {
            let p = if input.starts_with("./") {
                let cwd = cwd_opt.ok_or_else(|| {
                    eyre!("cannot fetch current working directory; supply an absolute path")
                })?;
                cwd.join(input.trim_start_matches("./"))
            } else {
                PathBuf::from(input)
            };
            if !p.is_absolute() && !p.exists() {
                return Err(eyre!(
                    "cannot fetch current working directory; supply an absolute path"
                ));
            }
            decide_from_path(p)?
        };

        match origin {
            Origin::UfoDisk => {
                let base = base_dir.unwrap();
                let model_dir = base.join(&model_name);
                if !model_dir.is_dir() {
                    return Err(eyre!(
                        "UFO model folder not found at {}",
                        model_dir.display()
                    ));
                }
                let base_path = canonical_dir_allow_abs_missing(&base)?;
                let restriction_name = resolve_ufo_restriction(
                    &base_path,
                    &model_name,
                    restriction_candidate.as_deref(),
                )?;
                Ok(ModelSpecification::UFOModelSpecification {
                    base_path,
                    model_name,
                    restriction_name,
                })
            }
            Origin::BuiltinJson => {
                let model_file = ensure_builtin_model(&model_name)?;
                let json_model = file_to_string(&model_file)?;
                let restriction_name = resolve_json_restriction_builtin_name(
                    &model_name,
                    restriction_candidate.as_deref(),
                )?;
                let json_restriction = match &restriction_name {
                    Some(r) => {
                        let rf = ensure_builtin_restriction(&model_name, r)?;
                        Some(file_to_string(&rf)?)
                    }
                    None => None,
                };
                Ok(ModelSpecification::JSONModelSpecification {
                    model_name,
                    restriction_name,
                    base_path: None,
                    json_model,
                    json_restriction,
                })
            }
            Origin::DiskJson => {
                let base = base_dir.unwrap();
                let model_path = disk_model_path_opt.unwrap();
                let json_model = fs::read_to_string(&model_path)
                    .wrap_err_with(|| format!("failed to read {}", model_path.display()))?;
                let restriction_name = resolve_json_restriction_disk_name(
                    &base,
                    &model_name,
                    restriction_candidate.as_deref(),
                )?;
                let json_restriction = match &restriction_name {
                    Some(r) => {
                        let rp = base.join(format!("restrict_{r}.json"));
                        Some(
                            fs::read_to_string(&rp)
                                .wrap_err_with(|| format!("failed to read {}", rp.display()))?,
                        )
                    }
                    None => None,
                };
                let base_abs = canonical_dir_allow_abs_missing(&base)?;
                Ok(ModelSpecification::JSONModelSpecification {
                    model_name,
                    restriction_name,
                    base_path: Some(base_abs.display().to_string()),
                    json_model,
                    json_restriction,
                })
            }
        }
    }
}

// ---------- helpers ----------

#[allow(clippy::type_complexity)]
fn decide_from_path(
    p: PathBuf,
) -> Result<(
    Origin,
    Option<PathBuf>,
    String,
    Option<String>,
    Option<PathBuf>,
)> {
    let token = last_segment_str(&p)?;
    if token.ends_with(".json") {
        let (name, r_opt) = split_name_and_restriction(token, false)?;
        let base = p
            .parent()
            .ok_or_else(|| eyre!("invalid JSON path"))?
            .to_path_buf();
        let (model_file, base_dir) = find_disk_json_model(&base, &name, r_opt.as_deref())?;
        Ok((
            Origin::DiskJson,
            Some(base_dir),
            name,
            r_opt,
            Some(model_file),
        ))
    } else {
        let (name, r_opt) = split_name_and_restriction(token, true)?;
        let base = p
            .parent()
            .ok_or_else(|| eyre!("invalid UFO path"))?
            .to_path_buf();
        Ok((Origin::UfoDisk, Some(base), name, r_opt, None))
    }
}

fn find_disk_json_model(
    base: &Path,
    model_name: &str,
    r: Option<&str>,
) -> Result<(PathBuf, PathBuf)> {
    let primary = base.join(assemble_json_filename(model_name, r));
    if primary.is_file() {
        return Ok((primary, base.to_path_buf()));
    }
    let alt = base.join(format!("{model_name}.json"));
    if r.is_some() && alt.is_file() {
        return Ok((alt, base.to_path_buf()));
    }
    Err(eyre!(
        "JSON model file not found at {} or {}",
        primary.display(),
        alt.display()
    ))
}

fn resolve_ufo_restriction(
    base: &Path,
    model_name: &str,
    given: Option<&str>,
) -> Result<Option<String>> {
    // If user requested "full", treat as None.
    if let Some(r) = given {
        if r.eq_ignore_ascii_case("full") {
            return Ok(None);
        }
    }
    match given {
        Some(r) => {
            let dat = base.join(model_name).join(format!("restrict_{r}.dat"));
            if dat.exists() {
                Ok(Some(r.to_string()))
            } else {
                Err(eyre!("missing {}", dat.display()))
            }
        }
        None => {
            let r = "default";
            let dat = base.join(model_name).join(format!("restrict_{r}.dat"));
            if dat.exists() {
                Ok(Some(r.to_string()))
            } else {
                Ok(None)
            }
        }
    }
}

fn resolve_json_restriction_builtin_name(
    model: &str,
    given: Option<&str>,
) -> Result<Option<String>> {
    // If user requested "full", treat as None.
    if let Some(r) = given {
        if r.eq_ignore_ascii_case("full") {
            return Ok(None);
        }
    }
    let exists = |r: &str| {
        BUILTIN_MODELS
            .get_file(format!("{model}/restrict_{r}.json"))
            .is_some()
    };
    match given {
        Some(r) => {
            if exists(r) {
                Ok(Some(r.to_string()))
            } else {
                Err(eyre!("builtin restriction not found for '{}'", r))
            }
        }
        None => {
            let r = "default";
            if exists(r) {
                Ok(Some(r.to_string()))
            } else {
                Ok(None)
            }
        }
    }
}

fn resolve_json_restriction_disk_name(
    base: &Path,
    _model: &str,
    given: Option<&str>,
) -> Result<Option<String>> {
    // If user requested "full", treat as None.
    if let Some(r) = given {
        if r.eq_ignore_ascii_case("full") {
            return Ok(None);
        }
    }
    let exists = |r: &str| base.join(format!("restrict_{r}.json")).exists();
    match given {
        Some(r) => {
            if exists(r) {
                Ok(Some(r.to_string()))
            } else {
                Err(eyre!("restriction file not found for '{}'", r))
            }
        }
        None => {
            let r = "default";
            if exists(r) {
                Ok(Some(r.to_string()))
            } else {
                Ok(None)
            }
        }
    }
}

fn split_name_and_restriction(token: &str, is_ufo_token: bool) -> Result<(String, Option<String>)> {
    let core = if !is_ufo_token && token.ends_with(".json") {
        Path::new(token)
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| eyre!("invalid JSON filename"))?
            .to_string()
    } else {
        token.to_string()
    };
    let parts: Vec<&str> = core.split('-').collect();
    match parts.len() {
        1 => Ok((parts[0].to_string(), None)),
        2 => {
            if parts[0].is_empty() || parts[1].is_empty() {
                return Err(eyre!("invalid model/restriction token '{}'", core));
            }
            Ok((parts[0].to_string(), Some(parts[1].to_string())))
        }
        _ => Err(eyre!("model name cannot contain hyphens; found '{}'", core)),
    }
}

fn assemble_json_filename(model_name: &str, restriction: Option<&str>) -> String {
    match restriction {
        Some(r) => format!("{model_name}-{r}.json"),
        None => format!("{model_name}.json"),
    }
}

fn builtin_has_model(model: &str) -> bool {
    BUILTIN_MODELS
        .get_file(format!("{model}/{model}.json"))
        .is_some()
}

fn ensure_builtin_model(model: &str) -> Result<File<'static>> {
    BUILTIN_MODELS
        .get_file(format!("{model}/{model}.json"))
        .ok_or_else(|| eyre!("builtin model '{model}/{model}.json' not found"))
        .cloned()
}

fn ensure_builtin_restriction(model: &str, r: &str) -> Result<File<'static>> {
    BUILTIN_MODELS
        .get_file(format!("{model}/restrict_{r}.json"))
        .ok_or_else(|| eyre!("builtin restriction '{model}/restrict_{r}.json' not found"))
        .cloned()
}

fn file_to_string(f: &File<'_>) -> Result<String> {
    f.contents_utf8()
        .map(|s| s.to_string())
        .ok_or_else(|| eyre!("builtin file '{}' is not valid UTF-8", f.path().display()))
}

fn last_segment_str(p: &Path) -> Result<&str> {
    p.file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| eyre!("invalid path segment"))
}

fn canonical_dir_allow_abs_missing(p: &Path) -> Result<PathBuf> {
    if p.exists() {
        p.canonicalize()
            .wrap_err("failed to canonicalize base path")
    } else if p.is_absolute() {
        Ok(p.to_path_buf())
    } else if let Some(parent) = p.parent() {
        let c = parent
            .canonicalize()
            .wrap_err("failed to canonicalize parent of base path")?;
        Ok(c.join(p.file_name().unwrap()))
    } else {
        Err(eyre!("invalid base path '{}'", p.display()))
    }
}
