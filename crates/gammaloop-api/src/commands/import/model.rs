use clap::Args;
use color_eyre::Result;
use colored::Colorize;
use eyre::{eyre, Context};
use gammalooprs::utils::F;
#[cfg(feature = "ufo_support")]
use pyo3::prelude::*;

use std::path::{Path, PathBuf};
#[cfg(any(feature = "cli", test))]
use std::{collections::BTreeMap, sync::OnceLock};

use gammalooprs::model::{InputParamCard, Model, ModelGammaLoopExt};
use include_dir::{include_dir, Dir, File};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::{env, fs};
use tracing::info;

use crate::state::State;

static BUILTIN_MODELS: Dir = include_dir!("$CARGO_MANIFEST_DIR/../../assets/models/json");
#[cfg(feature = "cli")]
static BUILTIN_MODEL_NAMES: OnceLock<Vec<String>> = OnceLock::new();
#[cfg(any(feature = "cli", test))]
static BUILTIN_MODEL_RESTRICTIONS: OnceLock<BTreeMap<String, Vec<String>>> = OnceLock::new();

#[derive(Debug, Args, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
/// Imports a built-in JSON model, a JSON model file, or a UFO model directory into the active state.
pub struct ImportModel {
    /// Built-in model name, JSON model path, or UFO directory to import.
    #[arg(value_hint = clap::ValueHint::AnyPath)]
    pub path: PathBuf,
    /// Simplify imported model expressions before storing the active model.
    #[arg(short = 's', long,
          action = clap::ArgAction::Set,     // take an optional bool value
          default_value_t = true,
          require_equals = true)] // --simplify-model=false
    pub simplify_model: bool,
}

#[cfg(not(feature = "ufo_support"))]
pub(crate) fn load_ufo_model(
    _path: &Path,
    _restriction_name: Option<String>,
    _simplify_model: bool,
) -> Result<(Model, InputParamCard<F<f64>>)> {
    Err(eyre!(
        "UFO support not enabled in gammaloop. Please recompile with the `ufo_support` feature."
    ))
}

/// Releases before 1.0.0 silently discard particle chemical potentials.
#[cfg(feature = "ufo_support")]
fn require_chemical_potential_loader(version: &str, python: &str) -> eyre::Result<()> {
    if version
        .split('.')
        .next()
        .and_then(|major| major.parse::<u64>().ok())
        .is_none_or(|major| major < 1)
    {
        return Err(eyre!(
            "UFO import requires ufo-model-loader>=1.0.0; found {version}. Upgrade the loader in Python environment {python}."
        ));
    }
    Ok(())
}

#[cfg(feature = "ufo_support")]
pub(crate) fn load_ufo_model(
    path: &std::path::Path,
    restriction_name: Option<String>,
    simplify_model: bool,
) -> eyre::Result<(Model, InputParamCard<F<f64>>)> {
    use feynkit_ufo::UfoLoader;

    // Interpreter initialization belongs to GammaLoop's application boundary;
    // the reusable UFO adapter only attaches to an interpreter supplied by its caller.
    Python::initialize();

    Python::attach(|py| {
        let sys = py.import("sys")?;
        let exe: String = sys.getattr("prefix")?.extract()?;

        let metadata = py.import("importlib.metadata")?;
        let loader_version: String = match metadata.call_method1("version", ("ufo-model-loader",)) {
            Ok(version) => version.extract()?,
            Err(e) if e.is_instance(py, &metadata.getattr("PackageNotFoundError")?) => {
                return Err(eyre!(
                    "UFO import requires ufo-model-loader>=1.0.0, but Python environment {exe} has no installed ufo-model-loader metadata. A source checkout on PYTHONPATH cannot be version-checked; install it, for example with `pip install -e <checkout>`."
                ));
            }
            Err(e) => {
                return Err(eyre!(
                    "Could not read the ufo-model-loader version: {e}\nPython: {exe}"
                ));
            }
        };
        require_chemical_potential_loader(&loader_version, &exe)?;

        // Resolution has already selected any default card; an explicit full
        // import must not let the Python loader select its default again.
        let loader = UfoLoader::new()
            .simplify_model(simplify_model)
            .restriction_name(restriction_name.as_deref().unwrap_or("full"));
        let loaded = loader
            .load(py, path)
            .wrap_err_with(|| format!("Failed to load UFO model from {}", path.display()))?;

        Ok((loaded.model, loaded.parameters.into()))
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
                        info!(
                            "Loading {} model {} from {}",
                            "JSON".blue(),
                            format!("{model_name}-{restriction_name}").green(),
                            json_path.yellow()
                        );
                    } else {
                        info!(
                            "Loading {} model {} from {}",
                            "JSON".blue(),
                            format!("{model_name}-full").green(),
                            json_path.yellow()
                        );
                    }
                } else if let Some(restriction_name) = &restriction_name {
                    info!(
                        "Loading {} {} model {}",
                        "built-in".yellow(),
                        "JSON".blue(),
                        format!("{model_name}-{restriction_name}").green(),
                    );
                } else {
                    info!(
                        "Loading {} {} model {}",
                        "built-in".yellow(),
                        "JSON".blue(),
                        format!("{model_name}-full").green(),
                    );
                }
                state.model = Model::from_json(&json_model)?;
                state.model = state.model.with_restriction(restriction_name.clone())?;
                state.model_parameters = if let Some(json_restriction) = json_restriction {
                    let mut param_card = InputParamCard::from_str(json_restriction, "json")?;
                    if self.simplify_model {
                        state.model.simplify(&mut param_card)?;
                    } else {
                        state.model.apply_param_card(&param_card)?;
                    }
                    param_card
                } else {
                    state.model.default_param_card()
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
                    info!(
                        "Loading {} model {} from {}",
                        "UFO".blue(),
                        format!("{model_name}-{restriction_name}").green(),
                        ufo_path_string.yellow()
                    );
                } else {
                    info!(
                        "Loading {} model {} from {}",
                        "UFO".blue(),
                        format!("{model_name}-full").green(),
                        ufo_path_string.yellow()
                    );
                }
                (state.model, state.model_parameters) =
                    load_ufo_model(&ufo_path, restriction_name.clone(), self.simplify_model)?;
                state.model = state.model.with_restriction(restriction_name.clone())?;
            }
        }

        Ok(())
    }
}

#[cfg(feature = "cli")]
pub(crate) fn builtin_json_model_names() -> &'static [String] {
    BUILTIN_MODEL_NAMES.get_or_init(scan_builtin_json_model_names)
}

#[cfg(any(feature = "cli", test))]
pub(crate) fn builtin_json_model_restriction_names(model: &str) -> Option<&'static [String]> {
    BUILTIN_MODEL_RESTRICTIONS
        .get_or_init(scan_builtin_json_model_restriction_names)
        .get(model)
        .map(Vec::as_slice)
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

#[cfg(feature = "cli")]
fn scan_builtin_json_model_names() -> Vec<String> {
    let mut names = BUILTIN_MODELS
        .dirs()
        .filter_map(|dir| dir.path().file_name())
        .map(|name| name.to_string_lossy().to_string())
        .collect::<Vec<_>>();
    names.sort();
    names.dedup();
    names
}

#[cfg(any(feature = "cli", test))]
fn scan_builtin_json_model_restriction_names() -> BTreeMap<String, Vec<String>> {
    BUILTIN_MODELS
        .dirs()
        .filter_map(|dir| {
            let model_name = dir.path().file_name()?.to_string_lossy().to_string();
            let mut restrictions = dir
                .files()
                .filter_map(|file| file.path().file_name())
                .filter_map(|name| name.to_str())
                .filter_map(|name| {
                    name.strip_prefix("restrict_")
                        .and_then(|name| name.strip_suffix(".json"))
                        .map(str::to_string)
                })
                .collect::<Vec<_>>();
            restrictions.sort();
            restrictions.dedup();
            Some((model_name, restrictions))
        })
        .collect()
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

#[cfg(test)]
mod tests {
    use super::*;
    use gammalooprs::model::{ParameterNature, UFOSymbol};
    use spenso::algebra::complex::Complex;
    #[cfg(feature = "ufo_support")]
    use symbolica::atom::AtomCore;
    use symbolica::domains::float::Complex as NativeComplex;
    #[cfg(feature = "ufo_support")]
    use symbolica::domains::rational::Rational;
    use symbolica::{atom::Atom, parse};

    #[test]
    fn json_restrictions_canonicalize_only_fixed_zero_masses() -> Result<()> {
        gammalooprs::initialisation::test_initialise()?;
        let full = Model::from_json(&file_to_string(&ensure_builtin_model("sm")?)?)?;
        let restrictions = builtin_json_model_restriction_names("sm").unwrap();
        assert_eq!(restrictions.len(), 12);
        // These explicit massless Dirac propagators remain an independent oracle
        // when every quark in the unrestricted model has a finite default mass.
        let zero_numerator =
            Atom::i() * parse!("UFO::{}::PSlash(UFO::{}::idx(1,2),UFO::{}::idx(1,1))");
        let zero_antiparticle_numerator =
            -Atom::i() * parse!("UFO::{}::PSlash(UFO::{}::idx(1,1),UFO::{}::idx(1,2))");
        let zero_denominator = parse!("UFO::{}::P(UFO::{}::idx(1,1))^2");
        for restriction in std::iter::once("full").chain(restrictions.iter().map(String::as_str)) {
            for simplify_model in [false, true] {
                let mut state = State::new_test();
                ImportModel {
                    path: format!("sm-{restriction}").into(),
                    simplify_model,
                }
                .run(&mut state)?;
                for (quark, mass, yukawa_input, yukawa, default) in [
                    ("u", "MU", "ymu", "yu", 0.0022),
                    ("d", "MD", "ymd", "yd", 0.0047),
                    ("s", "MS", "yms", "ys", 0.095),
                    ("c", "MC", "ymc", "yc", 1.27),
                ] {
                    let expected_value = match (restriction, quark) {
                        ("full", _) | ("thermal_s_mass", "s") | ("light_quarks_massless", "c") => {
                            default
                        }
                        ("c_mass", "c") => 1.55,
                        _ => 0.0,
                    };
                    let massless = simplify_model && expected_value == 0.0;
                    let expected_mass = if massless {
                        UFOSymbol::zero()
                    } else {
                        UFOSymbol::from(mass)
                    };
                    for parameter in [mass, yukawa_input] {
                        assert_eq!(
                            state.model.get_parameter(parameter).nature,
                            if massless {
                                ParameterNature::Internal
                            } else {
                                ParameterNature::External
                            },
                            "{restriction}: {parameter}, simplify={simplify_model}"
                        );
                        assert_eq!(
                            state.model.get_parameter(parameter).value,
                            Some(NativeComplex::new(expected_value, 0.0)),
                            "{restriction}: {parameter}, simplify={simplify_model}"
                        );
                    }
                    for (name, zero_numerator) in [
                        (quark.to_string(), &zero_numerator),
                        (format!("{quark}~"), &zero_antiparticle_numerator),
                    ] {
                        let particle_id = state.model.get_particle(&name);
                        let particle = state.model.particle_by_id(particle_id)?;
                        assert_eq!(
                            UFOSymbol::from(
                                state.model.parameter_by_id(particle.mass)?.name.as_str()
                            ),
                            expected_mass,
                            "{restriction}: {name}"
                        );
                        assert_eq!(
                            state.model.unresolved_particles()["QCD"].contains(&particle_id),
                            massless,
                            "{restriction}: {name}, simplify={simplify_model}"
                        );
                        let propagator = state.model.propagator(&format!("{name}_propFeynman"))?;
                        let reference = full.propagator(&format!("{name}_propFeynman"))?;
                        assert_eq!(
                            UFOSymbol::from(
                                state
                                    .model
                                    .parameter_by_id(
                                        state.model.particle_by_id(propagator.particle)?.mass
                                    )?
                                    .name
                                    .as_str()
                            ),
                            expected_mass
                        );
                        assert_eq!(
                            &propagator.numerator,
                            if massless {
                                zero_numerator
                            } else {
                                &reference.numerator
                            }
                        );
                        assert_eq!(
                            &propagator.denominator,
                            if massless {
                                &zero_denominator
                            } else {
                                &reference.denominator
                            }
                        );
                    }
                    let antiquark = format!("{quark}~");
                    for scalar in ["H", "G0"] {
                        assert_eq!(
                            state.model.vertex_rules().iter().any(|vertex| {
                                vertex
                                    .particles
                                    .iter()
                                    .map(|particle| {
                                        state.model.particle_by_id(*particle).unwrap().name.as_str()
                                    })
                                    .eq([antiquark.as_str(), quark, scalar])
                            }),
                            !massless,
                            "{restriction}: {quark} {scalar}, simplify={simplify_model}"
                        );
                    }
                    if !massless {
                        let higgs = state
                            .model
                            .vertex_rules()
                            .iter()
                            .find(|vertex| {
                                vertex
                                    .particles
                                    .iter()
                                    .map(|particle| {
                                        state.model.particle_by_id(*particle).unwrap().name.as_str()
                                    })
                                    .eq([antiquark.as_str(), quark, "H"])
                            })
                            .unwrap();
                        let coupling = *higgs.couplings.iter().flatten().flatten().next().unwrap();
                        // Pole masses and Yukawa inputs are independent and stay
                        // adjustable through runtime zeros without dropping vertices.
                        for (pole_mass, yukawa_mass) in [(0.0, 0.27), (0.19, 0.0)] {
                            state
                                .model_parameters
                                .insert(mass.into(), Complex::new_re(F(pole_mass)));
                            state
                                .model_parameters
                                .insert(yukawa_input.into(), Complex::new_re(F(yukawa_mass)));
                            state.model.apply_param_card(&state.model_parameters)?;
                            let ve = state.model.get_parameter("vev").value.unwrap().re;
                            assert_eq!(
                                state.model.get_parameter(mass).value,
                                Some(NativeComplex::new(pole_mass, 0.0))
                            );
                            assert_eq!(
                                UFOSymbol::from(
                                    state
                                        .model
                                        .parameter_by_id(state.model.particle(quark)?.mass)?
                                        .name
                                        .as_str()
                                ),
                                UFOSymbol::from(mass)
                            );
                            let y = state.model.get_parameter(yukawa).value.unwrap();
                            assert!(
                                (y.re - std::f64::consts::SQRT_2 * yukawa_mass / ve).abs()
                                    < 1.0e-15
                            );
                            assert_eq!(y.im, 0.0);
                            let value = state.model.coupling_by_id(coupling)?.value.unwrap();
                            assert!(value.re.abs() < 1.0e-15);
                            assert!((value.im + yukawa_mass / ve).abs() < 1.0e-15);
                        }
                        for parameter in [mass, yukawa_input] {
                            state
                                .model_parameters
                                .insert(parameter.into(), Complex::new_re(F(expected_value)));
                        }
                        state.model.apply_param_card(&state.model_parameters)?;
                    }
                }
                if restriction == "light_quarks_massless" {
                    if simplify_model {
                        assert_eq!(state.model.vertex_rules().len(), 153);
                        assert_eq!(state.model.couplings().len(), 108);
                    }
                    // This fixture removes only light-quark masses and Yukawas;
                    // heavy flavors, leptons, and CKM mixing remain adjustable.
                    let retained = [
                        ("MC", 1.27),
                        ("MB", 4.7),
                        ("MT", 173.0),
                        ("Me", 0.000511),
                        ("MM", 0.10566),
                        ("MTA", 1.777),
                        ("ymc", 1.27),
                        ("ymb", 4.2),
                        ("ymt", 164.5),
                        ("yme", 0.000511),
                        ("ymm", 0.10566),
                        ("ymtau", 1.777),
                        ("lamWS", 0.2253),
                        ("AWS", 0.808),
                        ("rhoWS", 0.132),
                        ("etaWS", 0.341),
                    ];
                    for (name, default) in retained {
                        let parameter = state.model.get_parameter(name);
                        assert_eq!(parameter.nature, ParameterNature::External, "{name}");
                        assert_eq!(
                            parameter.value,
                            Some(NativeComplex::new(default, 0.0)),
                            "{name}"
                        );
                        state
                            .model_parameters
                            .insert(name.into(), Complex::new_re(F(default * 1.25)));
                    }
                    state.model.apply_param_card(&state.model_parameters)?;
                    for (name, default) in retained {
                        assert_eq!(
                            state.model.get_parameter(name).value,
                            Some(NativeComplex::new(default * 1.25, 0.0)),
                            "{name} remains adjustable"
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[cfg(feature = "ufo_support")]
    #[test]
    fn ufo_import_rejects_loader_without_chemical_potential_support() {
        for version in ["0.1.8", "dev"] {
            let error = require_chemical_potential_loader(version, "python")
                .unwrap_err()
                .to_string();
            assert!(error.contains("ufo-model-loader>=1.0.0"), "{error}");
            assert!(error.contains(version), "{error}");
        }
        for version in ["1.0.0", "12.3"] {
            require_chemical_potential_loader(version, "python").unwrap();
        }
    }

    #[cfg(feature = "ufo_support")]
    #[test]
    fn ufo_and_json_restrictions_preserve_model_content() -> Result<()> {
        gammalooprs::initialisation::test_initialise()?;
        let assets = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../assets/models");
        let restrictions = builtin_json_model_restriction_names("sm").unwrap();
        for restriction in std::iter::once("full").chain(restrictions.iter().map(String::as_str)) {
            let path = assets.join(format!("ufo/sm-{restriction}"));
            let ModelSpecification::UFOModelSpecification {
                restriction_name, ..
            } = ModelSpecification::parse(path.to_str().unwrap())?
            else {
                panic!("an explicit UFO path must use the Python loader");
            };
            let (mut ufo, mut ufo_card) =
                load_ufo_model(&assets.join("ufo/sm"), restriction_name, true)?;
            let mut state = State::new_test();
            ImportModel {
                path: format!("sm-{restriction}").into(),
                simplify_model: true,
            }
            .run(&mut state)?;
            let mut json = state.model;
            let mut json_card = state.model_parameters;
            let restored = Model::from_json(&ufo.to_json()?)?;
            for name in ["muB", "muQ", "muLe", "muLmu", "muLtau"] {
                assert_eq!(
                    ufo.get_parameter(name).nature,
                    json.get_parameter(name).nature
                );
                assert_eq!(
                    ufo.get_parameter(name).value,
                    json.get_parameter(name).value
                );
            }
            assert_eq!(ufo.particles().len(), json.particles().len());
            for fresh in ufo.particles() {
                let name = fresh.name.as_str();
                let saved = json.particle(name)?;
                assert_eq!(fresh.charge, saved.charge, "{name}");
                assert_eq!(fresh.y_charge, saved.y_charge, "{name}");
                assert_eq!(fresh.y_charge_right, saved.y_charge_right, "{name}");
                assert_eq!(
                    fresh
                        .chemical_potential
                        .map(|id| &ufo.parameter_by_id(id).unwrap().name),
                    saved
                        .chemical_potential
                        .map(|id| &json.parameter_by_id(id).unwrap().name),
                    "{name}"
                );
                assert_eq!(
                    ufo.parameter_by_id(fresh.mass)?.name,
                    json.parameter_by_id(saved.mass)?.name,
                    "{restriction}: {name} pole mass"
                );
                let roundtrip = restored.particle(name)?;
                assert_eq!(
                    ufo.parameter_by_id(fresh.mass)?.name,
                    restored.parameter_by_id(roundtrip.mass)?.name,
                    "{name} pole mass roundtrip"
                );
                assert_eq!(fresh.charge, roundtrip.charge, "{name} roundtrip");
                assert_eq!(fresh.y_charge, roundtrip.y_charge, "{name} roundtrip");
                assert_eq!(
                    fresh.y_charge_right, roundtrip.y_charge_right,
                    "{name} roundtrip"
                );
                assert_eq!(
                    fresh
                        .chemical_potential
                        .map(|id| &ufo.parameter_by_id(id).unwrap().name),
                    roundtrip
                        .chemical_potential
                        .map(|id| &restored.parameter_by_id(id).unwrap().name),
                    "{name} roundtrip"
                );
            }
            assert_eq!(ufo.particle("u")?.charge, Rational::from((2, 3)));
            assert_eq!(ufo.particle("H")?.y_charge, None);
            for model in [&ufo, &json, &restored] {
                for (particle, parameter) in [
                    ("e-", "mue"),
                    ("e+", "minus_mue"),
                    ("ve", "muve"),
                    ("ve~", "minus_muve"),
                ] {
                    assert_eq!(
                        UFOSymbol::from(
                            model
                                .parameter_by_id(
                                    model.particle(particle)?.chemical_potential.unwrap()
                                )?
                                .name
                                .as_str()
                        ),
                        if restriction == "full" {
                            UFOSymbol::from(parameter)
                        } else {
                            UFOSymbol::zero()
                        },
                        "{restriction}: {particle} chemical potential"
                    );
                }
            }
            for (quark, mass, yukawa_input, default) in [
                ("u", "MU", "ymu", 0.0022),
                ("d", "MD", "ymd", 0.0047),
                ("s", "MS", "yms", 0.095),
                ("c", "MC", "ymc", 1.27),
            ] {
                let expected_value = match (restriction, quark) {
                    ("full", _) | ("thermal_s_mass", "s") | ("light_quarks_massless", "c") => {
                        default
                    }
                    ("c_mass", "c") => 1.55,
                    _ => 0.0,
                };
                let massless = expected_value == 0.0;
                for model in [&ufo, &json, &restored] {
                    for parameter in [mass, yukawa_input] {
                        assert_eq!(
                            model.get_parameter(parameter).nature,
                            if massless {
                                ParameterNature::Internal
                            } else {
                                ParameterNature::External
                            },
                            "{restriction}: {parameter}"
                        );
                        assert_eq!(
                            model.get_parameter(parameter).value,
                            Some(NativeComplex::new(expected_value, 0.0)),
                            "{restriction}: {parameter}"
                        );
                    }
                    for name in [quark.to_string(), format!("{quark}~")] {
                        let particle_id = model.get_particle(&name);
                        let particle = model.particle_by_id(particle_id)?;
                        assert_eq!(
                            UFOSymbol::from(model.parameter_by_id(particle.mass)?.name.as_str()),
                            if massless {
                                UFOSymbol::zero()
                            } else {
                                UFOSymbol::from(mass)
                            },
                            "{restriction}: {name}"
                        );
                        assert_eq!(
                            model.unresolved_particles()["QCD"].contains(&particle_id),
                            massless,
                            "{restriction}: {name} unresolved"
                        );
                        let propagator = model.propagator(&format!("{name}_propFeynman"))?;
                        let reference = json.propagator(&format!("{name}_propFeynman"))?;
                        assert_eq!(
                            model.particle_by_id(propagator.particle)?.mass,
                            particle.mass
                        );
                        assert_eq!(propagator.numerator, reference.numerator);
                        assert_eq!(propagator.denominator, reference.denominator);
                    }
                }
            }
            // Compare surviving interaction tensors, including the coupling slots:
            // merely matching vertex names misses a lost charged-Goldstone chirality.
            for model in [&ufo, &restored] {
                assert_eq!(
                    model.vertex_rules().len(),
                    json.vertex_rules().len(),
                    "{restriction}"
                );
                for vertex in model.vertex_rules() {
                    let reference = json.vertex_rule(&vertex.name)?;
                    assert!(vertex
                        .particles
                        .iter()
                        .map(|p| model.particle_by_id(*p).unwrap().name.as_str())
                        .eq(reference.particles.iter().map(|p| json
                            .particle_by_id(*p)
                            .unwrap()
                            .name
                            .as_str())));
                    let tensors = [(model, vertex), (&json, reference)].map(
                        |(model, vertex)| -> Result<Atom> {
                            let mut tensor = Atom::Zero;
                            for (color, row) in vertex.couplings.iter().enumerate() {
                                for (lorentz, coupling) in row.iter().enumerate() {
                                    if let Some(coupling) = coupling {
                                        tensor += &model.coupling_by_id(*coupling)?.expression
                                            * &vertex.color_structures[color]
                                            * &model
                                                .lorentz_structure_by_id(
                                                    vertex.lorentz_structures[lorentz],
                                                )?
                                                .structure;
                                    }
                                }
                            }
                            Ok(tensor.expand())
                        },
                    );
                    let [left, right] = tensors;
                    assert_eq!(left?, right?, "{restriction}: {}", vertex.name);
                }
            }
            if restriction != "full" && !restriction.starts_with("thermal") {
                assert_eq!(json.get_parameter("muB").nature, ParameterNature::Internal);
                assert_eq!(
                    json.get_parameter("muB").value,
                    Some(NativeComplex::new(0.0, 0.0))
                );
            } else {
                // A zero-density comparison must not freeze later finite-density runs.
                for mu in [0.0, 6.0] {
                    for (model, card) in [(&mut ufo, &mut ufo_card), (&mut json, &mut json_card)] {
                        assert_eq!(model.get_parameter("muB").nature, ParameterNature::External);
                        card.insert("muB".into(), Complex::new_re(F(mu)));
                        model.apply_param_card(card)?;
                        assert_eq!(
                            model.get_parameter("muu").value,
                            Some(NativeComplex::new(mu / 3.0, 0.0))
                        );
                        assert_eq!(
                            model.get_parameter("minus_muu").value,
                            Some(NativeComplex::new(-mu / 3.0, 0.0))
                        );
                    }
                }
            }
            if matches!(restriction, "full" | "thermal_s_mass") {
                for (model, card) in [(&mut ufo, &mut ufo_card), (&mut json, &mut json_card)] {
                    let higgs = model
                        .vertex_rules()
                        .iter()
                        .find(|vertex| {
                            vertex
                                .particles
                                .iter()
                                .map(|particle| {
                                    model.particle_by_id(*particle).unwrap().name.as_str()
                                })
                                .eq(["s~", "s", "H"])
                        })
                        .unwrap();
                    let coupling = *higgs.couplings.iter().flatten().flatten().next().unwrap();
                    for (mass, yukawa_mass) in [(0.0, 0.27), (0.19, 0.0)] {
                        card.insert("MS".into(), Complex::new_re(F(mass)));
                        card.insert("yms".into(), Complex::new_re(F(yukawa_mass)));
                        model.apply_param_card(card)?;
                        assert_eq!(
                            UFOSymbol::from(
                                model
                                    .parameter_by_id(model.particle("s")?.mass)?
                                    .name
                                    .as_str()
                            ),
                            UFOSymbol::from("MS")
                        );
                        assert_eq!(
                            model.get_parameter("MS").value,
                            Some(NativeComplex::new(mass, 0.0))
                        );
                        let vev = model.get_parameter("vev").value.unwrap().re;
                        let value = model.coupling_by_id(coupling)?.value.unwrap();
                        assert!(value.re.abs() < 1.0e-15);
                        assert!((value.im + yukawa_mass / vev).abs() < 1.0e-15);
                    }
                }
            }
            if restriction == "full" {
                assert_eq!(json.get_parameter("Me").nature, ParameterNature::External);
                assert_eq!(ufo.get_parameter("Me").nature, ParameterNature::External);
            }
        }
        Ok(())
    }

    #[cfg(feature = "ufo_support")]
    #[test]
    fn fresh_sm_ufo_preserves_virtual_gauge_and_covariant_cut_states() -> Result<()> {
        gammalooprs::initialisation::test_initialise()?;
        let assets = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../assets/models");
        let (imported, _) = load_ufo_model(&assets.join("ufo/sm"), None, false)?;
        let canonical = Model::from_path(assets.join("json/sm/sm.json"))?;
        assert_eq!(
            imported.covariant_cut_multiplets(),
            canonical.covariant_cut_multiplets()
        );
        for name in [
            "W+", "W-", "Z", "G+", "G-", "G0", "ghWp", "ghWp~", "ghWm", "ghWm~", "ghZ", "ghZ~",
        ] {
            let fresh = imported
                .propagators()
                .iter()
                .find(|p| imported.particle_by_id(p.particle).unwrap().name.as_str() == name)
                .unwrap();
            let saved = canonical
                .propagators()
                .iter()
                .find(|p| canonical.particle_by_id(p.particle).unwrap().name.as_str() == name)
                .unwrap();
            assert_eq!(fresh.numerator, saved.numerator, "{name} numerator");
            assert_eq!(fresh.denominator, saved.denominator, "{name} denominator");
            assert_eq!(
                imported
                    .particle_by_id(fresh.particle)
                    .unwrap()
                    .is_goldstone(),
                canonical
                    .particle_by_id(saved.particle)
                    .unwrap()
                    .is_goldstone(),
                "{name} Goldstone metadata"
            );
            assert_eq!(
                imported
                    .parameter_by_id(imported.particle_by_id(fresh.particle)?.mass)?
                    .name,
                canonical
                    .parameter_by_id(canonical.particle_by_id(saved.particle)?.mass)?
                    .name,
                "{name} pole mass"
            );
        }
        Ok(())
    }
}
