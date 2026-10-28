use dirs::home_dir;
use serde::{Deserialize, Deserializer, Serialize, de::DeserializeOwned};
use std::{
    env,
    io::Write,
    path::{Path, PathBuf},
    process::Command,
    sync::{
        Mutex, MutexGuard,
        atomic::{AtomicBool, Ordering},
    },
};

use color_eyre::{Result, Section};
use eyre::{Context, eyre};
use std::collections::BTreeMap;
use std::result::Result as StdResult;
use thiserror::Error;
use toml::ser::Error as TomlEncodeError;

use std::{fs, fs::File, io::Read};

#[derive(Error, Debug)]
pub enum SerdeFileError {
    #[error("File error: {0}")]
    FileError(#[from] std::io::Error),
    #[error("JSON parse error: {0}")]
    JsonParseError(#[from] serde_json::Error),
    #[error("YAML parse error: {0}")]
    YamlParseError(#[from] serde_yaml::Error),
    #[error("TOML parse error: {0}")]
    TomlParseError(#[from] toml::de::Error),
    #[error("Unknown file extension: {0}")]
    UnknownExtension(String),
    #[error("Could not determine file extension of file {0}")]
    NoExtension(String),
}

/// Preserve map entries as an ordered sequence of key-value pairs, including non-string keys.
pub(crate) mod map_as_pairs {
    use serde::{Deserialize, Deserializer};

    pub fn deserialize<'de, T, K, V, D>(deserializer: D) -> Result<T, D::Error>
    where
        D: Deserializer<'de>,
        T: FromIterator<(K, V)>,
        K: Deserialize<'de>,
        V: Deserialize<'de>,
    {
        let entries: Vec<_> = Deserialize::deserialize(deserializer)?;
        Ok(entries.into_iter().collect())
    }
}

pub(crate) fn deserialize_nonnegative_finite_f64<'de, D>(
    deserializer: D,
) -> std::result::Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let value = f64::deserialize(deserializer)?;
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(serde::de::Error::custom(
            "expected a finite, nonnegative floating-point value",
        ))
    }
}

fn current_git_commit(git_dir: &Path) -> Option<String> {
    let output = Command::new("git")
        .arg("rev-parse")
        .arg("HEAD")
        .current_dir(git_dir)
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }

    let commit = String::from_utf8(output.stdout).ok()?;
    let commit = commit.trim();
    if commit.is_empty() {
        None
    } else {
        Some(commit.to_owned())
    }
}

fn schema_ref() -> String {
    env::var("GAMMALOOP_SCHEMA_REF")
        .ok()
        .or_else(|| current_git_commit(Path::new(".")))
        .or_else(|| current_git_commit(Path::new(env!("CARGO_MANIFEST_DIR"))))
        .unwrap_or_else(|| format!("v{}", env!("CARGO_PKG_VERSION")))
}

/// Whether smart TOML serialization includes values equal to their defaults.
/// See [`SmartSerde::to_toml`] for an example of both modes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SerdeDefaults {
    /// Include default-valued fields so persisted settings record the full configuration.
    Include,
    /// Omit default-valued fields using their serialization predicates.
    Omit,
}

pub trait SmartSerde: Serialize + DeserializeOwned {
    /// Serialize TOML with an optional editor schema and explicit default visibility.
    /// Schema reference resolution precedes the defaults guard; only serialization
    /// holds the guard, and no output file is opened or written here.
    /// Callers must not already hold a [`ShowDefaultsGuard`].
    ///
    /// # Errors
    /// Returns the TOML encoder error if the value cannot be represented as TOML.
    ///
    /// ```
    /// # use gammalooprs::settings::RuntimeSettings;
    /// # use gammalooprs::utils::serde_utils::{SerdeDefaults, SmartSerde};
    /// # use toml::ser::Error as TomlEncodeError;
    /// let settings = RuntimeSettings::default();
    /// let persisted = settings.to_toml(SerdeDefaults::Include)?;
    /// let compact = settings.to_toml(SerdeDefaults::Omit)?;
    /// assert!(persisted.starts_with("#:schema "));
    /// assert!(persisted.len() > compact.len());
    /// # Ok::<(), TomlEncodeError>(())
    /// ```
    fn to_toml(&self, defaults: SerdeDefaults) -> StdResult<String, TomlEncodeError> {
        let mut text = self
            .schema_path()
            .map(|schema| format!("#:schema {}\n", schema.display()))
            .unwrap_or_default();
        let serialized = {
            let _guard = ShowDefaultsGuard::new(defaults == SerdeDefaults::Include);
            toml::to_string_pretty(self)?
        };
        text.push_str(&serialized);
        Ok(text)
    }

    /// Write JSON, YAML, or TOML according to the destination's extension.
    /// TOML always includes defaults and any editor schema reference; JSON and YAML
    /// follow the ambient defaults setting. Callers must not already hold a
    /// [`ShowDefaultsGuard`]: TOML serialization acquires this non-reentrant guard.
    /// Missing parent directories are created; existing files require `override_existing`.
    ///
    /// # Errors
    /// Returns filesystem, missing or unsupported extension, and serialization errors.
    /// File creation or truncation occurs before serialization, so failure may leave
    /// an incomplete file. Callers needing atomic publication must stage the output.
    ///
    /// ```no_run
    /// # use color_eyre::Report;
    /// # use gammalooprs::settings::RuntimeSettings;
    /// # use gammalooprs::utils::serde_utils::SmartSerde;
    /// RuntimeSettings::default().to_file("settings/runtime.toml", false)?;
    /// # Ok::<(), Report>(())
    /// ```
    fn to_file(&self, file_path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        if let Some(parent) = file_path.as_ref().parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent)?;
        }

        let mut f = if override_existing {
            File::create(file_path.as_ref())?
        } else {
            File::create_new(file_path.as_ref())?
        };

        if let Some(ext) = file_path.as_ref().extension() {
            if let Some(ext) = ext.to_str() {
                match ext {
                    "json" => serde_json::to_writer_pretty(f, self)
                        .map_err(|e| eyre!(format!("Error serializing json: {}", e)))?,
                    "yaml" | "yml" => serde_yaml::to_writer(f, self)
                        .map_err(|e| eyre!(format!("Error serializing yaml: {}", e)))?,
                    "toml" => {
                        let toml_string = self.to_toml(SerdeDefaults::Include)?;
                        f.write_all(toml_string.as_bytes())?;
                    }
                    _ => return Err(eyre!(format!("Unknown file extension: {}", ext))),
                }
            } else {
                return Err(eyre!(format!(
                    "Could not determine file extension of file {}",
                    file_path.as_ref().display()
                )));
            }
        } else {
            return Err(eyre!(format!(
                "Could not determine file extension of file {}",
                file_path.as_ref().display()
            )));
        }

        Ok(())
    }

    /// Loads and deserializes data from a file with typed error handling.
    ///
    /// This function differentiates between file errors (e.g., file not found)
    /// and parse errors (e.g., invalid JSON/YAML/TOML syntax) using `SerdeFileError`.
    ///
    /// # Arguments
    /// * `file_path` - Path to the file to load
    /// * `_name` - Name of the data type being loaded (for error messages, currently unused)
    ///
    /// # Returns
    /// * `Ok(Self)` - Successfully loaded and parsed data
    /// * `Err(SerdeFileError)` - Typed error indicating whether it was a file or parse error
    fn from_file_typed(file_path: impl AsRef<Path>) -> std::result::Result<Self, SerdeFileError> {
        let mut f = File::open(file_path.as_ref())?;

        if let Some(ext) = file_path.as_ref().extension() {
            if let Some(ext) = ext.to_str() {
                match ext {
                    "json" => Ok(serde_json::from_reader(f)?),
                    "yaml" | "yml" => Ok(serde_yaml::from_reader(f)?),
                    "toml" => {
                        let mut buf = String::new();
                        let _bytes = f.read_to_string(&mut buf)?;
                        Ok(toml::from_str(&buf)?)
                    }
                    _ => Err(SerdeFileError::UnknownExtension(ext.to_string())),
                }
            } else {
                Err(SerdeFileError::NoExtension(
                    file_path.as_ref().display().to_string(),
                ))
            }
        } else {
            Err(SerdeFileError::NoExtension(
                file_path.as_ref().display().to_string(),
            ))
        }
    }

    fn from_file(file_path: impl AsRef<Path>, name: &str) -> Result<Self> {
        Self::from_file_typed(file_path.as_ref()).map_err(|e| match e {
            SerdeFileError::FileError(io_err) => eyre::Report::from(io_err)
                .wrap_err(format!(
                    "Could not open {} file {}",
                    name,
                    file_path.as_ref().display()
                ))
                .suggestion("Does the path exist?"),
            SerdeFileError::JsonParseError(json_err) => eyre::Report::from(json_err)
                .wrap_err(format!("Error parsing {} json", name))
                .suggestion("Is it a correct json file?"),
            SerdeFileError::YamlParseError(yaml_err) => eyre::Report::from(yaml_err)
                .wrap_err(format!("Error parsing {} yaml", name))
                .suggestion("Is it a correct yaml file?"),
            SerdeFileError::TomlParseError(toml_err) => eyre::Report::from(toml_err)
                .wrap_err(format!("Error parsing {} toml", name))
                .suggestion("Is it a correct toml file?"),
            SerdeFileError::UnknownExtension(ext) => {
                eyre::eyre!("Unknown {} file extension: {}", name, ext)
                    .suggestion("Is it a .json, .yaml, or .toml file?")
            }
            SerdeFileError::NoExtension(path) => eyre::eyre!(
                "Could not determine file extension of {} file {}",
                name,
                path
            )
            .suggestion("Does the path exist?"),
        })
    }

    fn from_str(contents: String, format: &str, name: &str) -> Result<Self> {
        match format {
            "json" => serde_json::from_str(&contents)
                .map_err(|e| eyre!(format!("Error parsing {name} json: {}", e)))
                .suggestion("Is it a correct json file"),
            "yaml" | "yml" => serde_yaml::from_str(&contents)
                .map_err(|e| eyre!(format!("Error parsing {name} yaml: {}", e)))
                .suggestion("Is it a correct yaml file"),
            "toml" => toml::from_str(&contents)
                .map_err(|e| eyre!(format!("Error parsing {name} toml: {}", e)))
                .suggestion("Is it a correct toml file"),

            _ => Err(eyre!(format!("Unknown {name} file extension: {}", format)))
                .suggestion("Is it a .json or .yaml file?"),
        }
    }

    /// Return the editor schema reference, if this type defines one.
    /// The reference is included by [`Self::to_toml`], as shown in its example.
    fn schema_path(&self) -> Option<PathBuf> {
        None
    }
}

impl<T> SmartSerde for BTreeMap<String, (T, T)> where
    T: Clone + From<f64> + Serialize + DeserializeOwned
{
}

impl SmartSerde for Model {}
// impl SmartSerde for Schema {}

impl SmartSerde for GlobalSettings {
    fn schema_path(&self) -> Option<PathBuf> {
        Some(get_schema_folder().join("global.json"))
    }
}
impl SmartSerde for RuntimeSettings {
    fn schema_path(&self) -> Option<PathBuf> {
        Some(get_schema_folder().join("runtime.json"))
    }
}

/// Resolve the schema reference without creating directories or fetching schemas.
/// Uses `GAMMALOOP_SCHEMA_PATH` when set, otherwise the versioned online reference.
/// See [`SmartSerde::to_toml`] for its use in serialized settings.
pub fn get_schema_folder() -> PathBuf {
    env::var("GAMMALOOP_SCHEMA_PATH")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            let schema_ref = schema_ref();
            PathBuf::from(format!(
                "https://raw.githubusercontent.com/alphal00p/gammaloop/{schema_ref}/assets/schemas"
            ))
        })
}

/// Create the local schema export directory, honoring `GAMMALOOP_SCHEMA_PATH`.
///
/// # Errors
/// Fails if the home directory cannot be determined or the directory cannot be created.
///
/// ```no_run
/// # use color_eyre::Report;
/// # use gammalooprs::utils::serde_utils;
/// let directory = serde_utils::get_local_schema_folder()?;
/// assert!(directory.is_dir());
/// # Ok::<(), Report>(())
/// ```
pub fn get_local_schema_folder() -> Result<PathBuf> {
    let folder = match env::var("GAMMALOOP_SCHEMA_PATH") {
        Ok(path) => PathBuf::from(path),
        Err(_) => match home_dir() {
            Some(home) => home.join(".config").join("gammaloop").join("schemas"),
            None => {
                return Err(eyre!("Could not determine home directory"))
                    .with_suggestion(|| "Set the GAMMALOOP_SCHEMA_PATH environment variable");
            }
        },
    };

    if !folder.exists() {
        fs::create_dir_all(&folder).wrap_err_with(|| {
            format!(
                "Could not create schema folder at {}",
                folder.to_string_lossy()
            )
        })?;
    }

    Ok(folder)
}

use crate::{
    model::Model,
    settings::{GlobalSettings, RuntimeSettings},
    utils::F,
};

pub trait IsDefault {
    fn is_default(&self) -> bool;
}

pub static SHOWDEFAULTS: AtomicBool = AtomicBool::new(false);

pub(crate) fn show_defaults_helper(condition: bool) -> bool {
    if SHOWDEFAULTS.load(Ordering::Relaxed) {
        false
    } else {
        condition
    }
}

pub struct ShowDefaultsGuard {
    _lock: MutexGuard<'static, ()>,
    previous: bool,
}

impl ShowDefaultsGuard {
    pub fn new(show_defaults: bool) -> Self {
        static SHOWDEFAULTS_MUTEX: Mutex<()> = Mutex::new(());
        let lock = SHOWDEFAULTS_MUTEX
            .lock()
            .expect("SHOWDEFAULTS serialization mutex must not be poisoned");
        let previous = SHOWDEFAULTS.swap(show_defaults, Ordering::Relaxed);
        Self {
            _lock: lock,
            previous,
        }
    }
}

impl Drop for ShowDefaultsGuard {
    fn drop(&mut self) {
        SHOWDEFAULTS.store(self.previous, Ordering::Relaxed);
    }
}

impl<T: Default + PartialEq> IsDefault for T {
    fn is_default(&self) -> bool {
        show_defaults_helper(self == &T::default())
    }
}

pub fn is_default_pysecdec_relative_precision(val: &f64) -> bool {
    show_defaults_helper(*val == 1.0e-7_f64)
}

pub fn is_default_esurface_existence_threshold(val: &f64) -> bool {
    show_defaults_helper(*val == super::DEFAULT_ESURFACE_EXISTENCE_THRESHOLD)
}

pub fn is_default_vakint_normalization(val: &String) -> bool {
    show_defaults_helper(val == "MSbar")
}

pub fn is_one_string(val: &String) -> bool {
    show_defaults_helper(val == "1")
}

pub fn is_float<const D: i64>(val: &f64) -> bool {
    show_defaults_helper(*val == D as f64)
}
pub fn is_ffloat<const D: i64>(val: &F<f64>) -> bool {
    show_defaults_helper(*val == F(D as f64))
}

pub fn is_usize<const D: usize>(val: &usize) -> bool {
    show_defaults_helper(*val == D)
}

pub fn is_u64<const D: u64>(val: &u64) -> bool {
    show_defaults_helper(*val == D)
}

pub fn is_false(val: &bool) -> bool {
    show_defaults_helper(!*val)
}

pub fn is_true(val: &bool) -> bool {
    show_defaults_helper(*val)
}

pub fn is_default_input_rescaling(input_rescaling: &Vec<Vec<(f64, f64)>>) -> bool {
    show_defaults_helper(input_rescaling == &_default_input_rescaling())
}

pub fn is_default_shifts(shifts: &Vec<(f64, f64, f64, f64)>) -> bool {
    show_defaults_helper(shifts == &_default_shifts())
}

pub fn _default_input_rescaling() -> Vec<Vec<(f64, f64)>> {
    vec![vec![(0.0, 1.0); 3]; 15]
}
pub fn _default_shifts() -> Vec<(f64, f64, f64, f64)> {
    vec![(1.0, 0.0, 0.0, 0.0); 15]
}

pub fn is_default_form_path(form_path: &String) -> bool {
    show_defaults_helper(form_path == &_default_form_path())
}

pub fn _default_form_path() -> String {
    "form".to_string()
}

pub fn is_default_python_path(python_path: &String) -> bool {
    show_defaults_helper(python_path == &_default_python_path())
}

pub fn _default_python_path() -> String {
    "python3".to_string()
}

pub fn is_default_vakint_evaluation_methods(evaluation_methods: &Vec<String>) -> bool {
    show_defaults_helper(evaluation_methods == &_default_vakint_evaluation_methods())
}

pub fn _default_vakint_evaluation_methods() -> Vec<String> {
    vec![
        "alphaloop".to_string(),
        "matad".to_string(),
        "fmft".to_string(),
    ]
}

pub fn _default_stability_levels() -> Vec<crate::settings::runtime::StabilityLevelSetting> {
    vec![
        crate::settings::runtime::StabilityLevelSetting::default_double(),
        crate::settings::runtime::StabilityLevelSetting::default_quad(),
        crate::settings::runtime::StabilityLevelSetting::default_arb(),
    ]
}

pub fn is_default_stability_levels(
    levels: &Vec<crate::settings::runtime::StabilityLevelSetting>,
) -> bool {
    show_defaults_helper(levels == &_default_stability_levels())
}

pub fn _default_rotation_axis() -> Vec<crate::settings::runtime::RotationSetting> {
    vec![crate::settings::runtime::RotationSetting::EulerAngles {
        alpha: 0.1,
        beta: 0.2,
        gamma: 0.3,
    }]
}

pub fn is_default_rotation_axis(
    rotation_axis: &Vec<crate::settings::runtime::RotationSetting>,
) -> bool {
    show_defaults_helper(rotation_axis == &_default_rotation_axis())
}

#[cfg(test)]
mod tests {
    use crate::utils::{load_generic_model, test_utils::output_dir};
    use serde::{Deserialize, Serialize};
    use std::collections::BTreeMap;
    use std::fs;
    use std::path::PathBuf;

    use super::{IsDefault, SerdeDefaults, SerdeFileError, SmartSerde};

    #[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
    struct SettingsFixture {
        #[serde(default, skip_serializing_if = "IsDefault::is_default")]
        enabled: bool,
    }

    impl SmartSerde for SettingsFixture {
        fn schema_path(&self) -> Option<PathBuf> {
            Some(PathBuf::from("https://example.invalid/settings.json"))
        }
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct InvalidToml {
        values: BTreeMap<Vec<u8>, bool>,
    }

    impl SmartSerde for InvalidToml {}

    #[test]
    fn smart_toml_preserves_schema_and_explicit_default_visibility() {
        let settings = SettingsFixture::default();
        let included = settings.to_toml(SerdeDefaults::Include).unwrap();
        let omitted = settings.to_toml(SerdeDefaults::Omit).unwrap();
        for document in [&included, &omitted] {
            assert!(document.starts_with("#:schema https://example.invalid/settings.json\n"));
            assert_eq!(document.matches("#:schema").count(), 1);
            assert_eq!(
                toml::from_str::<SettingsFixture>(document).unwrap(),
                settings
            );
        }
        assert!(included.contains("enabled = false"));
        assert!(!omitted.contains("enabled"));

        let without_schema = BTreeMap::from([("entry".to_owned(), (1.0, 2.0))]);
        let document = without_schema.to_toml(SerdeDefaults::Include).unwrap();
        assert!(!document.contains("#:schema"));
    }

    #[test]
    fn smart_toml_releases_defaults_guard_on_encoder_error() {
        let invalid = InvalidToml {
            values: BTreeMap::from([(vec![1], true)]),
        };
        assert!(invalid.to_toml(SerdeDefaults::Include).is_err());
        let omitted = SettingsFixture::default()
            .to_toml(SerdeDefaults::Omit)
            .unwrap();
        assert!(!omitted.contains("enabled"));
    }

    #[test]
    fn test_file_vs_parse_errors() {
        // Test file not found error
        let result: Result<BTreeMap<String, (f64, f64)>, SerdeFileError> =
            SmartSerde::from_file_typed("/nonexistent/file.json");

        match result {
            Err(SerdeFileError::FileError(_)) => {
                // This is the expected file error
            }
            other => panic!("Expected FileError, got: {:?}", other),
        }

        // Test parse error with invalid JSON
        let temp_path = std::env::temp_dir().join("test_invalid.json");
        fs::write(&temp_path, "{ invalid json content").unwrap();

        let result: Result<BTreeMap<String, (f64, f64)>, SerdeFileError> =
            SmartSerde::from_file_typed(&temp_path);

        match result {
            Err(SerdeFileError::JsonParseError(_)) => {
                // This is the expected parse error
            }
            other => panic!("Expected JsonParseError, got: {:?}", other),
        }

        // Clean up
        let _ = fs::remove_file(&temp_path);

        // Test successful parsing
        let temp_path = std::env::temp_dir().join("test_valid.json");
        fs::write(&temp_path, r#"{"key": [1.0, 2.0]}"#).unwrap();

        let result: Result<BTreeMap<String, (f64, f64)>, SerdeFileError> =
            SmartSerde::from_file_typed(&temp_path);

        assert!(result.is_ok());
        let data = result.unwrap();
        assert_eq!(data.get("key"), Some(&(1.0, 2.0)));

        // Clean up
        let _ = fs::remove_file(&temp_path);
    }

    #[test]
    fn test_from_file_backward_compatibility() {
        // Test that the regular from_file function still works and returns color_eyre::Result
        let temp_path = std::env::temp_dir().join("test_compat.json");
        fs::write(&temp_path, r#"{"test": [1.0, 2.0]}"#).unwrap();

        let result: color_eyre::Result<BTreeMap<String, (f64, f64)>> =
            SmartSerde::from_file(&temp_path, "test");

        assert!(result.is_ok());
        let data = result.unwrap();
        assert_eq!(data.get("test"), Some(&(1.0, 2.0)));

        // Clean up
        let _ = fs::remove_file(&temp_path);
    }

    #[test]
    fn test_to_file_creates_missing_parent_directories() {
        let unique_dir = std::env::temp_dir().join(format!(
            "gammaloop-smart-serde-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let file_path = unique_dir.join("nested").join("data.json");

        let mut data = BTreeMap::new();
        data.insert("test".to_string(), (1.0, 2.0));

        data.to_file(&file_path, true).unwrap();

        assert!(file_path.exists());
        let roundtrip: BTreeMap<String, (f64, f64)> =
            SmartSerde::from_file(&file_path, "test").unwrap();
        assert_eq!(roundtrip.get("test"), Some(&(1.0, 2.0)));

        let settings = SettingsFixture::default();
        let settings_path = unique_dir.join("nested").join("settings.toml");
        settings.to_file(&settings_path, true).unwrap();
        assert_eq!(
            fs::read_to_string(&settings_path).unwrap(),
            settings.to_toml(SerdeDefaults::Include).unwrap()
        );
        assert_eq!(
            SettingsFixture::from_file_typed(&settings_path).unwrap(),
            settings
        );

        fs::remove_dir_all(&unique_dir).unwrap();
    }

    mod failing {
        use super::*;

        #[test]
        fn convert_models() {
            let name = "scalars";
            load_generic_model(name)
                .to_file(
                    output_dir().join(format!("gammaloop_models/{name}.json")),
                    true,
                )
                .unwrap();
        }
    }
}
