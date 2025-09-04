use dirs::home_dir;
use schemars::{schema_for, JsonSchema, Schema};
use serde::{de::DeserializeOwned, Serialize};
use std::{
    env,
    fs::OpenOptions,
    io::Write,
    path::PathBuf,
    sync::atomic::{AtomicBool, Ordering},
};
use toml::{ser::Buffer, Serializer};

use color_eyre::{Result, Section};
use eyre::{eyre, Context};

use std::{fs::File, io::Read, path::Path};

pub trait SmartSerde: Serialize + DeserializeOwned {
    fn to_file(&self, file_path: impl AsRef<Path>) -> Result<()> {
        let mut f = File::create(file_path.as_ref())?;

        if let Some(ext) = file_path.as_ref().extension() {
            if let Some(ext) = ext.to_str() {
                match ext {
                    "json" => serde_json::to_writer_pretty(f, self)
                        .map_err(|e| eyre!(format!("Error serializing json: {}", e)))?,
                    "yaml" | "yml" => serde_yaml::to_writer(f, self)
                        .map_err(|e| eyre!(format!("Error serializing yaml: {}", e)))?,
                    "toml" => {
                        let mut toml_string = if let Some(schema_path) = self.has_schema_path() {
                            let schema_path = schema_path?;
                            format!("#:schema {}\n", schema_path.display())
                        } else {
                            String::new()
                        };

                        toml_string.push_str(&toml::to_string_pretty(&self)?);

                        f.write(toml_string.as_bytes())?;
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

    fn from_file(file_path: impl AsRef<Path>, name: &str) -> Result<Self> {
        let mut f = File::open(file_path.as_ref())
            .wrap_err_with(|| {
                format!(
                    "Could not open {name} file {}",
                    file_path.as_ref().display()
                )
            })
            .suggestion("Does the path exist?")?;

        if let Some(ext) = file_path.as_ref().extension() {
            if let Some(ext) = ext.to_str() {
                match ext {
                    "json" => serde_json::from_reader(f)
                        .map_err(|e| eyre!(format!("Error parsing {name} json: {}", e)))
                        .suggestion("Is it a correct json file"),
                    "yaml" | "yml" => serde_yaml::from_reader(f)
                        .map_err(|e| eyre!(format!("Error parsing {name} yaml: {}", e)))
                        .suggestion("Is it a correct yaml file"),
                    "toml" => {
                        let mut buf = String::new();
                        let bytes = f.read_to_string(&mut buf)?;
                        toml::from_str(&buf)
                            .map_err(|e| eyre!(format!("Error parsing {name} toml: {}", e)))
                            .suggestion("Is it a correct toml file")
                    }

                    _ => Err(eyre!(format!("Unknown {name} file extension: {}", ext)))
                        .suggestion("Is it a .json or .yaml file?"),
                }
            } else {
                Err(eyre!(format!(
                    "Could not determine file extension of {name} file {}",
                    file_path.as_ref().display()
                )))
                .suggestion("Does the path exist?")
            }
        } else {
            Err(eyre!(format!(
                "Could not determine file extension of {name} file {}",
                file_path.as_ref().display()
            )))
            .suggestion("Does the path exist?")
        }
    }

    fn has_schema_path(&self) -> Option<Result<PathBuf>> {
        Option::None
    }
}

impl SmartSerde for SerializableModel {}
// impl SmartSerde for Schema {}

impl SmartSerde for RunHistory {
    fn has_schema_path(&self) -> Option<Result<PathBuf>> {
        Some(get_schema_folder().map(|f| f.join("runhistory.json")))
    }
}
impl SmartSerde for GlobalSettings {
    fn has_schema_path(&self) -> Option<Result<PathBuf>> {
        Some(get_schema_folder().map(|f| f.join("global.json")))
    }
}
impl SmartSerde for RuntimeSettings {
    fn has_schema_path(&self) -> Option<Result<PathBuf>> {
        Some(get_schema_folder().map(|f| f.join("runtime.json")))
    }
}

use crate::{
    cli::state::RunHistory,
    model::{Model, SerializableModel},
    settings::{GlobalSettings, RuntimeSettings},
    utils::F,
};

pub trait IsDefault {
    fn is_default(&self) -> bool;
}

pub static SHOWDEFAULTS: AtomicBool = AtomicBool::new(false);

fn show_defaults_helper(condition: bool) -> bool {
    if SHOWDEFAULTS.load(Ordering::Relaxed) {
        false
    } else {
        condition
    }
}

impl<T: Default + PartialEq> IsDefault for T {
    fn is_default(&self) -> bool {
        show_defaults_helper(self == &T::default())
    }
}

pub fn is_float<const D: i64>(val: &F<f64>) -> bool {
    show_defaults_helper(*val == F(D as f64))
}

pub fn is_usize<const D: usize>(val: &usize) -> bool {
    show_defaults_helper(*val == D)
}

pub fn is_u64<const D: u64>(val: &u64) -> bool {
    show_defaults_helper(*val == D)
}

pub fn is_false(val: &bool) -> bool {
    show_defaults_helper(*val == false)
}

pub fn is_true(val: &bool) -> bool {
    show_defaults_helper(*val == true)
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

pub fn _default_stability_levels() -> Vec<crate::settings::runtime::StabilityLevelSetting> {
    vec![
        crate::settings::runtime::StabilityLevelSetting::default_double(),
        crate::settings::runtime::StabilityLevelSetting::default_quad(),
    ]
}

pub fn is_default_stability_levels(
    levels: &Vec<crate::settings::runtime::StabilityLevelSetting>,
) -> bool {
    show_defaults_helper(levels == &_default_stability_levels())
}

pub fn _default_rotation_axis() -> Vec<crate::settings::runtime::RotationSetting> {
    vec![crate::settings::runtime::RotationSetting::default()]
}

pub fn is_default_rotation_axis(
    rotation_axis: &Vec<crate::settings::runtime::RotationSetting>,
) -> bool {
    show_defaults_helper(rotation_axis == &_default_rotation_axis())
}

fn get_schema_folder() -> Result<PathBuf> {
    let folder = match env::var("GAMMALOOP_SCHEMA_PATH") {
        Ok(path) => PathBuf::try_from(path)?,
        Err(_) => match home_dir() {
            Some(home) => home.join(".config").join("gammaloop").join("schemas"),
            None => {
                return Err(eyre!("Could not determine home directory"))
                    .with_suggestion(|| "Set the GAMMALOOP_SCHEMA_PATH environment variable");
            }
        },
    };

    if !folder.exists() {
        std::fs::create_dir_all(&folder).wrap_err_with(|| {
            format!(
                "Could not create schema folder at {}",
                folder.to_string_lossy()
            )
        })?;
    }

    Ok(folder)
}

pub fn write_schemas() -> Result<()> {
    let global_schema = schema_for!(GlobalSettings);
    let runtime_schema = schema_for!(RuntimeSettings);
    let runhistory_schema = schema_for!(RunHistory);
    let folder = get_schema_folder()?;

    let mut global_file = File::create(folder.join("global.json"))?;
    let mut runtime_file = File::create(folder.join("runtime.json"))?;
    let mut runhistory_file = File::create(folder.join("runhistory.json"))?;

    serde_json::to_writer_pretty(&mut global_file, &global_schema)
        .wrap_err("Could not write global schema")?;
    serde_json::to_writer_pretty(&mut runtime_file, &runtime_schema)
        .wrap_err("Could not write runtime schema")?;
    serde_json::to_writer_pretty(&mut runhistory_file, &runhistory_schema)
        .wrap_err("Could not write runhistory schema")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::utils::test_utils::{load_generic_model, output_dir};

    use super::SmartSerde;

    #[test]
    fn convert_models() {
        let name = "scalars";
        load_generic_model(name)
            .to_serializable()
            .to_file(output_dir().join(format!("gammaloop_models/{name}.json")))
            .unwrap();
    }
}
