use serde::{de::DeserializeOwned, Serialize};
use std::{
    io::Write,
    sync::atomic::{AtomicBool, Ordering},
};

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
                        f.write(toml::to_string_pretty(&self)?.as_bytes())?;
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
                        let mut buf = vec![];
                        f.read(&mut buf)?;
                        toml::from_slice(&buf)
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
}

impl<T: Serialize + DeserializeOwned> SmartSerde for T {}

use crate::utils::F;

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
