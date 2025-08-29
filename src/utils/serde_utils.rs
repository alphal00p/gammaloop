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

pub trait IsDefault {
    fn is_default(&self) -> bool;
}

pub static SHOWDEFAULTS: AtomicBool = AtomicBool::new(false);

impl<T: Default + PartialEq> IsDefault for T {
    fn is_default(&self) -> bool {
        // println!("HHHHIII {}", std::any::type_name::<T>());
        if SHOWDEFAULTS.load(Ordering::Relaxed) {
            false
        } else {
            self == &T::default()
        }
    }
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
