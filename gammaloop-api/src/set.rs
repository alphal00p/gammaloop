use std::{collections::HashMap, path::PathBuf, str::FromStr};

use clap::{command, ArgGroup, Args, Subcommand};
use color_eyre::Result;
use eyre::{eyre, Context};
use figment::{
    providers::{Format, Serialized},
    Figment,
};
use gammalooprs::{
    processes::ProcessCollection,
    settings::{GlobalSettings, RuntimeSettings},
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{
    generate::ProcessArgs,
    state::{State, SyncSettings},
};

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Set {
    /// Set the base directory for gammaloop state and logs
    BaseDir { path: PathBuf },
    /// Set GLOBAL settings
    Global {
        #[command(subcommand)]
        input: SetArgs,
    },

    /// Set DEFAULT RUNTIME settings
    DefaultRuntime {
        #[command(subcommand)]
        input: SetArgs,
    },

    /// Set settings for a PROCESS
    Process {
        #[command(subcommand)]
        input: SetArgs,
        #[command(flatten)]
        process: ProcessArgs,
    },
}

impl Set {
    pub fn run(
        &self,
        state: &mut State,
        global_settings: &mut GlobalSettings,
        default_runtime_settings: &mut RuntimeSettings,
        state_folder: &mut PathBuf,
    ) -> Result<()> {
        match self {
            Set::BaseDir { path } => {
                *state_folder = path.clone();
            }
            Self::Global { input } => {
                let fig = Figment::from(Serialized::defaults(&global_settings));

                *global_settings = input.merge_figment(fig)?.extract()?;
                global_settings.sync_settings()?;
            }
            Self::DefaultRuntime { input } => {
                let fig = Figment::from(Serialized::defaults(&default_runtime_settings));

                *default_runtime_settings = input.merge_figment(fig)?.extract()?;
            }
            Self::Process { input, process } => {
                let p = &mut state.process_list.processes[process.process_id];
                if let Some(name) = &process.name {
                    match &mut p.collection {
                        ProcessCollection::Amplitudes(a) => {
                            let Some(a) = a.get_mut(name) else {
                                return Err(eyre!(
                                    "No amplitude named '{}' in process id {}",
                                    name,
                                    process.process_id
                                ));
                            };
                            let Some(a) = &mut a.integrand else {
                                return Err(eyre!(
                                    "No integrand in amplitude '{}' in process id {}",
                                    name,
                                    process.process_id
                                ));
                            };

                            let settings = a.get_mut_settings();
                            let fig = Figment::from(Serialized::defaults(&settings));

                            *settings = input.merge_figment(fig)?.extract()?;
                        }
                        ProcessCollection::CrossSections(a) => {
                            // a[name].preprocess(&self.model, &global_settings.generation)?;
                        }
                    }
                } else {
                    match &mut p.collection {
                        ProcessCollection::Amplitudes(a) => {
                            for (_, amp) in a.iter_mut() {
                                if let Some(a) = &mut amp.integrand {
                                    let settings = a.get_mut_settings();
                                    let fig = Figment::from(Serialized::defaults(&settings));

                                    *settings = input.merge_figment(fig)?.extract()?;
                                };
                            }
                        }
                        ProcessCollection::CrossSections(a) => {
                            // a[name].preprocess(&self.model, &global_settings.generation)?;
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

// Shared input forms for global & default-runtime
#[derive(Subcommand, Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub enum SetArgs {
    /// Load from a settings file
    File { file: PathBuf },

    /// Set one or more dotted key-paths
    Kv {
        /// Any number of KEY=VALUE pairs
        #[arg(value_name = "KEY=VALUE", num_args = 1.., value_parser = KvPair::from_str)]
        pairs: Vec<KvPair>,
    },

    /// Use the stored settings file
    Stored,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct KvPair {
    pub key: String,
    pub value: String,
}

impl FromStr for KvPair {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (k, v) = s.split_once('=').ok_or("Expected KEY=VALUE")?;
        if k.trim().is_empty() {
            return Err("Empty KEY".into());
        }
        Ok(KvPair {
            key: k.trim().into(),
            value: v.trim().into(),
        })
    }
}
impl SetArgs {
    pub fn merge_figment(&self, fig: Figment) -> Result<Figment> {
        match self {
            SetArgs::File { file } => {
                let ext = file
                    .extension()
                    .and_then(|s| s.to_str())
                    .unwrap_or("")
                    .to_ascii_lowercase();
                match ext.as_str() {
                    "toml" => Ok(fig.merge(figment::providers::Toml::file(file))),
                    "yaml" | "yml" => Ok(fig.merge(figment::providers::Yaml::file(file))),
                    "json" => Ok(fig.merge(figment::providers::Json::file(file))),
                    _ => Err(color_eyre::eyre::eyre!(
                        "Unsupported settings file extension: {}",
                        ext
                    )),
                }
            }
            SetArgs::Kv { pairs } => {
                let mut figment = Figment::new();

                for KvPair { key, value } in pairs {
                    figment = figment.adjoin((key.clone(), infer_cli_value(value)?));
                }

                Ok(fig.merge(figment))
            }
            _ => Ok(fig),
        }
    }
}

use serde_json::Value as J;
use serde_yaml::Value as Y;

fn infer_cli_value(raw: &str) -> Result<J> {
    // 1) Try strict JSON (numbers/bools work; strings need quotes)
    if let Ok(v) = serde_json::from_str::<J>(raw) {
        return Ok(v);
    }
    // 2) YAML scalars: true/false, 42, 3.14, and bare words -> string
    if let Ok(v) = serde_yaml::from_str::<Y>(raw) {
        return yaml_to_json(v).with_context(|| format!("While trying yaml on :{raw}"));
    }
    // 3) Fallback: plain string
    Ok(J::String(raw.to_string()))
}

fn yaml_to_json(v: Y) -> Result<J> {
    match v {
        Y::Null => Ok(J::Null),
        Y::Bool(b) => Ok(J::Bool(b)),
        Y::Number(n) => {
            if let Some(i) = n.as_i64() {
                Ok(J::Number(i.into()))
            } else if let Some(f) = n.as_f64() {
                Ok(J::Number(serde_json::Number::from_f64(f).unwrap()))
            } else {
                Ok(J::String(n.to_string()))
            }
        }
        Y::String(s) => Ok(J::String(s)),
        other => Err(eyre!(
            "Unsupported complex YAML value in CLI key-value: {:?}",
            other
        )),
    }
}
