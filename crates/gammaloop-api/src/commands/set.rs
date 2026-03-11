use std::{path::PathBuf, str::FromStr};

use clap::Subcommand;
use color_eyre::Result;
use eyre::{eyre, Context, Report};
use figment::{
    providers::{Format, Serialized},
    Figment,
};
use gammalooprs::{
    model::UFOSymbol, processes::ProcessCollection, settings::RuntimeSettings, utils::F,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use tracing::warn;

use crate::{
    commands::generate::ProcessArgs,
    commands::Commands,
    state::{State, SyncSettings},
    CLISettings,
};

impl FromStr for Set {
    type Err = Report;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if let Commands::Set(cmd) = Commands::from_str(s)? {
            Ok(cmd)
        } else {
            Err(eyre!("Not a 'set' command"))
        }
    }
}

#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Set {
    /// Set the base directory for gammaloop state and logs
    BaseDir {
        #[arg(value_hint = clap::ValueHint::DirPath)]
        path: PathBuf,
    },
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

    /// Set Model parameters
    Model {
        /// Any number of KEY=VALUE pairs
        #[arg(value_name = "KEY=VALUE", num_args = 1.., value_parser = KvPair::from_str)]
        pairs: Vec<KvPair>,
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
        global_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<()> {
        match self {
            Set::Model { pairs } => {
                for KvPair { key, value } in pairs.iter() {
                    let value = json5::from_str::<Complex<F<f64>>>(value).with_context(|| {
                        format!("While parsing JSON value {value} for key '{key}'")
                    })?;

                    if let Some(p) = state
                        .model_parameters
                        .get_mut(&UFOSymbol::from(key.as_str()))
                    {
                        *p = value;
                        continue;
                    } else {
                        let possiblilities: Vec<String> = state
                            .model_parameters
                            .keys()
                            .map(|s| s.to_string())
                            .collect();
                        return Err(eyre!("No model parameter named '{key}'")).with_context(|| {
                            format!(
                                "Possible model parameters are: {}",
                                possiblilities.join(", ")
                            )
                        });
                    }
                }
                state.model_parameters.apply_to_model(&mut state.model)?;
            }
            Set::BaseDir { path } => {
                warn!(
                    "Ignoring base-dir change request to '{}': state folder is fixed for the current session ('{}')",
                    path.display(),
                    global_settings.state_folder.display()
                );
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
                let process_id = state.resolve_process_ref(process.process.as_ref())?;
                let default_runtime_template =
                    matches!(input, SetArgs::Defaults).then(|| default_runtime_settings.clone());
                let apply_runtime_settings = |settings: &mut RuntimeSettings| -> Result<()> {
                    if let Some(template) = default_runtime_template.as_ref() {
                        *settings = template.clone();
                    } else {
                        let fig = Figment::from(Serialized::defaults(&settings));
                        *settings = input.merge_figment(fig)?.extract()?;
                    }
                    Ok(())
                };

                if let Some(name) = &process.integrand_name {
                    let integrand = state.process_list.get_integrand_mut(process_id, name)?;

                    let settings = integrand.get_mut_settings();
                    apply_runtime_settings(settings)?;
                } else {
                    match &mut state.process_list.processes[process_id].collection {
                        ProcessCollection::Amplitudes(a) => {
                            for (_, amp) in a.iter_mut() {
                                if let Some(a) = &mut amp.integrand {
                                    let settings = a.get_mut_settings();
                                    apply_runtime_settings(settings)?;
                                };
                            }
                        }
                        ProcessCollection::CrossSections(a) => {
                            for (_, xs) in a.iter_mut() {
                                if let Some(a) = &mut xs.integrand {
                                    let settings = a.get_mut_settings();
                                    apply_runtime_settings(settings)?;
                                };
                            }
                            // a[name].preprocess(&self.model, &global_settings.generation)?;
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

// Shared input forms for set commands
#[derive(Subcommand, Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub enum SetArgs {
    /// Load from a settings file
    File {
        #[arg(value_hint = clap::ValueHint::FilePath)]
        file: PathBuf,
    },

    /// Load settings from TOML content passed as a CLI string
    String {
        #[arg(value_name = "TOML")]
        string: String,
    },

    /// Set one or more dotted key-paths
    Kv {
        /// Any number of KEY=VALUE pairs
        #[arg(value_name = "KEY=VALUE", num_args = 1.., value_parser = KvPair::from_str)]
        pairs: Vec<KvPair>,
    },

    /// Sync process runtime settings from the current default runtime settings
    Defaults,

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
                    // "yaml" | "yml" => Ok(fig.merge(figment::providers::Yaml::file(file))),
                    "json" => Ok(fig.merge(figment::providers::Json::file(file))),
                    _ => Err(color_eyre::eyre::eyre!(
                        "Unsupported settings file extension: {}",
                        ext
                    )),
                }
            }
            SetArgs::String { string } => Ok(fig.merge(figment::providers::Toml::string(string))),
            SetArgs::Kv { pairs } => {
                let mut figment = Figment::new();

                for KvPair { key, value } in pairs.iter() {
                    figment = figment.adjoin((key.clone(), infer_cli_value(value)?));
                }

                Ok(fig.merge(figment))
            }
            SetArgs::Stored => Ok(fig),
            SetArgs::Defaults => Err(eyre!("'defaults' is only supported for 'set process'")),
        }
    }
}

use serde_json::Value as J;
use serde_yaml::Value as Y;

fn infer_cli_value(raw: &str) -> Result<J> {
    // 1) Try strict JSON (numbers/bools work; strings need quotes)
    if let Ok(v) = json5::from_str::<J>(raw) {
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

#[cfg(test)]
mod test {
    use std::str::FromStr;

    use clap::Parser;
    use figment::{providers::Serialized, Figment};
    use gammalooprs::{
        model::{ParameterNature, UFOSymbol},
        settings::RuntimeSettings,
        utils::{test_utils::load_generic_model, F},
    };
    use serde::{Deserialize, Serialize};
    use spenso::algebra::complex::Complex;

    use crate::{
        state::{ProcessRef, State},
        CLISettings, Repl,
    };

    use super::{super::Commands, KvPair, Set, SetArgs};

    #[test]
    fn serialize_complex() {
        let s = Complex::new(F(1.), F(-2.));
        let j = serde_json::to_string(&s).unwrap();
        assert_eq!(j, r#"{"re":1.0,"im":-2.0}"#.to_string());
    }

    #[test]
    fn parse_set_process_defaults() {
        let cmd = Set::from_str("set process -p epem_a_tth -i LO defaults").unwrap();

        match cmd {
            Set::Process { input, process } => {
                assert_eq!(input, SetArgs::Defaults);
                assert_eq!(
                    process.process,
                    Some(ProcessRef::Unqualified("epem_a_tth".to_string()))
                );
                assert_eq!(process.integrand_name, Some("LO".to_string()));
            }
            other => panic!("Expected set process command, got {other:?}"),
        }
    }

    #[test]
    fn defaults_is_rejected_outside_set_process() {
        let err = SetArgs::Defaults.merge_figment(Figment::new()).unwrap_err();
        assert!(err
            .to_string()
            .contains("'defaults' is only supported for 'set process'"));
    }

    #[test]
    fn parse_set_process_string() {
        let repl = Repl::try_parse_from([
            "gammaloop",
            "set",
            "process",
            "-p",
            "epem_a_tth",
            "-i",
            "LO",
            "string",
            "alpha = 2",
        ])
        .unwrap();

        let cmd = match repl.command {
            Commands::Set(set) => set,
            other => panic!("Expected set command, got {other:?}"),
        };

        match cmd {
            Set::Process { input, process } => {
                assert_eq!(
                    input,
                    SetArgs::String {
                        string: "alpha = 2".to_string()
                    }
                );
                assert_eq!(
                    process.process,
                    Some(ProcessRef::Unqualified("epem_a_tth".to_string()))
                );
                assert_eq!(process.integrand_name, Some("LO".to_string()));
            }
            other => panic!("Expected set process command, got {other:?}"),
        }
    }

    #[test]
    fn parse_set_process_string_multiline() {
        let multiline = "alpha = 2\nbeta = true\n";
        let repl = Repl::try_parse_from([
            "gammaloop",
            "set",
            "process",
            "-p",
            "epem_a_tth",
            "-i",
            "LO",
            "string",
            multiline,
        ])
        .unwrap();

        let cmd = match repl.command {
            Commands::Set(set) => set,
            other => panic!("Expected set command, got {other:?}"),
        };

        match cmd {
            Set::Process { input, .. } => {
                assert_eq!(
                    input,
                    SetArgs::String {
                        string: multiline.to_string()
                    }
                );
            }
            other => panic!("Expected set process command, got {other:?}"),
        }
    }

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    struct MergeFixture {
        alpha: u64,
        beta: bool,
    }

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    struct CompileSettingsFixture {
        inline_asm: bool,
        optimization_level: String,
        fast_math: bool,
        unsafe_math: bool,
        custom: Vec<String>,
    }

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    struct NCoresFixture {
        feyngen: u64,
        generate: u64,
        compile: u64,
        integrate: u64,
    }

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    struct GenerationFixture {
        compile: CompileSettingsFixture,
    }

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    struct GlobalFixture {
        generation: GenerationFixture,
        n_cores: NCoresFixture,
    }

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    struct CliSettingsFixture {
        global: GlobalFixture,
    }

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    struct RunLikeFixture {
        cli_settings: CliSettingsFixture,
    }

    #[test]
    fn string_toml_merges_like_file_toml() {
        let defaults = MergeFixture {
            alpha: 1,
            beta: false,
        };
        let toml = "alpha = 123\nbeta = true\n";

        let tmp_name = format!(
            "gammaloop_set_args_test_{}_{}.toml",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        );
        let file_path = std::env::temp_dir().join(tmp_name);
        std::fs::write(&file_path, toml).unwrap();

        let base_a = Figment::from(Serialized::defaults(&defaults));
        let via_file: MergeFixture = SetArgs::File {
            file: file_path.clone(),
        }
        .merge_figment(base_a)
        .unwrap()
        .extract()
        .unwrap();

        let base_b = Figment::from(Serialized::defaults(&defaults));
        let via_string: MergeFixture = SetArgs::String {
            string: toml.to_string(),
        }
        .merge_figment(base_b)
        .unwrap()
        .extract()
        .unwrap();

        let _ = std::fs::remove_file(file_path);
        assert_eq!(via_string, via_file);
    }

    #[test]
    fn string_toml_supports_multiline_nested_tables() {
        let defaults = RunLikeFixture {
            cli_settings: CliSettingsFixture {
                global: GlobalFixture {
                    generation: GenerationFixture {
                        compile: CompileSettingsFixture {
                            inline_asm: false,
                            optimization_level: "O0".to_string(),
                            fast_math: false,
                            unsafe_math: false,
                            custom: vec!["-g".to_string()],
                        },
                    },
                    n_cores: NCoresFixture {
                        feyngen: 1,
                        generate: 1,
                        compile: 1,
                        integrate: 1,
                    },
                },
            },
        };

        let multiline = r#"
[cli_settings.global.generation.compile]
inline_asm = true
optimization_level = "O3"
fast_math = true
unsafe_math = true
custom = []

[cli_settings.global.n_cores]
feyngen = 10
generate = 1
compile = 10
integrate = 10
"#;

        let base = Figment::from(Serialized::defaults(&defaults));
        let merged: RunLikeFixture = SetArgs::String {
            string: multiline.to_string(),
        }
        .merge_figment(base)
        .unwrap()
        .extract()
        .unwrap();

        assert!(merged.cli_settings.global.generation.compile.inline_asm);
        assert_eq!(
            merged
                .cli_settings
                .global
                .generation
                .compile
                .optimization_level,
            "O3"
        );
        assert_eq!(merged.cli_settings.global.n_cores.feyngen, 10);
        assert_eq!(merged.cli_settings.global.n_cores.generate, 1);
        assert_eq!(merged.cli_settings.global.n_cores.compile, 10);
        assert_eq!(merged.cli_settings.global.n_cores.integrate, 10);
    }

    #[test]
    fn set_model_rejects_internal_parameters() {
        let model = load_generic_model("sm");
        let internal_param = model
            .parameters
            .values()
            .find(|param| param.nature == ParameterNature::Internal)
            .expect("SM model must contain internal parameters")
            .name
            .to_string();

        let mut state = State::new_test();
        state.model = model;
        state.model_parameters =
            gammalooprs::model::InputParamCard::default_from_model(&state.model);
        assert!(!state
            .model_parameters
            .contains_key(&UFOSymbol::from(internal_param.as_str())));

        let err = Set::Model {
            pairs: vec![KvPair {
                key: internal_param.clone(),
                value: "1.0".to_string(),
            }],
        }
        .run(
            &mut state,
            &mut CLISettings::default(),
            &mut RuntimeSettings::default(),
        )
        .unwrap_err();

        let err_text = format!("{err:?}");
        assert!(!err_text.is_empty());
    }
}
