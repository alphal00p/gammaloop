use std::{fs, path::PathBuf, str::FromStr};

use clap::Subcommand;
use color_eyre::Result;
use eyre::{eyre, Context, Report};
use figment::{
    providers::{Format, Serialized},
    Figment,
};
use gammalooprs::{
    model::{ParameterNature, ParameterType, UFOSymbol},
    processes::ProcessCollection,
    settings::RuntimeSettings,
    utils::F,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value as JsonValue;
use spenso::algebra::complex::Complex;
use tracing::warn;

use crate::{
    commands::generate::ProcessArgs,
    commands::Commands,
    state::{State, SyncSettings},
    tracing::{clear_file_log_filter_override_on_settings_change, file_log_boot_disabled_reason},
    CLISettings,
};

pub(crate) const MODEL_COMPLEX_VALUE_FORMAT_HINT: &str =
    "expects complex like 1.0, -1.0e13, 1.0+2.0i, or 1.0-2.0j";
pub(crate) const MODEL_REAL_VALUE_FORMAT_HINT: &str = "expects real like 1.0, -1.0, or 1.0e13";

pub(crate) fn model_value_format_hint(parameter_type: Option<ParameterType>) -> &'static str {
    match parameter_type {
        Some(ParameterType::Real) => MODEL_REAL_VALUE_FORMAT_HINT,
        _ => MODEL_COMPLEX_VALUE_FORMAT_HINT,
    }
}

fn parse_model_parameter_value(value: &str) -> Result<Complex<F<f64>>> {
    let value = value.trim();
    if value.is_empty() {
        return Err(eyre!("Model parameter value cannot be empty"));
    }
    if value.starts_with('{') || value.starts_with('[') {
        return Err(eyre!(
            "Legacy model-parameter syntax like '{{re:...,im:...}}' or '[re,im]' is deprecated; use {MODEL_COMPLEX_VALUE_FORMAT_HINT}"
        ));
    }

    let (re, im) = parse_complex_literal(value)?;
    Ok(Complex::new(F(re), F(im)))
}

fn parse_complex_literal(value: &str) -> Result<(f64, f64)> {
    let value = value.trim();
    if value.is_empty() {
        return Err(eyre!("Complex literal cannot be empty"));
    }

    if let Some(unit) = value.chars().last() {
        if matches!(unit, 'i' | 'j' | 'I' | 'J') {
            let body = &value[..value.len() - unit.len_utf8()];
            let body = body.trim();
            if body.is_empty() {
                return Err(eyre!(
                    "Missing coefficient before imaginary unit in '{value}'"
                ));
            }

            if let Some(index) = find_imaginary_separator(body) {
                let re = parse_float_literal(&body[..index], value)?;
                let im = parse_imaginary_literal(&body[index..], value)?;
                return Ok((re, im));
            }

            let im = parse_imaginary_literal(body, value)?;
            return Ok((0.0, im));
        }
    }

    Ok((parse_float_literal(value, value)?, 0.0))
}

fn find_imaginary_separator(value: &str) -> Option<usize> {
    let bytes = value.as_bytes();
    for index in (1..bytes.len()).rev() {
        let current = bytes[index] as char;
        if !matches!(current, '+' | '-') {
            continue;
        }
        let previous = bytes[index - 1] as char;
        if matches!(previous, 'e' | 'E') {
            continue;
        }
        return Some(index);
    }
    None
}

fn parse_float_literal(value: &str, full_value: &str) -> Result<f64> {
    value
        .trim()
        .parse::<f64>()
        .with_context(|| format!("Failed to parse real part of complex literal '{full_value}'"))
}

fn parse_imaginary_literal(value: &str, full_value: &str) -> Result<f64> {
    let value = value.trim();
    if value == "+" {
        return Ok(1.0);
    }
    if value == "-" {
        return Ok(-1.0);
    }
    value.parse::<f64>().with_context(|| {
        format!("Failed to parse imaginary part of complex literal '{full_value}'")
    })
}

fn validate_model_parameter_type(
    parameter_name: &str,
    parameter_type: ParameterType,
    value: &Complex<F<f64>>,
) -> Result<()> {
    if parameter_type == ParameterType::Real && value.im.0 != 0.0 {
        return Err(eyre!(
            "Model parameter '{parameter_name}' is Real and cannot be assigned an imaginary component; {}",
            MODEL_REAL_VALUE_FORMAT_HINT
        ));
    }
    Ok(())
}

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
        /// Any number of PARAM=COMPLEX pairs
        #[arg(value_name = "PARAM=COMPLEX", num_args = 1.., value_parser = KvPair::from_str)]
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
                    let parameter = state
                        .model
                        .get_parameter_opt(key)
                        .filter(|parameter| parameter.nature == ParameterNature::External);
                    if parameter.is_none() {
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
                    let parameter_type = parameter.unwrap().parameter_type.clone();
                    let value = parse_model_parameter_value(value).with_context(|| {
                        format!("While parsing model parameter value {value} for key '{key}'")
                    })?;
                    validate_model_parameter_type(key, parameter_type, &value)?;

                    if let Some(p) = state
                        .model_parameters
                        .get_mut(&UFOSymbol::from(key.as_str()))
                    {
                        *p = value;
                        continue;
                    }
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
                state.model_parameters.apply_to_model(&mut state.model)?;
            }
            Set::BaseDir { path } => {
                warn!(
                    "Ignoring base-dir change request to '{}': state folder is fixed for the current session ('{}')",
                    path.display(),
                    global_settings.state.folder.display()
                );
            }
            Self::Global { input } => {
                let fig = Figment::from(Serialized::defaults(&global_settings));
                let updates_display_directive = input.updates_global_display_directive()?;
                let updates_logfile_directive = input.updates_global_logfile_directive()?;

                if updates_logfile_directive {
                    if let Some(reason) = file_log_boot_disabled_reason() {
                        return Err(eyre!(
                            "Cannot change global.logfile_directive because this session was started with the logfile logger disabled ({reason})"
                        ));
                    }
                }

                *global_settings = input.merge_figment(fig)?.extract()?;
                if updates_display_directive {
                    crate::tracing::set_stderr_log_filter_override(None)?;
                }
                if updates_logfile_directive {
                    clear_file_log_filter_override_on_settings_change()?;
                }
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
    pub fn updates_global_display_directive(&self) -> Result<bool> {
        self.updates_global_path(&["global", "display_directive"], "global.display_directive")
    }

    pub fn updates_global_logfile_directive(&self) -> Result<bool> {
        self.updates_global_path(&["global", "logfile_directive"], "global.logfile_directive")
    }

    fn updates_global_path(&self, path: &[&str], label: &str) -> Result<bool> {
        let dotted_path = path.join(".");
        match self {
            SetArgs::Kv { pairs } => Ok(pairs.iter().any(|pair| pair.key == dotted_path)),
            SetArgs::String { string } => {
                let value: toml::Value = toml::from_str(string).with_context(|| {
                    format!("Failed parsing TOML payload while checking for {label}")
                })?;
                Ok(toml_contains_path(&value, path))
            }
            SetArgs::File { file } => {
                let ext = file
                    .extension()
                    .and_then(|s| s.to_str())
                    .unwrap_or("")
                    .to_ascii_lowercase();
                let contents = fs::read_to_string(file).with_context(|| {
                    format!("Failed reading settings file '{}'", file.display())
                })?;
                match ext.as_str() {
                    "toml" => {
                        let value: toml::Value = toml::from_str(&contents).with_context(|| {
                            format!(
                                "Failed parsing TOML settings file '{}' while checking for {label}",
                                file.display()
                            )
                        })?;
                        Ok(toml_contains_path(&value, path))
                    }
                    "json" => {
                        let value: JsonValue =
                            serde_json::from_str(&contents).with_context(|| {
                                format!(
                                "Failed parsing JSON settings file '{}' while checking for {label}",
                                file.display()
                            )
                            })?;
                        Ok(json_contains_path(&value, path))
                    }
                    _ => Ok(false),
                }
            }
            SetArgs::Stored | SetArgs::Defaults => Ok(false),
        }
    }

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

fn toml_contains_path(value: &toml::Value, path: &[&str]) -> bool {
    match path.split_first() {
        None => true,
        Some((head, tail)) => value
            .as_table()
            .and_then(|table| table.get(*head))
            .is_some_and(|next| toml_contains_path(next, tail)),
    }
}

fn json_contains_path(value: &JsonValue, path: &[&str]) -> bool {
    match path.split_first() {
        None => true,
        Some((head, tail)) => value
            .as_object()
            .and_then(|table| table.get(*head))
            .is_some_and(|next| json_contains_path(next, tail)),
    }
}

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
    use std::{str::FromStr, sync::Mutex};

    use clap::Parser;
    use figment::{providers::Serialized, Figment};
    use gammalooprs::{
        model::{ParameterNature, ParameterType, UFOSymbol},
        settings::RuntimeSettings,
        utils::{test_utils::load_generic_model, F},
    };
    use serde::{Deserialize, Serialize};
    use spenso::algebra::complex::Complex;

    use crate::{
        state::{ProcessRef, State},
        tracing::{get_stderr_log_filter, set_stderr_log_filter, set_stderr_log_filter_override},
        CLISettings, Repl,
    };

    use super::{
        super::Commands, model_value_format_hint, validate_model_parameter_type, KvPair, Set,
        SetArgs,
    };
    use super::{
        parse_model_parameter_value, MODEL_COMPLEX_VALUE_FORMAT_HINT, MODEL_REAL_VALUE_FORMAT_HINT,
    };

    #[test]
    fn serialize_complex() {
        let s = Complex::new(F(1.), F(-2.));
        let j = serde_json::to_string(&s).unwrap();
        assert_eq!(j, r#"{"re":1.0,"im":-2.0}"#.to_string());
    }

    #[test]
    fn parse_model_parameter_value_supports_real_literals() {
        let parsed = parse_model_parameter_value("-1.0e12").unwrap();
        assert_eq!(parsed, Complex::new(F(-1.0e12), F(0.0)));
    }

    #[test]
    fn parse_model_parameter_value_supports_i_and_j_suffixes() {
        let with_i = parse_model_parameter_value("-1.0e13-33.0e12i").unwrap();
        let with_j = parse_model_parameter_value("-1.0e13+33.0e12j").unwrap();
        let pure_imag = parse_model_parameter_value("-2.5j").unwrap();

        assert_eq!(with_i, Complex::new(F(-1.0e13), F(-33.0e12)));
        assert_eq!(with_j, Complex::new(F(-1.0e13), F(33.0e12)));
        assert_eq!(pure_imag, Complex::new(F(0.0), F(-2.5)));
    }

    #[test]
    fn parse_model_parameter_value_rejects_legacy_json_style() {
        let err = parse_model_parameter_value("{re:1.0,im:0.0}").unwrap_err();
        assert!(err.to_string().contains(MODEL_COMPLEX_VALUE_FORMAT_HINT));

        let err = parse_model_parameter_value("[1.0,0.0]").unwrap_err();
        assert!(err.to_string().contains(MODEL_COMPLEX_VALUE_FORMAT_HINT));
    }

    #[test]
    fn model_value_format_hint_depends_on_parameter_type() {
        assert_eq!(
            model_value_format_hint(Some(ParameterType::Real)),
            MODEL_REAL_VALUE_FORMAT_HINT
        );
        assert_eq!(
            model_value_format_hint(Some(ParameterType::Imaginary)),
            MODEL_COMPLEX_VALUE_FORMAT_HINT
        );
        assert_eq!(
            model_value_format_hint(None),
            MODEL_COMPLEX_VALUE_FORMAT_HINT
        );
    }

    #[test]
    fn validate_model_parameter_type_rejects_imaginary_part_for_real_parameters() {
        let err = validate_model_parameter_type(
            "alpha",
            ParameterType::Real,
            &Complex::new(F(1.0), F(2.0)),
        )
        .unwrap_err();
        assert!(err
            .to_string()
            .contains("cannot be assigned an imaginary component"));
        assert!(err.to_string().contains(MODEL_REAL_VALUE_FORMAT_HINT));
    }

    #[test]
    fn validate_model_parameter_type_accepts_complex_values_for_imaginary_parameters() {
        validate_model_parameter_type(
            "alpha",
            ParameterType::Imaginary,
            &Complex::new(F(1.0), F(2.0)),
        )
        .unwrap();
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

    #[test]
    fn set_model_accepts_python_style_complex_literals() {
        let mut state = State::new_test();
        state.model = load_generic_model("sm");
        state.model_parameters =
            gammalooprs::model::InputParamCard::default_from_model(&state.model);
        let parameter = state
            .model_parameters
            .keys()
            .next()
            .expect("expected at least one external model parameter")
            .to_string();

        Set::Model {
            pairs: vec![KvPair {
                key: parameter.clone(),
                value: "-1.0e13-33.0e12i".to_string(),
            }],
        }
        .run(
            &mut state,
            &mut CLISettings::default(),
            &mut RuntimeSettings::default(),
        )
        .unwrap();

        assert_eq!(
            state.model_parameters[&UFOSymbol::from(parameter.as_str())],
            Complex::new(F(-1.0e13), F(-33.0e12))
        );
    }

    #[test]
    fn set_model_rejects_imaginary_component_for_real_parameters() {
        let mut state = State::new_test();
        state.model = load_generic_model("sm");
        state.model_parameters =
            gammalooprs::model::InputParamCard::default_from_model(&state.model);
        let parameter = state
            .model
            .parameters
            .values()
            .find(|parameter| {
                parameter.nature == ParameterNature::External
                    && parameter.parameter_type == ParameterType::Real
            })
            .expect("expected at least one external real model parameter")
            .name
            .to_string();

        let err = Set::Model {
            pairs: vec![KvPair {
                key: parameter,
                value: "1.0+2.0i".to_string(),
            }],
        }
        .run(
            &mut state,
            &mut CLISettings::default(),
            &mut RuntimeSettings::default(),
        )
        .unwrap_err();

        assert!(err
            .to_string()
            .contains("cannot be assigned an imaginary component"));
    }

    #[test]
    fn updates_global_display_directive_detects_kv_string_and_file_inputs() {
        let tmp_name = format!(
            "gammaloop_set_display_directive_{}_{}.toml",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        );
        let file_path = std::env::temp_dir().join(tmp_name);
        std::fs::write(
            &file_path,
            "[global]\ndisplay_directive = \"warn\"\nlogfile_directive = \"off\"\n",
        )
        .unwrap();

        assert!(SetArgs::Kv {
            pairs: vec![KvPair {
                key: "global.display_directive".to_string(),
                value: "warn".to_string(),
            }]
        }
        .updates_global_display_directive()
        .unwrap());

        assert!(SetArgs::String {
            string: "[global]\ndisplay_directive = \"warn\"\n".to_string(),
        }
        .updates_global_display_directive()
        .unwrap());

        assert!(SetArgs::File {
            file: file_path.clone(),
        }
        .updates_global_display_directive()
        .unwrap());

        assert!(!SetArgs::Kv {
            pairs: vec![KvPair {
                key: "global.n_cores.generate".to_string(),
                value: "4".to_string(),
            }]
        }
        .updates_global_display_directive()
        .unwrap());

        let _ = std::fs::remove_file(file_path);
    }

    #[test]
    fn updates_global_logfile_directive_detects_kv_string_and_file_inputs() {
        let tmp_name = format!(
            "gammaloop_set_logfile_directive_{}_{}.toml",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        );
        let file_path = std::env::temp_dir().join(tmp_name);
        std::fs::write(
            &file_path,
            "[global]\ndisplay_directive = \"warn\"\nlogfile_directive = \"debug\"\n",
        )
        .unwrap();

        assert!(SetArgs::Kv {
            pairs: vec![KvPair {
                key: "global.logfile_directive".to_string(),
                value: "debug".to_string(),
            }]
        }
        .updates_global_logfile_directive()
        .unwrap());

        assert!(SetArgs::String {
            string: "[global]\nlogfile_directive = \"debug\"\n".to_string(),
        }
        .updates_global_logfile_directive()
        .unwrap());

        assert!(SetArgs::File {
            file: file_path.clone(),
        }
        .updates_global_logfile_directive()
        .unwrap());

        assert!(!SetArgs::Kv {
            pairs: vec![KvPair {
                key: "global.display_directive".to_string(),
                value: "warn".to_string(),
            }]
        }
        .updates_global_logfile_directive()
        .unwrap());

        let _ = std::fs::remove_file(file_path);
    }

    #[test]
    fn set_global_display_directive_clears_cli_stderr_override() {
        static TEST_MUTEX: Mutex<()> = Mutex::new(());
        let _guard = TEST_MUTEX.lock().unwrap();

        let mut state = State::new_test();
        let mut cli_settings = CLISettings::default();
        let mut runtime_settings = RuntimeSettings::default();

        set_stderr_log_filter("info").unwrap();
        set_stderr_log_filter_override(Some("gammaloop_api=debug,gammalooprs=debug".to_string()))
            .unwrap();

        Set::Global {
            input: SetArgs::Kv {
                pairs: vec![KvPair {
                    key: "global.display_directive".to_string(),
                    value: "warn".to_string(),
                }],
            },
        }
        .run(&mut state, &mut cli_settings, &mut runtime_settings)
        .unwrap();

        assert_eq!(cli_settings.global.display_directive, "warn");
        assert_eq!(get_stderr_log_filter(), "warn");
    }

    #[test]
    fn set_global_logfile_directive_errors_when_file_logger_was_boot_disabled() {
        static TEST_MUTEX: Mutex<()> = Mutex::new(());
        let _guard = TEST_MUTEX.lock().unwrap();

        let mut state = State::new_test();
        let mut cli_settings = CLISettings::default();
        let mut runtime_settings = RuntimeSettings::default();

        crate::tracing::set_file_log_filter("off").unwrap();
        crate::tracing::set_file_log_filter_override(Some(
            "gammaloop_api=off,gammalooprs=off".to_string(),
        ))
        .unwrap();
        crate::tracing::configure_file_log_boot_mode(true, Some("--logfile-level off")).unwrap();

        let err = Set::Global {
            input: SetArgs::Kv {
                pairs: vec![KvPair {
                    key: "global.logfile_directive".to_string(),
                    value: "debug".to_string(),
                }],
            },
        }
        .run(&mut state, &mut cli_settings, &mut runtime_settings)
        .unwrap_err();

        assert!(format!("{err:?}").contains("logfile logger disabled"));
        crate::tracing::configure_file_log_boot_mode(false, None).unwrap();
        crate::tracing::set_file_log_filter_override(None).unwrap();
    }
}
