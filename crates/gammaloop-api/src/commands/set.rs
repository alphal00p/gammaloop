use std::{collections::BTreeMap, fs, path::PathBuf, str::FromStr};

use clap::Subcommand;
use color_eyre::{Result, Section};
use eyre::{eyre, Context, Report};
use figment::{
    providers::{Format, Serialized},
    Figment,
};
use gammalooprs::{
    model::{ParameterNature, UFOSymbol},
    processes::ProcessCollection,
    settings::RuntimeSettings,
    utils::F,
};
use schemars::JsonSchema;
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use serde_json::Value as JsonValue;
use spenso::algebra::complex::Complex;
use tracing::warn;

use crate::{
    commands::generate::ProcessArgs,
    commands::process_settings::{
        observable_template, parse_quantity_kind, parse_selector_kind, quantity_template,
        selector_template,
    },
    commands::Commands,
    model_parameters::{
        external_model_parameter_type, parse_model_parameter_value, validate_model_parameter_type,
    },
    state::{State, SyncSettings},
    tracing::{clear_file_log_filter_override_on_settings_change, file_log_boot_disabled_reason},
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
        #[command(flatten)]
        target: ProcessArgs,
        /// Any number of PARAM=COMPLEX pairs, or 'defaults' for process-targeted resets
        #[arg(
            value_name = "defaults|PARAM=COMPLEX",
            num_args = 1..,
            value_parser = ModelSetValue::from_str
        )]
        values: Vec<ModelSetValue>,
    },

    /// Set settings for a PROCESS
    Process {
        #[command(subcommand)]
        input: ProcessSetArgs,
        #[command(flatten)]
        process: ProcessArgs,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub enum ModelSetValue {
    Defaults,
    Assignment(KvPair),
}

impl FromStr for ModelSetValue {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        if value.trim() == "defaults" {
            return Ok(Self::Defaults);
        }
        KvPair::from_str(value).map(Self::Assignment)
    }
}

fn resolve_generated_model_targets(
    state: &State,
    target: &ProcessArgs,
) -> Result<Option<Vec<(usize, String)>>> {
    match (&target.process, &target.integrand_name) {
        (None, None) => Ok(None),
        (None, Some(integrand_name)) => Ok(Some(vec![
            state.find_generated_integrand_ref_by_name(integrand_name)?
        ])),
        (Some(process), None) => {
            let process_id = state.resolve_process_ref(Some(process))?;
            let process = &state.process_list.processes[process_id];
            let generated = match &process.collection {
                ProcessCollection::Amplitudes(amplitudes) => amplitudes
                    .iter()
                    .filter(|(_, amplitude)| amplitude.integrand.is_some())
                    .map(|(name, _)| (process_id, name.clone()))
                    .collect::<Vec<_>>(),
                ProcessCollection::CrossSections(cross_sections) => cross_sections
                    .iter()
                    .filter(|(_, cross_section)| cross_section.integrand.is_some())
                    .map(|(name, _)| (process_id, name.clone()))
                    .collect::<Vec<_>>(),
            };
            if generated.is_empty() {
                return Err(eyre!(
                    "Process '{}' has no generated integrands. Per-integrand model parameters are only supported for generated integrands.",
                    process.definition.folder_name
                ));
            }
            Ok(Some(generated))
        }
        (Some(process), Some(integrand_name)) => {
            let process_id = state.resolve_process_ref(Some(process))?;
            let canonical_name = state.process_list.processes[process_id]
                .collection
                .find_integrand(Some(integrand_name.clone()))?;
            state
                .process_list
                .get_integrand(process_id, &canonical_name)
                .and_then(|resolved| resolved.require_generated().map(|_| resolved))
                .with_context(|| {
                    "Per-integrand model parameters are only supported for generated integrands"
                })?;
            Ok(Some(vec![(process_id, canonical_name)]))
        }
    }
}

#[allow(clippy::type_complexity)]
fn parse_model_assignments(
    state: &State,
    values: &[ModelSetValue],
) -> Result<Option<Vec<(String, Complex<F<f64>>)>>> {
    if values
        .iter()
        .any(|value| matches!(value, ModelSetValue::Defaults))
    {
        if values.len() != 1 {
            return Err(eyre!(
                "Model parameter updates cannot mix 'defaults' with explicit PARAM=VALUE assignments"
            ));
        }
        return Ok(None);
    }

    let assignments = values
        .iter()
        .map(|value| match value {
            ModelSetValue::Assignment(KvPair { key, value }) => {
                let parameter = state
                    .model
                    .get_parameter_opt(key)
                    .filter(|parameter| parameter.nature == ParameterNature::External);
                if parameter.is_none() {
                    let possibilities = state
                        .model_parameters
                        .keys()
                        .map(|s| s.to_string())
                        .collect::<Vec<_>>();
                    return Err(eyre!("No model parameter named '{key}'")).with_context(|| {
                        format!(
                            "Possible model parameters are: {}",
                            possibilities.join(", ")
                        )
                    });
                }
                let parameter_type = parameter.unwrap().parameter_type.clone();
                let value = parse_model_parameter_value(value).with_context(|| {
                    format!("While parsing model parameter value {value} for key '{key}'")
                })?;
                validate_model_parameter_type(key, parameter_type, &value)?;
                Ok((key.clone(), value))
            }
            ModelSetValue::Defaults => unreachable!("defaults handled above"),
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(Some(assignments))
}

fn apply_model_assignments(
    settings: &mut RuntimeSettings,
    default_model_parameters: &gammalooprs::model::InputParamCard<F<f64>>,
    assignments: &[(String, Complex<F<f64>>)],
) -> Result<()> {
    for (parameter_name, value) in assignments {
        let symbol = UFOSymbol::from(parameter_name.as_str());
        let default_value = default_model_parameters.get(&symbol).ok_or_else(|| {
            eyre!(
                "Model parameter '{parameter_name}' cannot be overridden because it is not present in the shared top-level model_parameters.json"
            )
        })?;

        if default_value == value {
            settings.model.external_parameters.remove(parameter_name);
        } else {
            settings
                .model
                .external_parameters
                .insert(parameter_name.clone(), (value.re, value.im));
        }
    }

    Ok(())
}

impl Set {
    pub fn run(
        &self,
        state: &mut State,
        global_settings: &mut CLISettings,
        default_runtime_settings: &mut RuntimeSettings,
    ) -> Result<()> {
        let runtime_model_validation = RuntimeModelValidationContext::from_state(state);
        match self {
            Set::Model { target, values } => {
                let targets = resolve_generated_model_targets(state, target)?;
                let assignments = parse_model_assignments(state, values)?;

                if let Some(targets) = targets {
                    let default_model_parameters = state.model_parameters.clone();
                    if let Some(assignments) = assignments.as_ref() {
                        for (process_id, integrand_name) in targets {
                            let settings = state
                                .process_list
                                .get_integrand_mut(process_id, &integrand_name)?
                                .get_mut_settings();
                            apply_model_assignments(
                                settings,
                                &default_model_parameters,
                                assignments,
                            )?;
                        }
                    } else {
                        for (process_id, integrand_name) in targets {
                            let settings = state
                                .process_list
                                .get_integrand_mut(process_id, &integrand_name)?
                                .get_mut_settings();
                            settings.model = default_runtime_settings.model.clone();
                        }
                    }
                } else {
                    let Some(assignments) = assignments else {
                        return Err(eyre!(
                            "The 'defaults' shortcut for set model requires a process or integrand target"
                        ));
                    };

                    for (parameter_name, value) in assignments {
                        if let Some(parameter) = state
                            .model_parameters
                            .get_mut(&UFOSymbol::from(parameter_name.as_str()))
                        {
                            *parameter = value;
                            continue;
                        }
                        let possibilities = state
                            .model_parameters
                            .keys()
                            .map(|s| s.to_string())
                            .collect::<Vec<_>>();
                        return Err(eyre!("No model parameter named '{parameter_name}'"))
                            .with_context(|| {
                                format!(
                                    "Possible model parameters are: {}",
                                    possibilities.join(", ")
                                )
                            });
                    }
                    state.model_parameters.apply_to_model(&mut state.model)?;
                }
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
                let merged = merge_runtime_settings_input(default_runtime_settings, input)?;
                runtime_model_validation.validate(&merged)?;
                *default_runtime_settings = merged;
            }
            Self::Process { input, process } => {
                let process_id = state.resolve_process_ref(process.process.as_ref())?;
                let apply_runtime_settings = |settings: &mut RuntimeSettings| -> Result<()> {
                    apply_process_set_args(
                        input,
                        settings,
                        default_runtime_settings,
                        Some(&runtime_model_validation),
                    )
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

#[derive(Subcommand, Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub enum ProcessSetArgs {
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

    /// Add a named quantity/observable/selector
    Add {
        #[command(subcommand)]
        target: ProcessAddTarget,
    },

    /// Update a named quantity/observable/selector
    Update {
        #[command(subcommand)]
        target: ProcessUpdateTarget,
    },

    /// Remove a named quantity/observable/selector
    Remove {
        #[command(subcommand)]
        target: ProcessRemoveTarget,
    },
}

#[derive(Subcommand, Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub enum ProcessAddTarget {
    Quantity {
        name: String,
        #[arg(value_parser = parse_quantity_kind)]
        kind: String,
        #[arg(value_name = "KEY=VALUE", num_args = 0.., value_parser = KvPair::from_str)]
        pairs: Vec<KvPair>,
    },
    Observable {
        name: String,
        #[arg(value_name = "KEY=VALUE", num_args = 0.., value_parser = KvPair::from_str)]
        pairs: Vec<KvPair>,
    },
    Selector {
        name: String,
        #[arg(value_parser = parse_selector_kind)]
        kind: String,
        #[arg(value_name = "KEY=VALUE", num_args = 0.., value_parser = KvPair::from_str)]
        pairs: Vec<KvPair>,
    },
}

#[derive(Subcommand, Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub enum ProcessUpdateTarget {
    Quantity {
        name: String,
        #[arg(value_name = "KEY=VALUE", num_args = 1.., value_parser = KvPair::from_str)]
        pairs: Vec<KvPair>,
    },
    Observable {
        name: String,
        #[arg(value_name = "KEY=VALUE", num_args = 1.., value_parser = KvPair::from_str)]
        pairs: Vec<KvPair>,
    },
    Selector {
        name: String,
        #[arg(value_name = "KEY=VALUE", num_args = 1.., value_parser = KvPair::from_str)]
        pairs: Vec<KvPair>,
    },
}

#[derive(Subcommand, Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub enum ProcessRemoveTarget {
    Quantity { name: String },
    Observable { name: String },
    Selector { name: String },
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

impl ProcessSetArgs {
    fn as_set_args(&self) -> Option<SetArgs> {
        match self {
            ProcessSetArgs::File { file } => Some(SetArgs::File { file: file.clone() }),
            ProcessSetArgs::String { string } => Some(SetArgs::String {
                string: string.clone(),
            }),
            ProcessSetArgs::Kv { pairs } => Some(SetArgs::Kv {
                pairs: pairs.clone(),
            }),
            ProcessSetArgs::Defaults => Some(SetArgs::Defaults),
            ProcessSetArgs::Stored => Some(SetArgs::Stored),
            ProcessSetArgs::Add { .. }
            | ProcessSetArgs::Update { .. }
            | ProcessSetArgs::Remove { .. } => None,
        }
    }
}

fn apply_process_set_args(
    input: &ProcessSetArgs,
    settings: &mut RuntimeSettings,
    default_runtime_settings: &RuntimeSettings,
    validation_context: Option<&RuntimeModelValidationContext>,
) -> Result<()> {
    match input {
        ProcessSetArgs::Add { target } => apply_process_add_target(settings, target),
        ProcessSetArgs::Update { target } => apply_process_update_target(settings, target),
        ProcessSetArgs::Remove { target } => apply_process_remove_target(settings, target),
        ProcessSetArgs::Defaults => {
            *settings = default_runtime_settings.clone();
            Ok(())
        }
        _ => {
            let input = input
                .as_set_args()
                .expect("non-mutation process set args should convert to SetArgs");
            *settings = merge_runtime_settings_input(settings, &input)?;
            if let Some(validation_context) = validation_context {
                validation_context.validate(settings)?;
            }
            Ok(())
        }
    }
}

fn apply_process_add_target(
    settings: &mut RuntimeSettings,
    target: &ProcessAddTarget,
) -> Result<()> {
    match target {
        ProcessAddTarget::Quantity { name, kind, pairs } => {
            let entry = merge_named_settings(
                &quantity_template(kind).ok_or_else(|| eyre!("Unknown quantity kind '{kind}'"))?,
                pairs,
                &format!("new quantity '{name}'"),
            )?;
            add_named_map_entry(&mut settings.quantities, name, entry, "quantity")
        }
        ProcessAddTarget::Observable { name, pairs } => {
            let entry = merge_named_settings(
                &observable_template(),
                pairs,
                &format!("new observable '{name}'"),
            )?;
            add_named_map_entry(&mut settings.observables, name, entry, "observable")
        }
        ProcessAddTarget::Selector { name, kind, pairs } => {
            let entry = merge_named_settings(
                &selector_template(kind).ok_or_else(|| eyre!("Unknown selector kind '{kind}'"))?,
                pairs,
                &format!("new selector '{name}'"),
            )?;
            add_named_map_entry(&mut settings.selectors, name, entry, "selector")
        }
    }
}

fn apply_process_update_target(
    settings: &mut RuntimeSettings,
    target: &ProcessUpdateTarget,
) -> Result<()> {
    match target {
        ProcessUpdateTarget::Quantity { name, pairs } => {
            update_named_map_entry(&mut settings.quantities, name, pairs, "quantity")
        }
        ProcessUpdateTarget::Observable { name, pairs } => {
            update_named_map_entry(&mut settings.observables, name, pairs, "observable")
        }
        ProcessUpdateTarget::Selector { name, pairs } => {
            update_named_map_entry(&mut settings.selectors, name, pairs, "selector")
        }
    }
}

fn apply_process_remove_target(
    settings: &mut RuntimeSettings,
    target: &ProcessRemoveTarget,
) -> Result<()> {
    match target {
        ProcessRemoveTarget::Quantity { name } => {
            remove_named_map_entry(&mut settings.quantities, name, "quantity")
        }
        ProcessRemoveTarget::Observable { name } => {
            remove_named_map_entry(&mut settings.observables, name, "observable")
        }
        ProcessRemoveTarget::Selector { name } => {
            remove_named_map_entry(&mut settings.selectors, name, "selector")
        }
    }
}

fn merge_named_settings<T>(base: &T, pairs: &[KvPair], label: &str) -> Result<T>
where
    T: Serialize + DeserializeOwned,
{
    let fig = Figment::from(Serialized::defaults(base));
    SetArgs::Kv {
        pairs: pairs.to_vec(),
    }
    .merge_figment(fig)?
    .extract()
    .with_context(|| format!("While building {label}"))
}

fn add_named_map_entry<T>(
    map: &mut BTreeMap<String, T>,
    name: &str,
    entry: T,
    label: &str,
) -> Result<()> {
    if map.contains_key(name) {
        return Err(eyre!("A {label} named '{name}' already exists"));
    }
    map.insert(name.to_string(), entry);
    Ok(())
}

fn update_named_map_entry<T>(
    map: &mut BTreeMap<String, T>,
    name: &str,
    pairs: &[KvPair],
    label: &str,
) -> Result<()>
where
    T: Clone + Serialize + DeserializeOwned,
{
    let current = map
        .get(name)
        .cloned()
        .ok_or_else(|| eyre!("No {label} named '{name}'"))?;
    let updated = merge_named_settings(&current, pairs, &format!("{label} '{name}'"))?;
    map.insert(name.to_string(), updated);
    Ok(())
}

fn remove_named_map_entry<T>(map: &mut BTreeMap<String, T>, name: &str, label: &str) -> Result<()> {
    if map.remove(name).is_none() {
        return Err(eyre!("No {label} named '{name}'"));
    }
    Ok(())
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

fn json_number_to_f64(value: &serde_json::Number) -> Result<f64> {
    value
        .as_f64()
        .ok_or_else(|| eyre!("Could not represent JSON number '{value}' as f64"))
}

fn parse_runtime_model_setting_value(value: &J) -> Result<(F<f64>, F<f64>)> {
    match value {
        J::Number(number) => Ok((F(json_number_to_f64(number)?), F(0.0))),
        J::String(string) => {
            let parsed = parse_model_parameter_value(string)?;
            Ok((parsed.re, parsed.im))
        }
        J::Array(entries) if entries.len() == 2 => {
            let parse_component = |entry: &J| -> Result<F<f64>> {
                match entry {
                    J::Number(number) => Ok(F(json_number_to_f64(number)?)),
                    other => Err(eyre!(
                        "Runtime model setting components must be numeric, found {other:?}"
                    )),
                }
            };
            Ok((parse_component(&entries[0])?, parse_component(&entries[1])?))
        }
        other => Err(eyre!(
            "Runtime model settings must be specified as a number, a complex string, or a two-element numeric array, found {other:?}"
        )),
    }
}

#[allow(clippy::type_complexity)]
fn extract_runtime_model_updates_from_json(
    value: &mut J,
) -> Result<BTreeMap<String, (F<f64>, F<f64>)>> {
    let Some(object) = value.as_object_mut() else {
        return Ok(BTreeMap::new());
    };

    let Some(model_value) = object.remove("model") else {
        return Ok(BTreeMap::new());
    };

    let J::Object(model_object) = model_value else {
        return Err(eyre!(
            "The runtime settings 'model' block must be a table/object"
        ));
    };

    model_object
        .into_iter()
        .map(|(parameter_name, value)| {
            parse_runtime_model_setting_value(&value)
                .map(|parsed| (parameter_name.clone(), parsed))
                .with_context(|| format!("While parsing runtime model setting '{parameter_name}'"))
        })
        .collect()
}

#[allow(clippy::type_complexity)]
fn extract_runtime_model_updates_from_toml(
    value: &mut toml::Value,
) -> Result<BTreeMap<String, (F<f64>, F<f64>)>> {
    let Some(table) = value.as_table_mut() else {
        return Ok(BTreeMap::new());
    };

    let Some(model_value) = table.remove("model") else {
        return Ok(BTreeMap::new());
    };

    let mut model_json = serde_json::json!({ "model": serde_json::to_value(model_value)? });
    extract_runtime_model_updates_from_json(&mut model_json)
}

#[derive(Clone)]
struct RuntimeModelValidationContext {
    model: gammalooprs::model::Model,
    model_parameters: gammalooprs::model::InputParamCard<F<f64>>,
}

impl RuntimeModelValidationContext {
    fn from_state(state: &State) -> Self {
        Self {
            model: state.model.clone(),
            model_parameters: state.model_parameters.clone(),
        }
    }

    fn validate(&self, settings: &RuntimeSettings) -> Result<()> {
        let mut possibilities = self
            .model_parameters
            .keys()
            .map(|symbol| symbol.to_string())
            .collect::<Vec<_>>();
        possibilities.sort();

        for (parameter_name, value) in &settings.model.external_parameters {
            let parameter_type = external_model_parameter_type(&self.model, parameter_name)
                .ok_or_else(|| eyre!("No model parameter named '{parameter_name}'"))
                .with_note(|| {
                    format!(
                        "Possible model parameters are: {}",
                        possibilities.join(", ")
                    )
                })?;

            let symbol = UFOSymbol::from(parameter_name.as_str());
            if !self.model_parameters.contains_key(&symbol) {
                return Err(eyre!(
                    "Model parameter '{parameter_name}' cannot be overridden because it is not present in the shared top-level model_parameters.json"
                ))
                .with_note(|| format!("Possible model parameters are: {}", possibilities.join(", ")));
            }

            validate_model_parameter_type(
                parameter_name,
                parameter_type,
                &Complex::new(value.0, value.1),
            )?;
        }

        Ok(())
    }
}

fn merge_runtime_settings_input(
    settings: &RuntimeSettings,
    input: &SetArgs,
) -> Result<RuntimeSettings> {
    let mut settings_without_model = settings.clone();
    let preserved_model = settings_without_model.model.clone();
    settings_without_model.model = Default::default();

    let mut merged_settings: Option<RuntimeSettings> = None;
    let mut model_updates = BTreeMap::new();

    match input {
        SetArgs::File { file } => {
            let ext = file
                .extension()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_ascii_lowercase();
            match ext.as_str() {
                "toml" => {
                    let raw = fs::read_to_string(file).with_context(|| {
                        format!("Trying to read runtime settings file {}", file.display())
                    })?;
                    let mut value = toml::from_str::<toml::Value>(&raw).with_context(|| {
                        format!(
                            "Trying to parse TOML runtime settings file {}",
                            file.display()
                        )
                    })?;
                    model_updates = extract_runtime_model_updates_from_toml(&mut value)?;
                    if value.as_table().is_some_and(|table| !table.is_empty()) {
                        let fig = Figment::from(Serialized::defaults(&settings_without_model));
                        merged_settings = Some(
                            fig.merge(figment::providers::Toml::string(&toml::to_string(&value)?))
                                .extract()?,
                        );
                    }
                }
                "json" => {
                    let raw = fs::read_to_string(file).with_context(|| {
                        format!("Trying to read runtime settings file {}", file.display())
                    })?;
                    let mut value = serde_json::from_str::<J>(&raw).with_context(|| {
                        format!(
                            "Trying to parse JSON runtime settings file {}",
                            file.display()
                        )
                    })?;
                    model_updates = extract_runtime_model_updates_from_json(&mut value)?;
                    if value.as_object().is_some_and(|object| !object.is_empty()) {
                        let fig = Figment::from(Serialized::defaults(&settings_without_model));
                        merged_settings = Some(
                            fig.merge(figment::providers::Json::string(&serde_json::to_string(
                                &value,
                            )?))
                            .extract()?,
                        );
                    }
                }
                _ => {
                    return Err(color_eyre::eyre::eyre!(
                        "Unsupported settings file extension: {}",
                        ext
                    ));
                }
            }
        }
        SetArgs::String { string } => {
            let mut value = toml::from_str::<toml::Value>(string)
                .with_context(|| "Trying to parse runtime settings TOML string")?;
            model_updates = extract_runtime_model_updates_from_toml(&mut value)?;
            if value.as_table().is_some_and(|table| !table.is_empty()) {
                let fig = Figment::from(Serialized::defaults(&settings_without_model));
                merged_settings = Some(
                    fig.merge(figment::providers::Toml::string(&toml::to_string(&value)?))
                        .extract()?,
                );
            }
        }
        SetArgs::Kv { pairs } => {
            let mut fig = Figment::from(Serialized::defaults(&settings_without_model));
            let mut saw_non_model_pair = false;
            for KvPair { key, value } in pairs {
                if let Some(parameter_name) = key.strip_prefix("model.") {
                    model_updates.insert(
                        parameter_name.to_string(),
                        parse_runtime_model_setting_value(&infer_cli_value(value)?).with_context(
                            || format!("While parsing runtime model setting '{parameter_name}'"),
                        )?,
                    );
                } else {
                    saw_non_model_pair = true;
                    fig = fig.adjoin((key.clone(), infer_cli_value(value)?));
                }
            }

            if saw_non_model_pair {
                merged_settings = Some(fig.extract()?);
            }
        }
        SetArgs::Stored => {}
        SetArgs::Defaults => unreachable!("defaults should be handled by the caller"),
    }

    let mut merged_settings = merged_settings.unwrap_or_else(|| settings.clone());
    merged_settings.model = preserved_model;
    for (parameter_name, value) in model_updates {
        merged_settings
            .model
            .external_parameters
            .insert(parameter_name, value);
    }

    Ok(merged_settings)
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
        Y::Sequence(values) => values
            .into_iter()
            .map(yaml_to_json)
            .collect::<Result<Vec<_>>>()
            .map(J::Array),
        Y::Mapping(entries) => {
            let mut object = serde_json::Map::with_capacity(entries.len());
            for (key, value) in entries {
                let key = match key {
                    Y::String(s) => s,
                    Y::Bool(b) => b.to_string(),
                    Y::Number(n) => n.to_string(),
                    Y::Null => "null".to_string(),
                    other => {
                        return Err(eyre!(
                            "Unsupported YAML mapping key in CLI key-value: {:?}",
                            other
                        ));
                    }
                };
                object.insert(key, yaml_to_json(value)?);
            }
            Ok(J::Object(object))
        }
        other => Err(eyre!(
            "Unsupported YAML value in CLI key-value: {:?}",
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
        model::{ParameterNature, ParameterType, UFOSymbol},
        observables::{FilterQuantity, QuantitySettings, SelectorDefinitionSettings},
        settings::RuntimeSettings,
        utils::{load_generic_model, F},
    };
    use serde::{Deserialize, Serialize};
    use spenso::algebra::complex::Complex;

    use crate::{
        model_parameters::{
            model_value_format_hint, parse_model_parameter_value, MODEL_COMPLEX_VALUE_FORMAT_HINT,
            MODEL_REAL_VALUE_FORMAT_HINT,
        },
        state::{ProcessRef, State},
        tracing::{get_stderr_log_filter, set_stderr_log_filter, set_stderr_log_filter_override},
        CLISettings, Repl,
    };

    use super::{
        super::Commands, apply_process_set_args, validate_model_parameter_type, KvPair,
        ModelSetValue, ProcessAddTarget, ProcessArgs, ProcessRemoveTarget, ProcessSetArgs,
        ProcessUpdateTarget, Set, SetArgs,
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
    fn merge_runtime_settings_input_preserves_existing_model_overrides() {
        let mut settings = RuntimeSettings::default();
        settings
            .model
            .external_parameters
            .insert("mass_scalar_2".to_string(), (F(2.0), F(0.0)));

        let merged = super::merge_runtime_settings_input(
            &settings,
            &SetArgs::Kv {
                pairs: vec![KvPair {
                    key: "integrator.n_start".to_string(),
                    value: "123".to_string(),
                }],
            },
        )
        .unwrap();

        assert_eq!(merged.integrator.n_start, 123);
        assert_eq!(
            merged.model.external_parameters.get("mass_scalar_2"),
            Some(&(F(2.0), F(0.0)))
        );
    }

    #[test]
    fn merge_runtime_settings_input_accepts_model_block_updates() {
        let merged = super::merge_runtime_settings_input(
            &RuntimeSettings::default(),
            &SetArgs::String {
                string: "[model]\nmass_scalar_2 = [1.0, 0.0]\n".to_string(),
            },
        )
        .unwrap();

        assert_eq!(
            merged.model.external_parameters.get("mass_scalar_2"),
            Some(&(F(1.0), F(0.0)))
        );
    }

    #[test]
    fn process_settings_reject_non_overridable_model_parameters() {
        let mut state = State::new_test();
        state.model = load_generic_model("sm");
        state.model_parameters =
            gammalooprs::model::InputParamCard::default_from_model(&state.model);
        let removable_parameter = state
            .model
            .parameters
            .values()
            .find(|parameter| {
                parameter.nature == ParameterNature::External
                    && state.model_parameters.contains_key(&parameter.name)
            })
            .expect("SM model must contain at least one overridable external parameter")
            .name
            .to_string();
        state
            .model_parameters
            .remove(&UFOSymbol::from(removable_parameter.as_str()));
        let validation_context = super::RuntimeModelValidationContext::from_state(&state);
        let mut settings = RuntimeSettings::default();

        let err = super::apply_process_set_args(
            &ProcessSetArgs::Kv {
                pairs: vec![KvPair {
                    key: format!("model.{removable_parameter}"),
                    value: "1.0".to_string(),
                }],
            },
            &mut settings,
            &RuntimeSettings::default(),
            Some(&validation_context),
        )
        .unwrap_err();

        assert!(err
            .to_string()
            .contains("cannot be overridden because it is not present"));
    }

    #[test]
    fn process_settings_reject_invalid_model_value_types() {
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
        let validation_context = super::RuntimeModelValidationContext::from_state(&state);
        let mut settings = RuntimeSettings::default();

        let err = super::apply_process_set_args(
            &ProcessSetArgs::Kv {
                pairs: vec![KvPair {
                    key: format!("model.{parameter}"),
                    value: "1.0+2.0i".to_string(),
                }],
            },
            &mut settings,
            &RuntimeSettings::default(),
            Some(&validation_context),
        )
        .unwrap_err();

        assert!(err
            .to_string()
            .contains("cannot be assigned an imaginary component"));
    }

    #[test]
    fn parse_set_process_defaults() {
        let cmd = Set::from_str("set process -p epem_a_tth -i LO defaults").unwrap();

        match cmd {
            Set::Process { input, process } => {
                assert_eq!(input, ProcessSetArgs::Defaults);
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
    fn parse_set_model_targeted_defaults() {
        let cmd = Set::from_str("set model -p epem_a_tth -i LO defaults").unwrap();

        match cmd {
            Set::Model { target, values } => {
                assert_eq!(
                    target.process,
                    Some(ProcessRef::Unqualified("epem_a_tth".to_string()))
                );
                assert_eq!(target.integrand_name, Some("LO".to_string()));
                assert_eq!(values, vec![ModelSetValue::Defaults]);
            }
            other => panic!("Expected set model command, got {other:?}"),
        }
    }

    #[test]
    fn parse_set_model_targeted_assignment() {
        let cmd = Set::from_str("set model -i LO alpha=1.0+2.0i").unwrap();

        match cmd {
            Set::Model { target, values } => {
                assert_eq!(target.process, None);
                assert_eq!(target.integrand_name, Some("LO".to_string()));
                assert_eq!(
                    values,
                    vec![ModelSetValue::Assignment(KvPair {
                        key: "alpha".to_string(),
                        value: "1.0+2.0i".to_string(),
                    })]
                );
            }
            other => panic!("Expected set model command, got {other:?}"),
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
                    ProcessSetArgs::String {
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
                    ProcessSetArgs::String {
                        string: multiline.to_string()
                    }
                );
            }
            other => panic!("Expected set process command, got {other:?}"),
        }
    }

    #[test]
    fn parse_set_process_add_quantity() {
        let cmd = Set::from_str(
            "set process -p epem_a_tth -i LO add quantity top_pt particle_scalar quantity=PT",
        )
        .unwrap();

        match cmd {
            Set::Process { input, process } => {
                assert_eq!(
                    input,
                    ProcessSetArgs::Add {
                        target: ProcessAddTarget::Quantity {
                            name: "top_pt".to_string(),
                            kind: "particle_scalar".to_string(),
                            pairs: vec![KvPair {
                                key: "quantity".to_string(),
                                value: "PT".to_string(),
                            }],
                        }
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
    fn parse_set_process_update_observable() {
        let cmd = Set::from_str(
            "set process -p epem_a_tth -i LO update observable top_pt_hist n_bins=100",
        )
        .unwrap();

        match cmd {
            Set::Process { input, .. } => {
                assert_eq!(
                    input,
                    ProcessSetArgs::Update {
                        target: ProcessUpdateTarget::Observable {
                            name: "top_pt_hist".to_string(),
                            pairs: vec![KvPair {
                                key: "n_bins".to_string(),
                                value: "100".to_string(),
                            }],
                        }
                    }
                );
            }
            other => panic!("Expected set process command, got {other:?}"),
        }
    }

    #[test]
    fn parse_set_process_remove_selector() {
        let cmd = Set::from_str("set process -p epem_a_tth -i LO remove selector top_cut").unwrap();

        match cmd {
            Set::Process { input, .. } => {
                assert_eq!(
                    input,
                    ProcessSetArgs::Remove {
                        target: ProcessRemoveTarget::Selector {
                            name: "top_cut".to_string(),
                        }
                    }
                );
            }
            other => panic!("Expected set process command, got {other:?}"),
        }
    }

    #[test]
    fn process_named_settings_mutations_update_runtime_maps() {
        let defaults = RuntimeSettings::default();
        let mut settings = RuntimeSettings::default();

        apply_process_set_args(
            &ProcessSetArgs::Add {
                target: ProcessAddTarget::Quantity {
                    name: "top_pt".to_string(),
                    kind: "particle_scalar".to_string(),
                    pairs: vec![
                        KvPair {
                            key: "pdgs".to_string(),
                            value: "[6,-6]".to_string(),
                        },
                        KvPair {
                            key: "quantity".to_string(),
                            value: "PT".to_string(),
                        },
                    ],
                },
            },
            &mut settings,
            &defaults,
            None,
        )
        .unwrap();

        let quantity = settings
            .quantities
            .get("top_pt")
            .expect("quantity should have been inserted");
        match quantity {
            QuantitySettings::ParticleScalar(particle) => {
                assert_eq!(particle.pdgs, vec![6, -6]);
                assert_eq!(particle.quantity, FilterQuantity::PT);
            }
            other => panic!("Expected particle-scalar quantity, got {other:?}"),
        }

        apply_process_set_args(
            &ProcessSetArgs::Add {
                target: ProcessAddTarget::Observable {
                    name: "top_pt_hist".to_string(),
                    pairs: vec![
                        KvPair {
                            key: "quantity".to_string(),
                            value: "top_pt".to_string(),
                        },
                        KvPair {
                            key: "x_max".to_string(),
                            value: "500.0".to_string(),
                        },
                        KvPair {
                            key: "n_bins".to_string(),
                            value: "50".to_string(),
                        },
                    ],
                },
            },
            &mut settings,
            &defaults,
            None,
        )
        .unwrap();
        assert_eq!(
            settings
                .observables
                .get("top_pt_hist")
                .expect("observable should have been inserted")
                .histogram
                .n_bins,
            50
        );

        apply_process_set_args(
            &ProcessSetArgs::String {
                string: r#"
[selectors.top_pt_cut]
quantity = "top_pt"
selector = "value_range"
entry_selection = "leading_only"
min = 10.0
max = 500.0
"#
                .to_string(),
            },
            &mut settings,
            &defaults,
            None,
        )
        .unwrap();
        let selector = settings
            .selectors
            .get("top_pt_cut")
            .expect("selector should have been inserted");
        assert_eq!(selector.quantity, "top_pt");
        match &selector.selector {
            SelectorDefinitionSettings::ValueRange(selector) => {
                assert_eq!(selector.min, 10.0);
                assert_eq!(selector.max, Some(500.0));
            }
            other => panic!("Expected value-range selector, got {other:?}"),
        }

        apply_process_set_args(
            &ProcessSetArgs::Update {
                target: ProcessUpdateTarget::Observable {
                    name: "top_pt_hist".to_string(),
                    pairs: vec![KvPair {
                        key: "n_bins".to_string(),
                        value: "80".to_string(),
                    }],
                },
            },
            &mut settings,
            &defaults,
            None,
        )
        .unwrap();
        assert_eq!(
            settings
                .observables
                .get("top_pt_hist")
                .expect("observable should still exist")
                .histogram
                .n_bins,
            80
        );

        apply_process_set_args(
            &ProcessSetArgs::Update {
                target: ProcessUpdateTarget::Selector {
                    name: "top_pt_cut".to_string(),
                    pairs: vec![KvPair {
                        key: "max".to_string(),
                        value: "250.0".to_string(),
                    }],
                },
            },
            &mut settings,
            &defaults,
            None,
        )
        .unwrap();
        match &settings
            .selectors
            .get("top_pt_cut")
            .expect("selector should still exist")
            .selector
        {
            SelectorDefinitionSettings::ValueRange(selector) => {
                assert_eq!(selector.max, Some(250.0));
            }
            other => panic!("Expected value-range selector, got {other:?}"),
        }

        apply_process_set_args(
            &ProcessSetArgs::Remove {
                target: ProcessRemoveTarget::Quantity {
                    name: "top_pt".to_string(),
                },
            },
            &mut settings,
            &defaults,
            None,
        )
        .unwrap();
        assert!(!settings.quantities.contains_key("top_pt"));
    }

    #[test]
    fn infer_cli_value_accepts_yaml_sequences() {
        let parsed = super::infer_cli_value("[alphaloop, matad]").unwrap();
        assert_eq!(parsed, serde_json::json!(["alphaloop", "matad"]));
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
            target: ProcessArgs {
                process: None,
                integrand_name: None,
            },
            values: vec![ModelSetValue::Assignment(KvPair {
                key: internal_param.clone(),
                value: "1.0".to_string(),
            })],
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
            .find(|a| {
                matches!(
                    state
                        .model
                        .get_parameter(a.0.get_stripped_name())
                        .parameter_type,
                    ParameterType::Imaginary
                )
            })
            .expect("expected at least one external model parameter")
            .to_string();

        Set::Model {
            target: ProcessArgs {
                process: None,
                integrand_name: None,
            },
            values: vec![ModelSetValue::Assignment(KvPair {
                key: parameter.clone(),
                value: "-1.0e13-33.0e12i".to_string(),
            })],
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
            target: ProcessArgs {
                process: None,
                integrand_name: None,
            },
            values: vec![ModelSetValue::Assignment(KvPair {
                key: parameter,
                value: "1.0+2.0i".to_string(),
            })],
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
        let _guard = crate::LOG_TEST_MUTEX
            .lock()
            .unwrap_or_else(|err| err.into_inner());

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
        let _guard = crate::LOG_TEST_MUTEX
            .lock()
            .unwrap_or_else(|err| err.into_inner());

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
