//! Neutral wire types for source-generated command, settings, and topology data.

use serde::{Deserialize, Serialize};
use serde_json::Value;

pub const GENERATED_REFERENCE_SCHEMA: u32 = 3;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct GeneratedSource {
    pub kind: String,
    pub path: String,
    pub symbol: String,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct GammaLoopReference {
    pub schema_version: u32,
    pub sources: Vec<GeneratedSource>,
    pub commands: Vec<CliCommand>,
    pub settings: Vec<SettingReference>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct CliCommand {
    pub path: String,
    pub name: String,
    pub usage: String,
    pub about: String,
    pub hidden: bool,
    pub generated_help: bool,
    pub short_flag: Option<char>,
    pub long_flag: Option<String>,
    pub aliases: Vec<CliAlias>,
    pub arguments: Vec<CliArgument>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum CliAliasKind {
    Name,
    ShortFlag,
    LongFlag,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct CliAlias {
    pub name: String,
    pub kind: CliAliasKind,
    pub visible: bool,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum CliArgumentAction {
    Set,
    Append,
    SetTrue,
    SetFalse,
    Count,
    Help,
    HelpShort,
    HelpLong,
    Version,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct CliValueArity {
    pub min: usize,
    pub max: Option<usize>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct CliArgument {
    pub id: String,
    pub long: Option<String>,
    pub short: Option<char>,
    pub aliases: Vec<CliAlias>,
    pub action: CliArgumentAction,
    pub arity: CliValueArity,
    pub takes_values: bool,
    pub value_required: bool,
    pub value_names: Vec<String>,
    pub help: String,
    pub required: bool,
    pub positional: bool,
    pub index: Option<usize>,
    pub hidden: bool,
    pub global: bool,
    pub inherited: bool,
    pub exclusive: bool,
    pub require_equals: bool,
    pub value_delimiter: Option<char>,
    pub value_terminator: Option<String>,
    pub conflicts_with: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub requires: Vec<String>,
    pub defaults: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub default_missing_values: Vec<String>,
    pub possible_values: Vec<String>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct SettingReference {
    pub path: String,
    pub value_type: String,
    pub description: String,
    pub required: bool,
    pub default: Option<Value>,
    pub possible_values: Vec<String>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct VakintReference {
    pub schema_version: u32,
    pub sources: Vec<GeneratedSource>,
    pub dependencies: Vec<DependencyReference>,
    pub topologies: Vec<TopologyReference>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct DependencyReference {
    pub name: String,
    pub minimum_version: String,
    pub source_symbol: String,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct TopologyReference {
    pub name: String,
    pub loops: usize,
    pub propagator_slots: usize,
}
