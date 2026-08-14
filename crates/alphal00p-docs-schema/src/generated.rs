//! Neutral wire types for source-generated command, settings, and topology data.

use serde::{Deserialize, Serialize};
use serde_json::Value;

pub const GENERATED_REFERENCE_SCHEMA: u32 = 1;

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
    pub about: String,
    pub aliases: Vec<String>,
    pub arguments: Vec<CliArgument>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct CliArgument {
    pub id: String,
    pub long: Option<String>,
    pub short: Option<char>,
    pub value_names: Vec<String>,
    pub help: String,
    pub required: bool,
    pub defaults: Vec<String>,
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
