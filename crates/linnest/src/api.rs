use figment::{providers::Serialized, Figment, Profile};
use linnet::parser::set::DotGraphSet;
use serde::{Deserialize, Serialize};

use crate::{CBORTypstGraph, TypstGraph};

#[cfg(target_arch = "wasm32")]
use wasm_minimal_protocol::*;

fn decode_cbor_map(arg: &[u8]) -> Result<ciborium::Value, String> {
    ciborium::de::from_reader(arg).map_err(|e| format!("Failed to deserialize CBOR value: {}", e))
}

fn decode_dot_string(arg: &[u8]) -> Result<&str, String> {
    let dot_string = match std::str::from_utf8(arg) {
        Ok(s) => s,
        Err(_) => return Err("Invalid UTF-8".to_string()),
    };
    Ok(dot_string)
}

pub fn parse_dot_graphs_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let dot_string = decode_dot_string(arg)?;
    let dots = DotGraphSet::from_string(dot_string).map_err(|a| a.to_string())?;
    dots.to_rkyv_bytes::<4096>().map(|bytes| bytes.to_vec())
}

pub fn layout_parsed_graphs_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let cbor_map = decode_cbor_map(arg2)?;
    let figment = Figment::from(Serialized::from(cbor_map.clone(), Profile::Default));
    let mut dots = unsafe { rkyv::from_bytes_unchecked::<DotGraphSet>(arg) }
        .map_err(|e| format!("Failed to deserialize archived dot graphs: {}", e))?;

    for global_data in &mut dots.global_data {
        global_data.set_figment(figment.clone());
    }

    let mut graphs = Vec::new();
    for g in dots.into_iter() {
        let mut typst_graph = TypstGraph::from_dot(g, &Figment::new());

        typst_graph.layout();
        graphs.push((
            typst_graph.to_cbor(),
            typst_graph.to_dot_graph().debug_dot(),
        ));
    }

    let mut buffer = Vec::new();
    ciborium::ser::into_writer(&TypstOutput { graphs, cbor_map }, &mut buffer)
        .map_err(|a| a.to_string())?;
    Ok(buffer)
}

pub fn layout_graph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let parsed = parse_dot_graphs_bytes(arg)?;
    layout_parsed_graphs_bytes(&parsed, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn parse_graph(arg: &[u8]) -> Result<Vec<u8>, String> {
    parse_dot_graphs_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn layout_parsed_graph(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    layout_parsed_graphs_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn layout_graph(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    layout_graph_bytes(arg, arg2)
}

#[derive(Debug, Serialize, Deserialize)]
pub struct TypstOutput {
    pub(crate) graphs: Vec<(CBORTypstGraph, String)>,
    pub(crate) cbor_map: ciborium::Value,
}
