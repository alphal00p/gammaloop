use figment::{providers::Serialized, Figment, Profile};
use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};
use linnet::parser::{set::DotGraphSet, DotGraph};

use crate::{
    graph_api::{decode_graph_bytes_list, encode_cbor},
    TypstGraph,
};

#[cfg(target_arch = "wasm32")]
use crate::{__ToResult, __send_result_to_host, __write_args_to_buffer};
#[cfg(target_arch = "wasm32")]
use wasm_minimal_protocol::*;

fn decode_cbor_map(arg: &[u8]) -> Result<ciborium::Value, String> {
    ciborium::de::from_reader(arg).map_err(|e| format!("Failed to deserialize CBOR value: {e}"))
}

fn take_layout_subgraph(cbor_map: &mut ciborium::Value) -> Result<Option<String>, String> {
    let ciborium::Value::Map(entries) = cbor_map else {
        return Ok(None);
    };

    let mut subgraph = None;
    entries.retain(|(key, value)| {
        let is_subgraph = matches!(key, ciborium::Value::Text(key) if key == "subgraph");
        if is_subgraph {
            subgraph = Some(value.clone());
        }
        !is_subgraph
    });

    let Some(value) = subgraph else {
        return Ok(None);
    };

    match value {
        ciborium::Value::Null => Ok(None),
        ciborium::Value::Text(label) if label == "none" => Ok(None),
        ciborium::Value::Text(label) => Ok(Some(label)),
        other => Err(format!(
            "layout subgraph must be a base62 string or none, got {other:?}"
        )),
    }
}

fn decode_dot_string(arg: &[u8]) -> Result<&str, String> {
    std::str::from_utf8(arg).map_err(|_| "Invalid UTF-8".to_string())
}

pub fn parse_dot_graphs_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let dot_string = decode_dot_string(arg)?;
    let dots = DotGraphSet::from_string(dot_string).map_err(|err| err.to_string())?;
    let graph_bytes = dots.into_graph_bytes_set::<4096>()?.graphs;
    encode_cbor(&graph_bytes)
}

pub fn layout_parsed_graph_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let mut cbor_map = decode_cbor_map(arg2)?;
    let subgraph_label = take_layout_subgraph(&mut cbor_map)?;
    let figment = Figment::from(Serialized::from(cbor_map, Profile::Default));

    let mut graph = unsafe { rkyv::from_bytes_unchecked::<DotGraph>(arg) }
        .map_err(|err| format!("Failed to deserialize archived dot graph: {err}"))?;
    graph.global_data.set_figment(figment.clone());
    let subgraph = subgraph_label
        .as_deref()
        .map(|label| {
            SuBitGraph::from_base62(label, graph.n_hedges()).ok_or_else(|| {
                format!(
                    "Invalid subgraph label {label:?} for graph with {} half edges",
                    graph.n_hedges()
                )
            })
        })
        .transpose()?;

    let mut typst_graph = TypstGraph::from_dot(graph, &figment);
    typst_graph.layout_with_subgraph(subgraph.as_ref())?;
    typst_graph
        .to_dot_graph()
        .to_rkyv_bytes::<4096>()
        .map(|bytes| bytes.to_vec())
}

pub fn layout_parsed_graphs_bytes(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    let graph_bytes = decode_graph_bytes_list(arg)?;
    let mut graphs = Vec::with_capacity(graph_bytes.len());
    for graph_bytes in graph_bytes {
        graphs.push(layout_parsed_graph_bytes(&graph_bytes, arg2)?);
    }

    encode_cbor(&graphs)
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
pub fn graph_from_spec(arg: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_from_spec_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_with_payloads(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_with_payloads_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn layout_parsed_graph(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    layout_parsed_graph_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn layout_graph(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    layout_graph_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_info(arg: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_info_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_dot(arg: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_dot_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_nodes(arg: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_nodes_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_nodes_of(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_nodes_of_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_nodes_of_subgraph(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_nodes_of_archived_subgraph_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_edges(arg: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_edges_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_edges_of(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_edges_of_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_edges_of_subgraph(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_edges_of_archived_subgraph_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_subgraph(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_subgraph_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_archived_subgraph(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_archived_subgraph_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_compass_subgraph(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_compass_subgraph_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_archived_compass_subgraph(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_archived_compass_subgraph_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn subgraph_label(arg: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::subgraph_label_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn subgraph_hedges(arg: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::subgraph_hedges_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn subgraph_contains_hedge(arg: &[u8], arg2: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::subgraph_contains_hedge_bytes(arg, arg2)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_cycle_basis(arg: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_cycle_basis_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_spanning_forests(arg: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_spanning_forests_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_join_by_edge_key(arg: &[u8], arg2: &[u8], arg3: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_join_by_edge_key_bytes(arg, arg2, arg3)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn graph_join_by_hedge_key(arg: &[u8], arg2: &[u8], arg3: &[u8]) -> Result<Vec<u8>, String> {
    crate::graph_api::graph_join_by_hedge_key_bytes(arg, arg2, arg3)
}
