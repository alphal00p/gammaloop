//! Construct inspection graphs through FeynKit's authenticated saved-DOT path.

use std::collections::BTreeMap;

use feynkit_graph::{EdgeId, FeynmanDiagram};
use feynkit_model::Model;
use serde_json::Value;

use super::InspectionError;
use crate::{VakintError, rustred_evaluation::data::LoopData};

pub(super) fn master_graph(
    data: &LoopData,
    metadata: &Value,
) -> Result<FeynmanDiagram, VakintError> {
    let master = metadata["master_id"].as_str().unwrap_or("unknown master");
    let unavailable = |detail: String| InspectionError::UnavailableGraph {
        master: master.into(),
        detail,
    };
    let read_text = |path: &str| {
        data.read(path)
            .map_err(|error| unavailable(error.to_string()))
            .and_then(|bytes| {
                String::from_utf8(bytes).map_err(|error| unavailable(error.to_string()))
            })
    };
    let model_path = metadata["model"]
        .as_str()
        .ok_or_else(|| unavailable("no authenticated FeynKit model is declared".into()))?;
    let graph_path = metadata["files"]["graph.dot"]["path"]
        .as_str()
        .ok_or_else(|| unavailable("no master-specific DOT payload is declared".into()))?;
    let model = Model::from_json(&read_text(model_path)?)
        .map_err(|error| unavailable(error.to_string()))?;
    let diagram = FeynmanDiagram::from_dot(model, &read_text(graph_path)?)
        .map_err(|error| unavailable(error.to_string()))?;

    // The integral definition and the native graph are independently
    // authenticated payloads. They must describe the same physical powers;
    // replacing one with the other would conceal inconsistent data.
    if let Some(rows) = metadata["graph_definition"]["edge_powers"].as_array() {
        let mut powers = BTreeMap::new();
        for row in rows {
            let edge = row[0].as_u64().and_then(|id| usize::try_from(id).ok());
            let power = row[1].as_i64().and_then(|power| i32::try_from(power).ok());
            let (Some(edge), Some(power)) = (edge, power) else {
                return Err(unavailable("invalid physical edge power".into()).into());
            };
            if power <= 0 || powers.insert(EdgeId(edge), power).is_some() {
                return Err(
                    unavailable("physical edge powers must be positive and unique".into()).into(),
                );
            }
        }
        if powers
            .keys()
            .copied()
            .ne(diagram.edges().map(|(id, _, _)| id))
        {
            return Err(unavailable(
                "integral powers do not cover exactly the graph's physical edges".into(),
            )
            .into());
        }
        if diagram
            .edges()
            .any(|(id, _, edge)| edge.propagator_power != powers[&id])
        {
            return Err(unavailable(
                "master DOT propagator powers disagree with its integral definition".into(),
            )
            .into());
        }
    }
    diagram
        .validate()
        .map_err(|error| unavailable(error.to_string()))?;
    Ok(diagram)
}
