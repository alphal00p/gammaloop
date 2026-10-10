//! DOT import: a finalized diagram together with the input data it does not store.

use super::*;

/// DOT attributes of a graph and of each of its vertices and edges, after defaults and
/// invisible-node attributes have been merged into every statement.
///
/// The maps are indexed by the finalized [`VertexId`] and [`EdgeId`]. An application can
/// read its own attributes from them, for example evaluator settings that FeynKit does
/// not interpret. Sewing merges the two legs of a cross-section initial state into one
/// edge; the incoming leg's values take precedence.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DotAttributes {
    pub graph: BTreeMap<String, String>,
    pub vertices: Vec<BTreeMap<String, String>>,
    pub edges: Vec<BTreeMap<String, String>>,
}

/// The local numerator fragments and the projector a DOT input wrote explicitly.
///
/// Every other fragment is one, and an omitted projector is one, until the model's
/// Feynman rules fill them in. The annotated dialect of [`FeynmanDiagram::to_dot`]
/// records a finalized diagram, so all of its factors are explicit.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ExplicitFactors {
    pub vertices: BTreeSet<VertexId>,
    pub edges: BTreeSet<EdgeId>,
    pub projector: bool,
}

/// A diagram read from DOT, with the attributes and explicit factors of its input.
#[derive(Debug, Clone)]
pub struct DotDiagram {
    pub diagram: FeynmanDiagram,
    pub attributes: DotAttributes,
    pub explicit: ExplicitFactors,
}

impl DotDiagram {
    /// Parse compact model-aware DOT or the annotated dialect emitted by
    /// [`FeynmanDiagram::to_dot`].
    pub fn from_dot(model: impl Into<Arc<Model>>, input: &str) -> Result<Self, DiagramError> {
        let parsed: DotGraph = DotGraph::from_string(input)
            .map_err(|error| DiagramError::DotParse(error.to_string()))?;
        Self::from_parsed(model.into(), parsed)
    }

    /// Parse every graph of one DOT document, in document order.
    ///
    /// Each graph may use either dialect. Annotated exports must carry the
    /// fingerprint of `model`.
    pub fn from_dot_set(
        model: impl Into<Arc<Model>>,
        input: &str,
    ) -> Result<Vec<Self>, DiagramError> {
        let model = model.into();
        let parsed = DotGraphSet::from_string(input)
            .map_err(|error| DiagramError::DotParse(error.to_string()))?;
        parsed
            .set
            .into_iter()
            .zip(parsed.global_data)
            .map(|(graph, global_data)| {
                Self::from_parsed(Arc::clone(&model), DotGraph { graph, global_data })
            })
            .collect()
    }

    fn from_parsed(model: Arc<Model>, parsed: DotGraph) -> Result<Self, DiagramError> {
        if [
            "model_fingerprint",
            "feynkit_name",
            "half_edge_order",
            "loop_momentum_basis",
            "cuts",
            "topology_threshold_candidates",
        ]
        .iter()
        .any(|key| parsed.global_data.statements.contains_key(*key))
        {
            return FeynmanDiagram::from_annotated_dot(model, parsed);
        }
        // A model particle index is meaningful only under the model fingerprint of the
        // annotated dialect. Older GammaLoop exports wrote one without it.
        if let Some((_, edge, _)) = parsed
            .iter_edges()
            .find(|(_, _, edge)| edge.data.statements.contains_key("particle_id"))
        {
            return Err(DiagramError::InvalidDotAttribute {
                target: format!("edge {}", edge.0),
                attribute: "particle_id",
                value: "a model particle index needs the model_fingerprint of FeynKit's \
                        annotated DOT; name the particle with `particle` or `pdg`"
                    .into(),
            });
        }
        FeynmanDiagram::from_compact_dot(model, parsed)
    }
}
