use crate::{
    graph::LoopMomentumBasis,
    uv::{
        ApproximationType,
        approx::{ForestNodeLike, UVCtx},
    },
};
use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Local3DLoopRescaling {
    FullSubgraph,
    ReducedSubgraph,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub(crate) enum Local3DConceptualBranch {
    #[serde(rename = "U")]
    U,
    #[serde(rename = "S")]
    S,
    #[serde(rename = "US")]
    US,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub(crate) struct Local3DSignedBranch {
    pub branch: Local3DConceptualBranch,
    pub coefficient: i8,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Local3DMaterialization {
    DirectU,
    /// The exact linear identity H(X) = S(X) + U(X-S(X)).
    FactorizedSoft,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub(crate) struct Local3DRouteSignature {
    pub edge: usize,
    pub loop_signature: String,
    pub external_signature: String,
}

/// Atom-free record of one direct operation on an orientation-local CFF.
///
/// `active_subgraph` is the accumulated local sector, while
/// `rescaled_subgraph` is the part acted on by this step. They differ during
/// disconnected replay. `canonical_route_*` records the component route and
/// `route_*` records the exact local LMB selected by the same code path that
/// constructs the Laurent projection.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub(crate) struct Local3DProjectionStep {
    pub topo_order: usize,
    pub current_component: String,
    pub given_component: String,
    pub active_subgraph: String,
    pub rescaled_subgraph: String,
    pub scheme: ApproximationType,
    pub dod: i32,
    pub rescaling: Local3DLoopRescaling,
    pub canonical_route_loop_edges: Vec<usize>,
    pub canonical_route_external_edges: Vec<usize>,
    pub canonical_route_signatures: Vec<Local3DRouteSignature>,
    pub route_loop_edges: Vec<usize>,
    pub route_external_edges: Vec<usize>,
    pub route_signatures: Vec<Local3DRouteSignature>,
    pub conceptual_branches: Vec<Local3DSignedBranch>,
    pub materialization: Local3DMaterialization,
}

/// One ordered projection history contributing to a local node. Nested steps
/// are stored child-to-parent; disconnected and local/integrated replay
/// sectors remain distinct paths.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub(crate) struct Local3DProjectionPath {
    pub steps: Vec<Local3DProjectionStep>,
}

impl Local3DProjectionPath {
    fn with_step(&self, step: &Local3DProjectionStep) -> Self {
        let mut path = self.clone();
        path.steps.push(step.clone());
        path
    }
}

pub(in crate::uv) fn paths_with_step(
    paths: &[Local3DProjectionPath],
    step: &Local3DProjectionStep,
) -> Vec<Local3DProjectionPath> {
    if paths.is_empty() {
        vec![Local3DProjectionPath {
            steps: vec![step.clone()],
        }]
    } else {
        paths.iter().map(|path| path.with_step(step)).collect()
    }
}

impl Local3DLoopRescaling {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::uv) fn projection_step<S: ForestNodeLike>(
        self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        active_subgraph: &SuBitGraph,
        rescaled_subgraph: Option<&SuBitGraph>,
        canonical_lmb: &LoopMomentumBasis,
        lmb: &LoopMomentumBasis,
    ) -> Local3DProjectionStep {
        let rescaled_subgraph = rescaled_subgraph.unwrap_or_else(|| current.subgraph());
        let route_signatures = |lmb: &LoopMomentumBasis, subgraph: &SuBitGraph| {
            ctx.graph
                .paired_edges(subgraph)
                .into_iter()
                .map(|edge| {
                    let signature = &lmb.edge_signatures[edge];
                    Local3DRouteSignature {
                        edge: usize::from(edge),
                        loop_signature: signature.internal.to_string(),
                        external_signature: signature.external.to_string(),
                    }
                })
                .collect::<Vec<_>>()
        };
        let (conceptual_branches, materialization) =
            if current.renormalization_scheme() == ApproximationType::IR && current.dod() > 0 {
                (
                    vec![
                        Local3DSignedBranch {
                            branch: Local3DConceptualBranch::U,
                            coefficient: 1,
                        },
                        Local3DSignedBranch {
                            branch: Local3DConceptualBranch::S,
                            coefficient: 1,
                        },
                        Local3DSignedBranch {
                            branch: Local3DConceptualBranch::US,
                            coefficient: -1,
                        },
                    ],
                    Local3DMaterialization::FactorizedSoft,
                )
            } else {
                (
                    vec![Local3DSignedBranch {
                        branch: Local3DConceptualBranch::U,
                        coefficient: 1,
                    }],
                    Local3DMaterialization::DirectU,
                )
            };

        Local3DProjectionStep {
            topo_order: current.topo_order(),
            current_component: current.subgraph().string_label(),
            given_component: given.subgraph().string_label(),
            active_subgraph: active_subgraph.string_label(),
            rescaled_subgraph: rescaled_subgraph.string_label(),
            scheme: current.renormalization_scheme(),
            dod: current.dod(),
            rescaling: self,
            canonical_route_loop_edges: canonical_lmb
                .loop_edges
                .iter()
                .map(|edge| usize::from(*edge))
                .collect(),
            canonical_route_external_edges: canonical_lmb
                .ext_edges
                .iter()
                .map(|edge| usize::from(*edge))
                .collect(),
            canonical_route_signatures: route_signatures(canonical_lmb, current.subgraph()),
            route_loop_edges: lmb
                .loop_edges
                .iter()
                .map(|edge| usize::from(*edge))
                .collect(),
            route_external_edges: lmb
                .ext_edges
                .iter()
                .map(|edge| usize::from(*edge))
                .collect(),
            route_signatures: route_signatures(lmb, rescaled_subgraph),
            conceptual_branches,
            materialization,
        }
    }
}
