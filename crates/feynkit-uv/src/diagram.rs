//! UV combinatorics of FeynKit diagrams.

use feynkit_graph::{
    DiagramEdge, DiagramError, DiagramVertex, EdgeId, FeynmanDiagram, LoopMomentumBasis,
    MomentumSignature,
};
use feynkit_model::{Model, Particle};
use linnet::half_edge::involution::EdgeIndex;
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetOps};
use linnet::half_edge::{HedgeGraph, NoData};

use crate::UvGraph;

/// Spacetime dimension of FeynKit's UV power counting.
const DIMENSION: i32 = 4;

/// FeynKit diagrams enumerate UV subgraphs among their internal propagators,
/// classify with the hard-loop degree of
/// [`FeynmanDiagram::uv_degree_of_divergence_of`], prune with the larger of that
/// degree and the local superficial degree of
/// [`FeynmanDiagram::superficial_degree_of_divergence_of`], and route spinneys
/// with FeynKit's [`LoopMomentumBasis`].
///
/// Particle identities are model-local, so particle queries accept only the
/// diagram's own model (by fingerprint) and fail for any other.
impl UvGraph for FeynmanDiagram {
    type Edge = DiagramEdge;
    type Vertex = DiagramVertex;
    type HedgeData = NoData;
    type Basis = LoopMomentumBasis;
    type Error = DiagramError;

    fn hedge_graph(&self) -> &HedgeGraph<DiagramEdge, DiagramVertex> {
        self.underlying()
    }

    /// Internal propagators only: dummy attachments and external legs, including
    /// sewn initial states, never become UV hard.
    fn uv_domain(&self, subgraph: &SuBitGraph) -> SuBitGraph {
        subgraph.intersection(&self.internal_subgraph())
    }

    fn is_dummy_edge(&self, edge: EdgeIndex) -> bool {
        self.underlying()[edge].is_dummy
    }

    /// Dummy attachments carry no particle for classification.
    fn edge_pdg_code(&self, edge: EdgeIndex, model: &Model) -> Result<Option<isize>, DiagramError> {
        if self.is_dummy_edge(edge) {
            return Ok(None);
        }
        let pdg_code = self.edge_particle(edge, model)?.pdg_code;
        Ok(Some(
            pdg_code
                .try_into()
                .expect("validated particle PDG codes must fit in isize"),
        ))
    }

    /// Dummy attachments are never massive for classification.
    fn is_massive_edge(&self, edge: EdgeIndex, model: &Model) -> Result<bool, DiagramError> {
        if self.is_dummy_edge(edge) {
            return Ok(false);
        }
        Ok(!self.edge_particle(edge, model)?.is_massless(model))
    }

    /// The larger of the superficial local degree and the classifying degree, so
    /// the pruning bound holds by construction.
    fn local_degree_of_divergence(&self, subgraph: &SuBitGraph) -> Result<i32, DiagramError> {
        let superficial = self.superficial_degree_of_divergence_of(subgraph, DIMENSION)?;
        Ok(superficial.max(self.uv_degree_of_divergence_of(subgraph, DIMENSION)?))
    }

    /// Decides the sign without the classifying degree when the superficial
    /// degree already passes.
    fn passes_pruning(&self, subgraph: &SuBitGraph) -> Result<bool, DiagramError> {
        Ok(
            self.superficial_degree_of_divergence_of(subgraph, DIMENSION)? >= 0
                || self.uv_degree_of_divergence_of(subgraph, DIMENSION)? >= 0,
        )
    }

    fn degree_of_divergence(&self, subgraph: &SuBitGraph) -> Result<i32, DiagramError> {
        self.uv_degree_of_divergence_of(subgraph, DIMENSION)
    }

    fn empty_basis(&self) -> LoopMomentumBasis {
        LoopMomentumBasis {
            tree_edges: Vec::new(),
            loop_edges: Vec::new(),
            external_edges: Vec::new(),
            dependent_externals: Vec::new(),
            edge_signatures: self
                .underlying()
                .iter_edges()
                .map(|(_, edge, _)| (EdgeId(edge.0), MomentumSignature::default()))
                .collect(),
        }
    }

    fn compatible_sub_basis(
        &self,
        subgraph: &InternalSubGraph,
        parent: &LoopMomentumBasis,
    ) -> Result<LoopMomentumBasis, DiagramError> {
        self.compatible_momentum_basis_of(&subgraph.filter, parent)
    }
}

/// Particle lookup bound to the diagram's own model.
trait DiagramParticleExt {
    fn edge_particle<'m>(
        &self,
        edge: EdgeIndex,
        model: &'m Model,
    ) -> Result<&'m Particle, DiagramError>;
}

impl DiagramParticleExt for FeynmanDiagram {
    /// The particle of `edge`, resolved in `model` after checking that `model` is
    /// the diagram's own: particle identifiers are only meaningful there.
    fn edge_particle<'m>(
        &self,
        edge: EdgeIndex,
        model: &'m Model,
    ) -> Result<&'m Particle, DiagramError> {
        let own = self.model().fingerprint();
        let supplied = model.fingerprint();
        if own != supplied {
            return Err(DiagramError::ModelFingerprintMismatch {
                serialized: own,
                actual: supplied,
            });
        }
        Ok(model.particle_by_id(self.underlying()[edge].particle)?)
    }
}
