//! The owned structural endpoint of generation, before symbolic Feynman rules.

use feynkit_graph::{
    DiagramError, DiagramId, EdgeEndpoints, EdgeId, ExternalLeg, FeynmanDiagram, LoopMomentumBasis,
    VertexId, VertexSlot,
};
use feynkit_model::{Model, ParticleId, VertexRuleId};

use super::{FermionSigns, GenerationReport};

/// Generated amplitude diagrams with exact weights and a complete momentum routing.
/// No numerator-based removal or grouping is performed.
#[derive(Debug, Clone)]
pub struct RoutedGenerationResult {
    pub diagrams: Vec<SignedRoutedDiagram>,
    pub report: GenerationReport,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct DiagramWeight {
    pub(super) sign: i8,
    pub(super) multiplicity: u64,
}

impl DiagramWeight {
    pub(super) fn new(signs: FermionSigns, multiplicity: u64) -> Self {
        let negative = signs.internal_loop_negative
            ^ (signs.include_external_ordering && signs.external_ordering_negative)
            ^ (signs.include_antifermion_spin_sum && signs.antifermion_spin_sum_negative);
        Self {
            sign: if negative { -1 } else { 1 },
            multiplicity,
        }
    }
}

/// A finalized physical topology, including its Grassmann sign and momentum basis.
///
/// The exact diagram weight is `sign * multiplicity / automorphism_order`.
/// Local Lorentz/color numerators and external-state projectors have not been
/// constructed. The private geometry owner shares finalization and validation
/// with full Feynman diagrams without exposing unset fragments as physical data.
#[derive(Debug, Clone)]
pub struct SignedRoutedDiagram {
    geometry: FeynmanDiagram,
    weight: DiagramWeight,
}

/// Structural metadata of one interaction vertex.
#[derive(Debug, Clone, Copy)]
pub struct RoutedVertex<'a> {
    pub name: &'a str,
    pub interaction: Option<VertexRuleId>,
}

/// Structural metadata of one edge in its stored source-to-target orientation.
#[derive(Debug, Clone, Copy)]
pub struct RoutedEdge<'a> {
    pub particle: ParticleId,
    pub directed: bool,
    pub external: Option<&'a ExternalLeg>,
    pub source_slot: VertexSlot,
    pub target_slot: VertexSlot,
}

impl SignedRoutedDiagram {
    pub(super) fn new(geometry: FeynmanDiagram, weight: DiagramWeight) -> Self {
        Self { geometry, weight }
    }

    pub fn model(&self) -> &Model {
        self.geometry.model()
    }

    pub fn name(&self) -> Option<&str> {
        self.geometry.name()
    }

    pub fn id(&self) -> Result<DiagramId, DiagramError> {
        self.geometry.id()
    }

    pub fn sign(&self) -> i8 {
        self.weight.sign
    }

    pub fn multiplicity(&self) -> u64 {
        self.weight.multiplicity
    }

    /// The denominator of the native diagram's exact symmetry weight.
    pub fn automorphism_order(&self) -> u64 {
        self.geometry.symmetry_factor()
    }

    pub fn vertices(&self) -> impl Iterator<Item = (VertexId, RoutedVertex<'_>)> {
        self.geometry.vertices().map(|(id, data)| {
            (
                id,
                RoutedVertex {
                    name: &data.name,
                    interaction: data.interaction,
                },
            )
        })
    }

    pub fn edges(&self) -> impl Iterator<Item = (EdgeId, EdgeEndpoints, RoutedEdge<'_>)> {
        self.geometry.edges().map(|(id, ends, data)| {
            (
                id,
                ends,
                RoutedEdge {
                    particle: data.particle,
                    directed: data.directed,
                    external: data.external.as_ref(),
                    source_slot: data.source_slot(),
                    target_slot: data.target_slot(),
                },
            )
        })
    }

    /// Independent loop and external carriers, dependent external momenta, and
    /// exact signed coefficients for every edge in its source-to-target frame.
    pub fn loop_momentum_basis(&self) -> &LoopMomentumBasis {
        self.geometry.loop_momentum_basis()
    }

    pub fn validate(&self) -> Result<(), DiagramError> {
        self.geometry.validate()
    }
}
