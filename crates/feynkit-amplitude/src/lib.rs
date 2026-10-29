//! Diagram-backed symbolic amplitudes, independent of any integration runtime.
//!
//! Expressions are amputated, unintegrated Feynman-rule operators. External
//! wavefunctions are represented by labeled open tensor ports. Squaring keeps
//! the two copies distinct until explicit spin and color completeness sums.

#![forbid(unsafe_code)]

mod color;
mod spin;

pub use color::{ColorRepresentation, ColorSum, ColorSumError};
pub use spin::{AxialReference, SpinSum, SpinSumError};

use feynkit_graph::{
    DiagramError, ExternalState, FeynmanDiagram, expressions::evaluate_overall_factor, symbols,
};
use feynkit_model::{Model, ModelError, ParameterType, ParticleId};
use idenso::{
    CookMode, CookSettings,
    dirac::GammaSimplifySettings,
    representations::Bispinor,
    tensor::{SymbolicTensor, inference::TensorInferenceError},
};
use linnet::half_edge::involution::HedgePair;
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        abstract_index::AbstractIndex,
        dimension::Dimension,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        representation::{LibrarySlot, Minkowski},
        slot::{DualSlotTo, IsAbstractSlot, ParseableAind},
    },
};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};
use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::id::{Pattern, Replacement};
use symbolica::{function, symbol};
use thiserror::Error;

/// Invalid diagram collection or unsupported external-state operation.
#[derive(Debug, Error)]
pub enum AmplitudeError {
    #[error("an amplitude needs at least one diagram")]
    Empty,
    #[error("diagram '{0}' belongs to a different model")]
    DifferentModel(String),
    #[error("diagram '{0}' has incompatible external states or tensor ports")]
    DifferentExternals(String),
    #[error("diagram '{0}' is already sewn; construct an amplitude from unsewn diagrams")]
    SewnDiagram(String),
    #[error("cannot identify external tensor ports in diagram '{0}': {1}")]
    ExternalPorts(String, String),
    #[error("external leg {0} does not exist")]
    UnknownLeg(usize),
    #[error("{kind} states of external leg {index} have already been summed")]
    AlreadySummed { kind: &'static str, index: usize },
    #[error("invalid tensor expression: {0}")]
    Tensor(String),
    #[error(transparent)]
    Diagram(#[from] DiagramError),
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error(transparent)]
    Spin(#[from] SpinSumError),
    #[error(transparent)]
    Color(#[from] ColorSumError),
}

impl From<TensorInferenceError> for AmplitudeError {
    fn from(error: TensorInferenceError) -> Self {
        Self::Tensor(error.to_string())
    }
}

impl From<idenso::tensor::composition::TensorCompositionError> for AmplitudeError {
    fn from(error: idenso::tensor::composition::TensorCompositionError) -> Self {
        Self::Tensor(error.to_string())
    }
}

impl From<idenso::AdjointError> for AmplitudeError {
    fn from(error: idenso::AdjointError) -> Self {
        Self::Tensor(error.to_string())
    }
}

/// Choices shared by every term. Model-declared real parameters and physical
/// momenta are real automatically; extra scalar assumptions are explicit.
#[derive(Clone, Debug)]
pub struct AmplitudeOptions {
    pub dimension: Dimension,
    pub real: Vec<Atom>,
}

impl Default for AmplitudeOptions {
    fn default() -> Self {
        Self {
            dimension: Dimension::Concrete(4),
            real: Vec::new(),
        }
    }
}

/// A physical external state, independent of a diagram's half-edge numbering.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AmplitudeLeg {
    pub index: usize,
    pub particle: ParticleId,
    pub state: ExternalState,
    pub slots: Vec<LibrarySlot<AbstractIndex>>,
    pub tensor_index: AbstractIndex,
}

impl AmplitudeLeg {
    /// Physical incoming/outgoing momentum, with the external leg's stable label.
    pub fn momentum(&self) -> Atom {
        symbols::external_momentum().call(self.index)
    }

    /// Common bare label for this leg's spin and color ports.
    pub fn index_atom(&self) -> Atom {
        self.tensor_index.to_atom()
    }

    fn color_slot(&self) -> Option<LibrarySlot<AbstractIndex>> {
        self.slots
            .iter()
            .copied()
            .find(|s| s.rep_name() != Bispinor {}.into() && s.rep_name() != Minkowski {}.into())
    }
}

/// A coherent sum of weighted, amputated diagram operators.
///
/// Source diagrams and individual terms remain available. Overall factors are
/// applied exactly once; the diagnostic automorphism order is not another weight.
#[derive(Clone, Debug)]
pub struct Amplitude {
    diagrams: Vec<Arc<FeynmanDiagram>>,
    terms: Vec<SymbolicTensor<PartialStructure>>,
    legs: Vec<AmplitudeLeg>,
    options: AmplitudeOptions,
    conjugated: bool,
}

impl Amplitude {
    pub fn from_diagram(diagram: impl Into<Arc<FeynmanDiagram>>) -> Result<Self, AmplitudeError> {
        Self::from_diagrams([diagram.into()])
    }

    pub fn from_diagrams(
        diagrams: impl IntoIterator<Item = Arc<FeynmanDiagram>>,
    ) -> Result<Self, AmplitudeError> {
        Self::new(diagrams, AmplitudeOptions::default())
    }

    pub fn new(
        diagrams: impl IntoIterator<Item = Arc<FeynmanDiagram>>,
        options: AmplitudeOptions,
    ) -> Result<Self, AmplitudeError> {
        let diagrams: Vec<_> = diagrams.into_iter().collect();
        let first = diagrams.first().ok_or(AmplitudeError::Empty)?;
        let fingerprint = first.model().fingerprint();
        let mut terms = Vec::with_capacity(diagrams.len());
        let mut legs = Vec::new();
        for diagram in &diagrams {
            diagram.validate()?;
            if !Arc::ptr_eq(&first.model_arc(), &diagram.model_arc())
                && diagram.model().fingerprint() != fingerprint
            {
                return Err(AmplitudeError::DifferentModel(diagram.name().to_owned()));
            }
            let (term, external) = Self::diagram_term(diagram, &options)?;
            if terms.is_empty() {
                legs = external;
            } else if legs != external {
                return Err(AmplitudeError::DifferentExternals(
                    diagram.name().to_owned(),
                ));
            }
            terms.push(term);
        }
        Ok(Self {
            diagrams,
            terms,
            legs,
            options,
            conjugated: false,
        })
    }

    pub fn diagrams(&self) -> &[Arc<FeynmanDiagram>] {
        &self.diagrams
    }
    pub fn terms(&self) -> &[SymbolicTensor<PartialStructure>] {
        &self.terms
    }
    pub fn legs(&self) -> &[AmplitudeLeg] {
        &self.legs
    }
    pub fn model(&self) -> &Model {
        self.diagrams[0].model()
    }
    pub fn dimension(&self) -> Dimension {
        self.options.dimension
    }
    pub fn is_conjugated(&self) -> bool {
        self.conjugated
    }

    pub fn expression(&self) -> Result<SymbolicTensor<PartialStructure>, AmplitudeError> {
        let (first, rest) = self.terms.split_first().ok_or(AmplitudeError::Empty)?;
        rest.iter()
            .try_fold(first.clone(), |sum, term| Ok(sum.try_add(term)?))
    }

    /// Open tensor slots in physical external-leg order, retaining graph identities.
    pub fn structure(&self) -> PartialStructure {
        PartialStructure::from_logical_slots(
            self.legs
                .iter()
                .flat_map(|leg| &leg.slots)
                .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind))),
        )
    }

    /// Dirac adjunction preserves physical leg identities, including across
    /// diagrams with different open-fermion-chain pairings. Color conjugation
    /// and reversal of complex scalar coefficients are delegated to Idenso.
    pub fn conjugate(&self) -> Result<Self, AmplitudeError> {
        let mut result = self.clone();
        for leg in &mut result.legs {
            leg.slots = leg.slots.iter().map(|slot| slot.dual()).collect();
            leg.slots.sort();
        }
        let interface = result.structure();
        let settings = GammaSimplifySettings {
            gamma0: true,
            evaluate_traces: false,
            ..GammaSimplifySettings::default()
        };
        result.terms = self
            .terms
            .iter()
            .map(|term| {
                let adjoint = term.dirac_adjoint(true)?;
                let adjoint = self.apply_reality(&adjoint)?;
                // Physical leg labels and their declared order belong to the
                // amplitude. The shared scheduler owns gamma/metric cleanup.
                SymbolicTensor::from_declared_parts(adjoint.expression().clone(), interface.clone())
                    .and_then(|value| {
                        value.simplify_algebra(&idenso::tensor::AlgebraSettings {
                            gamma: Some(settings),
                            epsilon: settings.output == idenso::dirac::GammaOutput::Reduced,
                            ..Default::default()
                        })
                    })
                    .and_then(|value| {
                        value.contract(idenso::tensor::ContractSettings {
                            collect_chains: false,
                            collect_traces: false,
                            ..Default::default()
                        })
                    })
                    .map_err(|error| AmplitudeError::Tensor(error.to_string()))
            })
            .collect::<Result<_, AmplitudeError>>()?;
        result.conjugated = !self.conjugated;
        Ok(result)
    }

    pub fn squared(&self) -> Result<SquaredAmplitude, AmplitudeError> {
        // Keep a consistent ket/bra convention even when called on the adjoint.
        let amplitude = if self.conjugated {
            self.conjugate()?
        } else {
            self.clone()
        };
        let adjoint = amplitude.conjugate()?;
        // Loop variables in two independent loop integrals must not be identified.
        let loop_index = symbol!("feynkit_amplitude::loop_");
        let args = symbol!("feynkit_amplitude::loop_args___");
        let bra = adjoint.expression()?.wrap_indices(Self::bra(), false)?;
        let mapped = bra
            .expression()
            .replace(function!(symbols::loop_momentum(), loop_index, args))
            .with(function!(
                symbols::loop_momentum(),
                Self::bra().call(loop_index),
                args
            ));
        let bra = bra.with_checked_expression(mapped)?;
        let expression = amplitude.expression()?.multiply(&bra)?;
        Ok(SquaredAmplitude {
            amplitude,
            expression,
            spin_summed: BTreeSet::new(),
            color_summed: BTreeSet::new(),
        })
    }

    fn bra() -> Symbol {
        symbol!("feynkit_amplitude::bra")
    }

    fn apply_reality(
        &self,
        expression: &SymbolicTensor<PartialStructure>,
    ) -> Result<SymbolicTensor<PartialStructure>, AmplitudeError> {
        let conjugate = symbol!("spenso::conj");
        let mut rules = Vec::new();
        for parameter in self
            .model()
            .parameters()
            .iter()
            .filter(|p| p.parameter_type == ParameterType::Real)
        {
            let real = Atom::var(symbol!(&format!("UFO::{}", parameter.name)));
            // Model and user assumptions identify exact expressions. Their
            // names may end in underscores without becoming pattern wildcards.
            rules.push(Replacement::new(
                Pattern::Literal(conjugate.call(&real)),
                Pattern::Literal(real),
            ));
        }
        for real in &self.options.real {
            rules.push(Replacement::new(
                Pattern::Literal(conjugate.call(real)),
                Pattern::Literal(real.clone()),
            ));
        }
        let args = symbol!("feynkit_amplitude::args___");
        for head in [symbols::external_momentum(), symbols::loop_momentum()] {
            let vector = head.call(args);
            rules.push(Replacement::new(
                conjugate.call(&vector).to_pattern(),
                vector,
            ));
        }
        // The metric is real. Compact products of the known physical momentum
        // families are real even though arbitrary user tensor products need not be.
        let left = symbol!("feynkit_amplitude::left___");
        let right = symbol!("feynkit_amplitude::right___");
        for p in [symbols::external_momentum(), symbols::loop_momentum()] {
            for q in [symbols::external_momentum(), symbols::loop_momentum()] {
                let dot = function!(SPENSO_TAG.dot, p.call(left), q.call(right));
                rules.push(Replacement::new(conjugate.call(&dot).to_pattern(), dot));
            }
        }
        let dimension = self.options.dimension.to_symbolic();
        rules.push(Replacement::new(
            Pattern::Literal(conjugate.call(&dimension)),
            Pattern::Literal(dimension),
        ));
        Ok(expression.with_checked_expression(expression.expression().replace_multiple(rules))?)
    }

    fn diagram_term(
        diagram: &FeynmanDiagram,
        options: &AmplitudeOptions,
    ) -> Result<(SymbolicTensor<PartialStructure>, Vec<AmplitudeLeg>), AmplitudeError> {
        if !diagram.cuts().is_empty() {
            return Err(AmplitudeError::SewnDiagram(diagram.name().to_owned()));
        }
        let numerator = diagram
            .numerator()
            .multiply(diagram.numerator_prefactor())?
            .with_lorentz_dimension(options.dimension)?;
        let ports = numerator
            .structure()
            .logical_slots()
            .into_iter()
            .map(|slot| {
                let PartialIndex::Explicit(index) = slot.aind else {
                    return Err(AmplitudeError::ExternalPorts(
                        diagram.name().to_owned(),
                        "physical external ports need explicit leg identities".into(),
                    ));
                };
                Ok(slot.rep().slot::<AbstractIndex, _>(index))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let mut legs = BTreeMap::new();
        let mut identified = 0;
        for (pair, _, edge) in diagram.underlying().iter_edges() {
            let Some(external) = &edge.data.external else {
                continue;
            };
            let HedgePair::Unpaired { hedge, .. } = pair else {
                return Err(AmplitudeError::SewnDiagram(diagram.name().to_owned()));
            };
            let old_index = symbols::hedge_index().call((hedge.0, 1));
            let mut leg = AmplitudeLeg {
                index: external.index,
                particle: edge.data.particle,
                state: external.state,
                slots: Vec::new(),
                tensor_index: AbstractIndex::try_from(old_index.as_view())
                    .map_err(|e| AmplitudeError::Tensor(e.to_string()))?,
            };
            for slot in ports
                .iter()
                .filter(|s| cooking.uncook(s.aind.to_atom().as_view()) == old_index)
            {
                leg.tensor_index = slot.aind;
                leg.slots.push(*slot);
                identified += 1;
            }
            leg.slots.sort();
            let particle = diagram.model().particle_by_id(leg.particle)?;
            let expected = usize::from(particle.spin != 1) + usize::from(particle.color != 1);
            if leg.slots.len() != expected || legs.insert(leg.index, leg).is_some() {
                return Err(AmplitudeError::ExternalPorts(
                    diagram.name().to_owned(),
                    "missing, duplicate, or unsupported external slots".into(),
                ));
            }
        }
        if identified != ports.len() {
            return Err(AmplitudeError::ExternalPorts(
                diagram.name().to_owned(),
                "unassigned open tensor slots".into(),
            ));
        }
        let basis = diagram.loop_momentum_basis();
        let mut momentum_rules = Vec::new();
        let args = symbol!("feynkit_amplitude::momentum_args___");
        for (position, edge) in basis.external_edges.iter().enumerate() {
            let external = diagram
                .edges()
                .find(|(id, _, _)| id == edge)
                .and_then(|(_, _, e)| e.external.as_ref())
                .expect("validated external basis edge");
            momentum_rules.push(Replacement::new(
                function!(symbols::external_momentum(), position, args).to_pattern(),
                function!(symbols::external_momentum(), external.index, args),
            ));
        }
        let denominator = diagram.denominator_of_in_dimension(
            &diagram.internal_subgraph(),
            &BTreeMap::new(),
            options.dimension,
        )?;
        let [a, b, c, d] = ["den_edge_", "den_momentum_", "den_mass_", "den_value_"]
            .map(|s| symbol!(&format!("feynkit_amplitude::{s}")));
        let denominator = denominator
            .replace(function!(symbols::denominator(), a, b, c, d))
            .with(d);
        let factor = evaluate_overall_factor(diagram.overall_factor().as_view());
        let expression = diagram
            .model()
            .expand_couplings(&(factor * numerator.expression() / denominator));
        let expression = basis
            .route_expression(&expression)
            .replace_multiple(momentum_rules);
        let numerator = numerator.with_checked_expression(expression)?;
        let legs = legs.into_values().collect::<Vec<_>>();
        let interface = PartialStructure::from_logical_slots(
            legs.iter()
                .flat_map(|leg| &leg.slots)
                .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind))),
        );
        let numerator =
            SymbolicTensor::from_declared_parts(numerator.expression().clone(), interface)?;
        Ok((numerator, legs))
    }
}

/// A coherent amplitude times its adjoint, with independent tensor ports.
#[derive(Clone, Debug)]
pub struct SquaredAmplitude {
    amplitude: Amplitude,
    expression: SymbolicTensor<PartialStructure>,
    spin_summed: BTreeSet<usize>,
    color_summed: BTreeSet<usize>,
}

impl SquaredAmplitude {
    pub fn from_diagram(diagram: impl Into<Arc<FeynmanDiagram>>) -> Result<Self, AmplitudeError> {
        Amplitude::from_diagram(diagram)?.squared()
    }
    pub fn amplitude(&self) -> &Amplitude {
        &self.amplitude
    }
    pub fn expression(&self) -> &SymbolicTensor<PartialStructure> {
        &self.expression
    }
    pub fn spin_summed(&self) -> &BTreeSet<usize> {
        &self.spin_summed
    }
    pub fn color_summed(&self) -> &BTreeSet<usize> {
        &self.color_summed
    }

    /// Close selected physical spin ports. Repeated sums are errors, preventing
    /// accidental repeated averaging. References/spin vectors use unindexed momenta.
    pub fn sum_spins(
        &self,
        legs: &[usize],
        average_initial: bool,
        references: &BTreeMap<usize, Atom>,
        spin_vectors: &BTreeMap<usize, Atom>,
    ) -> Result<Self, AmplitudeError> {
        let mut result = self.clone();
        for &index in references.keys().chain(spin_vectors.keys()) {
            if !legs.contains(&index) {
                return Err(AmplitudeError::UnknownLeg(index));
            }
        }
        for &index in legs {
            let leg = self.leg(index)?;
            if !result.spin_summed.insert(index) {
                return Err(AmplitudeError::AlreadySummed {
                    kind: "spin",
                    index,
                });
            }
            let particle = self.amplitude.model().particle_by_id(leg.particle)?;
            let sum = SpinSum::new(particle, self.amplitude.model())?
                .with_dimension(&self.amplitude.dimension().to_symbolic())?
                .averaged(average_initial && leg.state == ExternalState::Incoming);
            let [ket, bra] = [
                leg.index_atom(),
                leg.tensor_index.scoped(Amplitude::bra()).to_atom(),
            ];
            let column = (leg.state == ExternalState::Incoming) != particle.is_antiparticle();
            let indices = if column { [ket, bra] } else { [bra, ket] };
            // Scalars contribute one; unsupported spins still report an error.
            let completeness = SymbolicTensor::infer(sum.expression(
                &leg.momentum(),
                indices,
                references.get(&index),
                spin_vectors.get(&index),
            )?)?;
            result.expression = result.expression.multiply(&completeness)?;
        }
        Ok(result)
    }

    /// Close selected color ports using their actual dual representations.
    pub fn sum_colors(
        &self,
        legs: &[usize],
        average_initial: bool,
    ) -> Result<Self, AmplitudeError> {
        let mut result = self.clone();
        for &index in legs {
            let leg = self.leg(index)?;
            if !result.color_summed.insert(index) {
                return Err(AmplitudeError::AlreadySummed {
                    kind: "color",
                    index,
                });
            }
            let particle = self.amplitude.model().particle_by_id(leg.particle)?;
            ColorSum::new(particle)?;
            if let Some(slot) = leg.color_slot() {
                let mut left = slot.dual();
                left.aind = leg.tensor_index;
                let mut right = slot;
                right.aind = leg.tensor_index.scoped(Amplitude::bra());
                let completeness =
                    SymbolicTensor::infer(ETS.metric(left.to_atom(), right.to_atom()))?;
                result.expression = result.expression.multiply(&completeness)?;
                if average_initial && leg.state == ExternalState::Incoming {
                    result.expression = result.expression.with_checked_expression(
                        result.expression.expression() / Atom::num(particle.color.unsigned_abs()),
                    )?;
                }
            }
        }
        Ok(result)
    }

    fn leg(&self, index: usize) -> Result<&AmplitudeLeg, AmplitudeError> {
        self.amplitude
            .legs
            .iter()
            .find(|l| l.index == index)
            .ok_or(AmplitudeError::UnknownLeg(index))
    }
}
