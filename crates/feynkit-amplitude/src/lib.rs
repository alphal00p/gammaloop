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
    IndexTooling, dirac::GammaSimplifySettings, representations::Bispinor, tensor::SymbolicTensor,
};
use linnet::half_edge::involution::HedgePair;
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        abstract_index::AbstractIndex,
        dimension::Dimension,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        representation::{LibrarySlot, Minkowski, RepName},
        slot::{DualSlotTo, IsAbstractSlot, ParseableAind},
    },
};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    function,
    id::Replacement,
    symbol,
};
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
    #[error("invalid squared-amplitude pairing: {0}")]
    InvalidPairing(String),
    #[error(transparent)]
    Diagram(#[from] DiagramError),
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error(transparent)]
    Spin(#[from] SpinSumError),
    #[error(transparent)]
    Color(#[from] ColorSumError),
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

/// Resolved class relation for a constrained amplitude square.
///
/// Classes identify already-selected, weighted terms; they do not infer orders from symbolic
/// expressions. Ordered pairs include both orientations of every off-diagonal interference.
/// Class identifiers remain stable when selection removes every term of a class: that class
/// contributes zero, and its surviving partner does not become unrestricted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SquaredPairing {
    class_count: usize,
    term_classes: Vec<usize>,
    allowed: BTreeSet<(usize, usize)>,
}

impl SquaredPairing {
    /// Validate class references and symmetry. Duplicate pairs denote the same relation entry
    /// and are counted once. The amplitude constructor separately checks the number of terms.
    pub fn new(
        class_count: usize,
        term_classes: impl IntoIterator<Item = usize>,
        allowed_ordered_pairs: impl IntoIterator<Item = (usize, usize)>,
    ) -> Result<Self, AmplitudeError> {
        let term_classes: Vec<_> = term_classes.into_iter().collect();
        if let Some(&class) = term_classes.iter().find(|&&class| class >= class_count) {
            return Err(AmplitudeError::InvalidPairing(format!(
                "term class {class} is outside {class_count} classes"
            )));
        }
        let allowed: BTreeSet<_> = allowed_ordered_pairs.into_iter().collect();
        for &(left, right) in &allowed {
            if left >= class_count || right >= class_count {
                return Err(AmplitudeError::InvalidPairing(format!(
                    "pair ({left}, {right}) is outside {class_count} classes"
                )));
            }
            if !allowed.contains(&(right, left)) {
                return Err(AmplitudeError::InvalidPairing(format!(
                    "pair ({left}, {right}) has no reverse pair"
                )));
            }
        }
        Ok(Self {
            class_count,
            term_classes,
            allowed,
        })
    }

    pub fn class_count(&self) -> usize {
        self.class_count
    }

    pub fn term_classes(&self) -> &[usize] {
        &self.term_classes
    }

    pub fn allowed_pairs(&self) -> &BTreeSet<(usize, usize)> {
        &self.allowed
    }

    fn class_sums(&self, terms: &[Atom]) -> BTreeMap<usize, Atom> {
        let mut classes = BTreeMap::<usize, Vec<&Atom>>::new();
        for (&class, term) in self.term_classes.iter().zip(terms) {
            classes.entry(class).or_default().push(term);
        }
        classes
            .into_iter()
            .map(|(class, terms)| (class, Atom::add_many(terms)))
            .collect()
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
    terms: Vec<Atom>,
    legs: Vec<AmplitudeLeg>,
    options: AmplitudeOptions,
    conjugated: bool,
    pairing: Option<SquaredPairing>,
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
                return Err(AmplitudeError::DifferentModel(
                    diagram.name().unwrap_or("<unnamed>").to_owned(),
                ));
            }
            let (term, external) = Self::diagram_term(diagram, &options)?;
            if terms.is_empty() {
                legs = external;
            } else if legs != external {
                return Err(AmplitudeError::DifferentExternals(
                    diagram.name().unwrap_or("<unnamed>").to_owned(),
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
            pairing: None,
        })
    }

    /// Attach a resolved collection policy. Raw-list constructors remain unrestricted.
    /// Terms and their classes must have the same order after selection or grouping.
    pub fn with_pairing(mut self, pairing: SquaredPairing) -> Result<Self, AmplitudeError> {
        if pairing.term_classes.len() != self.terms.len() {
            return Err(AmplitudeError::InvalidPairing(format!(
                "{} term classes for {} amplitude terms",
                pairing.term_classes.len(),
                self.terms.len()
            )));
        }
        self.pairing = Some(pairing);
        Ok(self)
    }

    pub fn pairing(&self) -> Option<&SquaredPairing> {
        self.pairing.as_ref()
    }

    pub fn diagrams(&self) -> &[Arc<FeynmanDiagram>] {
        &self.diagrams
    }
    pub fn terms(&self) -> &[Atom] {
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

    /// The coherent retained amplitude. A transferred pairing constrains only its square.
    pub fn expression(&self) -> Atom {
        Atom::add_many(&self.terms)
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
                let adjoint = term
                    .dirac_adjoint::<AbstractIndex>(true)
                    .map_err(|e| AmplitudeError::Tensor(e.to_string()))?;
                let adjoint = self.apply_reality(&adjoint);
                // Physical leg labels and their declared order belong to the
                // amplitude. The shared scheduler owns gamma/metric cleanup.
                SymbolicTensor::checked_parts(adjoint, interface.clone())
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
                    .map(SymbolicTensor::into_expression)
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
        let bra = |expression: Atom| {
            expression
                .wrap_indices(Self::bra())
                .replace(function!(symbols::loop_momentum(), loop_index, args))
                .with(function!(
                    symbols::loop_momentum(),
                    Self::bra().call(loop_index),
                    args
                ))
        };
        let expression = if let Some(pairing) = &amplitude.pairing {
            // Terms already include their overall factors. Sum each class before multiplying,
            // and retain both ordered cross terms rather than assuming a real interference.
            let ket = pairing.class_sums(&amplitude.terms);
            let bra: BTreeMap<_, _> = pairing
                .class_sums(&adjoint.terms)
                .into_iter()
                .map(|(class, expression)| (class, bra(expression)))
                .collect();
            Atom::add_many(
                pairing
                    .allowed
                    .iter()
                    .filter_map(|(left, right)| Some(ket.get(left)? * bra.get(right)?)),
            )
        } else {
            amplitude.expression() * bra(adjoint.expression())
        };
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

    fn apply_reality(&self, expression: &Atom) -> Atom {
        let conjugate = symbol!("spenso::conj");
        let mut rules = Vec::new();
        for parameter in self
            .model()
            .parameters()
            .iter()
            .filter(|p| p.parameter_type == ParameterType::Real)
        {
            let real = Atom::var(symbol!(&format!("UFO::{}", parameter.name)));
            rules.push(Replacement::new(conjugate.call(&real).to_pattern(), real));
        }
        for real in &self.options.real {
            rules.push(Replacement::new(
                conjugate.call(real).to_pattern(),
                real.clone(),
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
            conjugate.call(&dimension).to_pattern(),
            dimension,
        ));
        expression.replace_multiple(rules)
    }

    fn diagram_term(
        diagram: &FeynmanDiagram,
        options: &AmplitudeOptions,
    ) -> Result<(Atom, Vec<AmplitudeLeg>), AmplitudeError> {
        if !diagram.cuts().is_empty() {
            return Err(AmplitudeError::SewnDiagram(
                diagram.name().unwrap_or("<unnamed>").to_owned(),
            ));
        }
        let dimension = options.dimension.to_symbolic();
        let index = symbol!("feynkit_amplitude::index_");
        let numerator = (diagram.numerator() * diagram.numerator_prefactor())
            .replace(Minkowski {}.new_rep(4).pattern(Atom::var(index)))
            .with(Minkowski {}.to_symbolic([dimension, Atom::var(index)]));
        let ports = numerator
            .list_dangling::<AbstractIndex>()
            .map_err(|e| AmplitudeError::Tensor(e.to_string()))?
            .iter()
            .map(|a| {
                LibrarySlot::<AbstractIndex>::try_from(a.as_view())
                    .map_err(|e| AmplitudeError::Tensor(e.to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut legs = BTreeMap::new();
        let mut identified = 0;
        for (pair, _, edge) in diagram.underlying().iter_edges() {
            let Some(external) = &edge.data.external else {
                continue;
            };
            let HedgePair::Unpaired { hedge, .. } = pair else {
                return Err(AmplitudeError::SewnDiagram(
                    diagram.name().unwrap_or("<unnamed>").to_owned(),
                ));
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
            for slot in ports.iter().filter(|s| s.aind.to_atom() == old_index) {
                leg.slots.push(*slot);
                identified += 1;
            }
            leg.slots.sort();
            let particle = diagram.model().particle_by_id(leg.particle)?;
            let expected = usize::from(particle.spin != 1) + usize::from(particle.color != 1);
            if leg.slots.len() != expected || legs.insert(leg.index, leg).is_some() {
                return Err(AmplitudeError::ExternalPorts(
                    diagram.name().unwrap_or("<unnamed>").to_owned(),
                    "missing, duplicate, or unsupported external slots".into(),
                ));
            }
        }
        if identified != ports.len() {
            return Err(AmplitudeError::ExternalPorts(
                diagram.name().unwrap_or("<unnamed>").to_owned(),
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
            .expand_couplings(&(factor * numerator / denominator));
        let expression = basis
            .route_expression(&expression)
            .replace_multiple(momentum_rules);
        Ok((expression, legs.into_values().collect()))
    }
}

/// The resolved allowed interference of an amplitude and its adjoint, with independent ports.
/// Without a pairing policy this is the full coherent product.
#[derive(Clone, Debug)]
pub struct SquaredAmplitude {
    amplitude: Amplitude,
    expression: Atom,
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
    pub fn expression(&self) -> &Atom {
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
            result.expression *= sum.expression(
                &leg.momentum(),
                indices,
                references.get(&index),
                spin_vectors.get(&index),
            )?;
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
                result.expression *= ETS.metric(left.to_atom(), right.to_atom());
                if average_initial && leg.state == ExternalState::Incoming {
                    result.expression /= Atom::num(particle.color.unsigned_abs());
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

#[cfg(test)]
mod pairing_tests {
    use super::*;
    use symbolica::parse;

    // Reuse the scalar DOT boundary from the existing generator amplitude tests. Each source
    // numerator stays factorized; only derived scalar squares are expanded for exact comparison.
    fn scalars(terms: &[(&str, i64)]) -> Amplitude {
        let model = Arc::new(
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap(),
        );
        Amplitude::from_diagrams(terms.iter().map(|(term, factor)| {
            Arc::new(
                FeynmanDiagram::from_dot(
                    Arc::clone(&model),
                    &format!(
                        r#"digraph {{
                        a [num="amplitude_pairing_test::{term}"];
                        ext [style=invis];
                        ext -> a [particle="H"];
                        a -> ext [particle="H"];
                        a -> ext [particle="H"];
                    }}"#
                    ),
                )
                .unwrap()
                .with_overall_factor(Atom::num(*factor)),
            )
        }))
        .unwrap()
    }

    #[test]
    fn pairing_validates_references_symmetry_and_term_count() {
        for invalid in [
            SquaredPairing::new(2, [2], []),
            SquaredPairing::new(2, [0], [(0, 2), (2, 0)]),
            SquaredPairing::new(2, [0, 1], [(0, 1)]),
            SquaredPairing::new(0, [0], []),
        ] {
            assert!(matches!(invalid, Err(AmplitudeError::InvalidPairing(_))));
        }
        let pairing = SquaredPairing::new(2, [0, 1], [(0, 1), (1, 0), (0, 1)]).unwrap();
        assert_eq!(pairing.class_count(), 2);
        assert_eq!(pairing.term_classes(), &[0, 1]);
        assert_eq!(pairing.allowed_pairs().len(), 2);
        assert!(matches!(
            scalars(&[("a", 1)]).with_pairing(pairing),
            Err(AmplitudeError::InvalidPairing(_))
        ));
        assert!(SquaredPairing::new(0, [], []).is_ok());
        // The existing empty-amplitude error remains independent of the pairing relation.
        assert!(matches!(
            Amplitude::from_diagrams([]),
            Err(AmplitudeError::Empty)
        ));
    }

    #[test]
    fn constrained_square_preserves_oriented_complex_interference_and_class_sums() {
        // The two a occurrences and b share one class; c is its partner. Two d occurrences
        // share a diagonal class. Every overall factor and ordered interference occurs once.
        let amplitude = scalars(&[("a", 2), ("a", 5), ("b", 3), ("c", 7), ("d", 11), ("d", 2)])
            .with_pairing(
                SquaredPairing::new(4, [1, 1, 1, 3, 2, 2], [(1, 3), (3, 1), (2, 2)]).unwrap(),
            )
            .unwrap();
        let expected = parse!(
            "49*amplitude_pairing_test::a*spenso::conj(amplitude_pairing_test::c)\
             +21*amplitude_pairing_test::b*spenso::conj(amplitude_pairing_test::c)\
             +49*amplitude_pairing_test::c*spenso::conj(amplitude_pairing_test::a)\
             +21*amplitude_pairing_test::c*spenso::conj(amplitude_pairing_test::b)\
             +169*amplitude_pairing_test::d*spenso::conj(amplitude_pairing_test::d)"
        );
        let square = amplitude.squared().unwrap();
        assert!(
            (square.expression() - &expected).expand().is_zero(),
            "every scalar coefficient, including forbidden diagonals, must agree exactly"
        );
        assert_eq!(square.amplitude().pairing(), amplitude.pairing());
        let coherent = parse!(
            "7*amplitude_pairing_test::a+3*amplitude_pairing_test::b\
             +7*amplitude_pairing_test::c+13*amplitude_pairing_test::d"
        );
        assert_eq!(amplitude.expression(), coherent);
        let adjoint = amplitude.conjugate().unwrap();
        assert_eq!(adjoint.pairing(), amplitude.pairing());
        assert_eq!(adjoint.clone().pairing(), amplitude.pairing());
        assert_eq!(adjoint.conjugate().unwrap().expression(), coherent);
        assert!(
            (adjoint.squared().unwrap().expression() - &expected)
                .expand()
                .is_zero()
        );
    }

    #[test]
    fn a_removed_partner_contributes_zero_without_restoring_unrestricted_pairing() {
        let relation = [(1, 3), (3, 1)];
        let amplitude = scalars(&[("a", 2)])
            .with_pairing(SquaredPairing::new(4, [1], relation).unwrap())
            .unwrap();
        assert!(!amplitude.expression().is_zero());
        assert!(amplitude.squared().unwrap().expression().is_zero());
        assert_eq!(amplitude.pairing().unwrap().allowed_pairs().len(), 2);
        let empty_relation = scalars(&[("a", 2)])
            .with_pairing(SquaredPairing::new(1, [0], []).unwrap())
            .unwrap();
        assert!(empty_relation.squared().unwrap().expression().is_zero());
    }

    #[test]
    fn raw_lists_remain_unrestricted_and_keep_full_interference() {
        let amplitude = scalars(&[("a", 2), ("b", 3)]);
        assert!(amplitude.pairing().is_none());
        let expected = parse!(
            "(2*amplitude_pairing_test::a+3*amplitude_pairing_test::b)\
             *(2*spenso::conj(amplitude_pairing_test::a)\
             +3*spenso::conj(amplitude_pairing_test::b))"
        );
        assert!(
            (amplitude.squared().unwrap().expression() - &expected)
                .expand()
                .is_zero()
        );
        let constrained = amplitude
            .with_pairing(SquaredPairing::new(2, [0, 1], [(0, 0), (0, 1), (1, 0), (1, 1)]).unwrap())
            .unwrap();
        assert!(
            (constrained.squared().unwrap().expression() - &expected)
                .expand()
                .is_zero()
        );
    }

    #[test]
    fn constrained_bra_uses_the_existing_independent_loop_scope() {
        let model = Arc::new(
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap(),
        );
        let diagram = FeynmanDiagram::from_dot(
            model,
            r#"digraph triangle {
            ext [style=invis];
            ext -> a [id=0, particle="H"];
            b -> ext [id=1, particle="H"];
            c -> ext [id=2, particle="H"];
            a -> b [id=3, particle="H", lmb_id=0];
            b -> c [id=4, particle="H"];
            c -> a [id=5, particle="H"];
        }"#,
        )
        .unwrap();
        let amplitude = Amplitude::from_diagram(diagram).unwrap();
        let unrestricted = amplitude.squared().unwrap();
        let constrained = amplitude
            .with_pairing(SquaredPairing::new(1, [0], [(0, 0)]).unwrap())
            .unwrap()
            .squared()
            .unwrap();
        assert_eq!(constrained.expression(), unrestricted.expression());
        assert!(
            constrained
                .expression()
                .contains(parse!("gammalooprs::K(0,spenso::mink(4))"))
        );
        assert!(constrained.expression().contains(parse!(
            "gammalooprs::K(feynkit_amplitude::bra(0),spenso::mink(4))"
        )));
    }
}
