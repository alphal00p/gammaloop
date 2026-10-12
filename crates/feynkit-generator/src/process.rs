use std::{collections::BTreeMap, fmt, str::FromStr};
#[cfg(test)]
use std::{collections::BTreeSet, sync::Arc};

use feynkit_model::{Model, ModelError, ModelFingerprint, ParticleId, VertexRuleId};
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Particle name or PDG selector accepted by process definitions.
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
#[serde(rename_all = "snake_case", tag = "kind", content = "value")]
pub enum ParticleSelector {
    Id {
        particle: ParticleId,
        model: ModelFingerprint,
    },
    Name(String),
    Pdg(i64),
}

impl From<&str> for ParticleSelector {
    fn from(value: &str) -> Self {
        Self::Name(value.to_owned())
    }
}

impl From<String> for ParticleSelector {
    fn from(value: String) -> Self {
        Self::Name(value)
    }
}

impl From<i64> for ParticleSelector {
    fn from(value: i64) -> Self {
        Self::Pdg(value)
    }
}

impl fmt::Display for ParticleSelector {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Id { particle, model } => {
                write!(formatter, "particle#{}@{model}", particle.index())
            }
            Self::Name(name) => formatter.write_str(name),
            Self::Pdg(pdg) => pdg.fmt(formatter),
        }
    }
}

impl ParticleSelector {
    pub fn by_id(model: &Model, particle: ParticleId) -> Result<Self, SelectorError> {
        model.particle_by_id(particle)?;
        Ok(Self::Id {
            particle,
            model: model.fingerprint(),
        })
    }

    pub fn resolve(&self, model: &Model) -> Result<ParticleId, SelectorError> {
        match self {
            Self::Id {
                particle,
                model: selector_model,
            } => {
                let target_model = model.fingerprint();
                if *selector_model != target_model {
                    return Err(SelectorError::ModelMismatch {
                        kind: "particle",
                        selector_model: *selector_model,
                        target_model,
                    });
                }
                model.particle_by_id(*particle)?;
                Ok(*particle)
            }
            Self::Name(name) => Ok(model.particle_id(name)?),
            Self::Pdg(pdg) => Ok(model.particle_id_by_pdg(*pdg)?),
        }
    }
}

/// Vertex-rule selector accepted at import boundaries and resolved before generation.
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
#[serde(rename_all = "snake_case", tag = "kind", content = "value")]
pub enum VertexSelector {
    Id {
        vertex: VertexRuleId,
        model: ModelFingerprint,
    },
    Name(String),
}

impl From<&str> for VertexSelector {
    fn from(value: &str) -> Self {
        Self::Name(value.to_owned())
    }
}

impl From<String> for VertexSelector {
    fn from(value: String) -> Self {
        Self::Name(value)
    }
}

impl fmt::Display for VertexSelector {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Id { vertex, model } => {
                write!(formatter, "vertex#{}@{model}", vertex.index())
            }
            Self::Name(name) => formatter.write_str(name),
        }
    }
}

impl VertexSelector {
    pub fn by_id(model: &Model, vertex: VertexRuleId) -> Result<Self, SelectorError> {
        model.vertex_rule_by_id(vertex)?;
        Ok(Self::Id {
            vertex,
            model: model.fingerprint(),
        })
    }

    pub fn resolve(&self, model: &Model) -> Result<VertexRuleId, SelectorError> {
        match self {
            Self::Id {
                vertex,
                model: selector_model,
            } => {
                let target_model = model.fingerprint();
                if *selector_model != target_model {
                    return Err(SelectorError::ModelMismatch {
                        kind: "vertex rule",
                        selector_model: *selector_model,
                        target_model,
                    });
                }
                model.vertex_rule_by_id(*vertex)?;
                Ok(*vertex)
            }
            Self::Name(name) => Ok(model.vertex_rule_id(name)?),
        }
    }
}

#[derive(Debug, Error)]
pub enum SelectorError {
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error(
        "{kind} selector belongs to model {selector_model}, but generation uses model {target_model}"
    )]
    ModelMismatch {
        kind: &'static str,
        selector_model: ModelFingerprint,
        target_model: ModelFingerprint,
    },
}

#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
#[serde(rename_all = "snake_case")]
pub enum GenerationType {
    Amplitude,
    CrossSection,
}

bincode::impl_borrow_decode!(GenerationType);

impl fmt::Display for GenerationType {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Amplitude => "Amplitude",
            Self::CrossSection => "Cross-section",
        })
    }
}

impl FromStr for GenerationType {
    type Err = ProcessError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "amplitude" => Ok(Self::Amplitude),
            "cross_section" => Ok(Self::CrossSection),
            value => Err(ProcessError::InvalidGenerationType(value.to_owned())),
        }
    }
}

#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
/// External states, independent of calculation kind, order, and query restrictions.
pub struct Process {
    incoming: Vec<ParticleSelector>,
    outgoing_alternatives: Vec<Vec<ParticleSelector>>,
}

bincode::impl_borrow_decode!(Process);

#[derive(Debug, Clone, Error, PartialEq, Eq)]
pub enum ProcessError {
    #[error("invalid generation type '{0}'")]
    InvalidGenerationType(String),
    #[error("a process must specify at least one final-state alternative")]
    MissingFinalState,
    #[error("amplitude generation accepts exactly one final-state alternative")]
    MultipleAmplitudeFinalStates,
    #[error("invalid loop range {minimum}..={maximum}")]
    InvalidLoopRange { minimum: usize, maximum: usize },
    #[error(
        "Final-state PDG {member} overlaps the covariant cut multiplet of requested physical PDG {physical}; separate physical and diagnostic requests so their observable labels remain unambiguous. For imported covariant partner graphs, supply the physical channel with --process-spec"
    )]
    CovariantCutOverlap { physical: i64, member: i64 },
    #[error(
        "Particle veto removes PDG {member} from the covariant cut multiplet of physical PDG {physical}; retain the complete vector/Goldstone/ghost sector for a physical cross section"
    )]
    CovariantCutVeto { physical: i64, member: i64 },
}

impl Process {
    /// Request unsewn compact diagrams with an explicit total loop range.
    pub fn amplitude(self, loops: crate::Loops) -> crate::Request {
        crate::Request::amplitude(self, loops)
    }

    /// Request sewn forward compact diagrams with an explicit forward loop range.
    pub fn cross_section(self, loops: crate::Loops) -> crate::Request {
        crate::Request::cross_section(self, loops)
    }

    pub fn new<I, O, PI, PO>(incoming: I, outgoing: O) -> Self
    where
        I: IntoIterator<Item = PI>,
        O: IntoIterator<Item = PO>,
        PI: Into<ParticleSelector>,
        PO: Into<ParticleSelector>,
    {
        Self {
            incoming: incoming.into_iter().map(Into::into).collect(),
            outgoing_alternatives: vec![outgoing.into_iter().map(Into::into).collect()],
        }
    }

    pub fn with_final_state_alternatives<I, O, P>(
        mut self,
        alternatives: I,
    ) -> Result<Self, ProcessError>
    where
        I: IntoIterator<Item = O>,
        O: IntoIterator<Item = P>,
        P: Into<ParticleSelector>,
    {
        let alternatives: Vec<_> = alternatives
            .into_iter()
            .map(|outgoing| outgoing.into_iter().map(Into::into).collect())
            .collect();
        if alternatives.is_empty() {
            return Err(ProcessError::MissingFinalState);
        }
        self.outgoing_alternatives = alternatives;
        Ok(self)
    }

    /// Validate every external selector against one concrete model.
    pub fn validate_in(&self, model: &Model) -> Result<(), crate::GenerationError> {
        self.validate()?;
        for selector in self
            .incoming
            .iter()
            .chain(self.outgoing_alternatives.iter().flatten())
        {
            selector.resolve(model)?;
        }
        Ok(())
    }

    pub fn incoming(&self) -> &[ParticleSelector] {
        &self.incoming
    }

    pub fn outgoing_alternatives(&self) -> &[Vec<ParticleSelector>] {
        &self.outgoing_alternatives
    }

    /// Resolve explicit final-state alternatives to PDG codes.
    pub fn outgoing_pdgs(&self, model: &Model) -> Result<Vec<Vec<i64>>, crate::GenerationError> {
        self.outgoing_alternatives
            .iter()
            .map(|state| {
                state
                    .iter()
                    .map(|selector| Ok(model.particle_by_id(selector.resolve(model)?)?.pdg_code))
                    .collect()
            })
            .collect()
    }

    /// Observable labels for the complete covariant partner sets requested by
    /// physical vector states. Explicit diagnostic states retain their identity.
    pub fn covariant_cut_representatives(
        &self,
        model: &Model,
        generation_type: GenerationType,
    ) -> Result<BTreeMap<i64, i64>, crate::GenerationError> {
        if generation_type != GenerationType::CrossSection {
            return Ok(BTreeMap::new());
        }
        let states = self.outgoing_pdgs(model)?;
        Ok(model
            .covariant_cut_multiplets()
            .iter()
            .filter(|(physical, _)| states.iter().any(|state| state.contains(physical)))
            .flat_map(|(&physical, members)| members.iter().map(move |&member| (member, physical)))
            .collect())
    }

    /// Restrictions on a physical cross section must retain each declared covariant partner.
    pub fn validate_covariant_cut_restrictions(
        &self,
        model: &Model,
        restrictions: &crate::Restrict,
        kind: GenerationType,
    ) -> Result<(), crate::GenerationError> {
        let representatives = self.covariant_cut_representatives(model, kind)?;
        let restrictions = restrictions.resolve(model)?;
        for (member, physical) in representatives {
            let id = model.particle_id_by_pdg(member)?;
            if !restrictions.particle(id) {
                return Err(ProcessError::CovariantCutVeto { physical, member }.into());
            }
        }
        Ok(())
    }

    pub fn covariant_cut_states(
        &self,
        model: &Model,
        generation_type: GenerationType,
    ) -> Result<Vec<Vec<i64>>, crate::GenerationError> {
        let requested = self.outgoing_pdgs(model)?;
        if generation_type != GenerationType::CrossSection {
            return Ok(requested);
        }
        let representatives = self.covariant_cut_representatives(model, generation_type)?;
        for &member in requested.iter().flatten() {
            if let Some(&physical) = representatives.get(&member)
                && member != physical
            {
                return Err(ProcessError::CovariantCutOverlap { physical, member }.into());
            }
        }
        Ok(model.covariant_cut_states(&requested)?)
    }

    pub fn validate(&self) -> Result<(), ProcessError> {
        if self.outgoing_alternatives.is_empty() {
            return Err(ProcessError::MissingFinalState);
        }
        Ok(())
    }
}

/// Construct validated process definitions without coupling the model crate to generation.
pub trait ModelProcessExt {
    fn process<I, O, PI, PO>(
        &self,
        incoming: I,
        outgoing: O,
    ) -> Result<Process, crate::GenerationError>
    where
        I: IntoIterator<Item = PI>,
        O: IntoIterator<Item = PO>,
        PI: Into<ParticleSelector>,
        PO: Into<ParticleSelector>;
}

impl ModelProcessExt for Model {
    fn process<I, O, PI, PO>(
        &self,
        incoming: I,
        outgoing: O,
    ) -> Result<Process, crate::GenerationError>
    where
        I: IntoIterator<Item = PI>,
        O: IntoIterator<Item = PO>,
        PI: Into<ParticleSelector>,
        PO: Into<ParticleSelector>,
    {
        let process = Process::new(incoming, outgoing);
        process.validate_in(self)?;
        Ok(process)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn particle_selection_is_the_complementary_species_veto() {
        use crate::{BuiltDiagramSet, Loops, NumeratorOptions, Restrict, Run, TopologyRules};
        let model = Arc::new(
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap(),
        );
        let process = model.process(["e-", "e+"], ["a", "a"]).unwrap();
        let request = process
            .amplitude(Loops::total(0..=0))
            .with_topology(TopologyRules::new().with_self_loops(false));
        let selected = request.clone().with_restrict(
            Restrict::new().only_particles([ParticleSelector::from("e-"), 22_i64.into()]),
        );
        selected.validate_in(&model).unwrap();
        let allowed = selected.restrictions().resolve(&model).unwrap();
        let vetoes = model
            .particles()
            .iter()
            .map(|particle| model.particle_id(&particle.name).unwrap())
            .filter(|&particle| !allowed.particle(particle))
            .map(|particle| ParticleSelector::by_id(&model, particle).unwrap())
            .collect::<Vec<_>>();
        let vetoed = request
            .clone()
            .with_restrict(Restrict::new().exclude_particles(vetoes));
        let run = Run::default();
        let generated = selected
            .generate(model.clone(), &run)
            .unwrap()
            .build(&NumeratorOptions::default(), &run)
            .unwrap();
        assert_eq!(generated.len(), 2);
        let ids = |result: BuiltDiagramSet| {
            result
                .iter()
                .map(|diagram| diagram.id().unwrap())
                .collect::<BTreeSet<_>>()
        };
        assert_eq!(
            ids(generated),
            ids(vetoed
                .generate(model.clone(), &run)
                .unwrap()
                .build(&NumeratorOptions::default(), &run)
                .unwrap())
        );
        let restored: crate::Request =
            serde_json::from_str(&serde_json::to_string(&selected).unwrap()).unwrap();
        assert_eq!(restored, selected);
        // Restriction clauses now intersect: an explicit empty allow-list permits none,
        // including when it is combined with an empty exclusion list.
        let empty = request.clone().with_restrict(
            Restrict::new()
                .only_particles(Vec::<ParticleSelector>::new())
                .exclude_particles(Vec::<ParticleSelector>::new()),
        );
        assert!(empty.generate(model.clone(), &run).unwrap().is_empty());
        let empty = request
            .clone()
            .with_restrict(Restrict::new().only_particles(Vec::<ParticleSelector>::new()));
        assert!(empty.generate(model.clone(), &run).unwrap().is_empty());
        assert_eq!(request.generate(model.clone(), &run).unwrap().len(), 2);
        assert!(
            request
                .with_restrict(Restrict::new().only_particles(["missing"]))
                .validate_in(&model)
                .is_err()
        );
    }

    #[test]
    fn model_process_validates_references_and_round_trips_only_process_state() {
        use crate::{Loops, Restrict};
        let model =
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
        assert!(model.process(["missing"], ["a"]).is_err());
        let process = model.process(["e-", "e+"], ["a", "a"]).unwrap();
        process.validate_in(&model).unwrap();
        let request = process
            .clone()
            .amplitude(Loops::total(0..=0))
            .with_restrict(
                Restrict::new()
                    .exclude_particles(["t"])
                    .only_vertices(["V_98"]),
            );
        request.validate_in(&model).unwrap();
        let definition = serde_json::to_value(&process).unwrap();
        assert!(definition.get("loop_count").is_none());
        assert!(definition.get("symmetrize_final").is_none());
        assert_eq!(
            serde_json::from_value::<Process>(definition).unwrap(),
            process
        );
        assert!(
            process
                .amplitude(Loops::total(0..=0))
                .with_restrict(Restrict::new().only_vertices(["missing"]))
                .validate_in(&model)
                .is_err()
        );
    }

    #[test]
    fn process_veto_cannot_remove_a_physical_covariant_partner() {
        use crate::{Loops, Restrict};
        let model =
            Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
        let process = model.process(["H"], ["W+", "W-"]).unwrap();
        process.validate_in(&model).unwrap();
        let request = process
            .cross_section(Loops::forward(1..=1))
            .with_restrict(Restrict::new().exclude_particles(["ghWp"]));
        let error = request.validate_in(&model).unwrap_err();
        assert!(matches!(
            error,
            crate::GenerationError::Process(ProcessError::CovariantCutVeto { .. })
        ));
    }

    #[test]
    fn accepts_vacuum_processes_and_empty_cross_section_alternatives() {
        assert!(
            Process::new(Vec::<i64>::new(), Vec::<i64>::new())
                .validate()
                .is_ok()
        );
        assert!(
            Process::new([1_i64], [1_i64])
                .with_final_state_alternatives([Vec::<i64>::new(), vec![1]])
                .is_ok()
        );
        assert!(
            Process::new(Vec::<i64>::new(), Vec::<i64>::new())
                .validate()
                .is_ok()
        );
    }
}
