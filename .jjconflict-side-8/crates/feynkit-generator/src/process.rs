use std::{fmt, ops::RangeInclusive, str::FromStr};

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
pub struct Process {
    generation_type: GenerationType,
    incoming: Vec<ParticleSelector>,
    outgoing_alternatives: Vec<Vec<ParticleSelector>>,
    loop_count: RangeInclusive<usize>,
    symmetrize_initial: bool,
    symmetrize_final: bool,
    symmetrize_left_right: bool,
    symmetrize_external_fermions: bool,
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
}

impl Process {
    pub fn amplitude<I, O, PI, PO>(incoming: I, outgoing: O) -> Self
    where
        I: IntoIterator<Item = PI>,
        O: IntoIterator<Item = PO>,
        PI: Into<ParticleSelector>,
        PO: Into<ParticleSelector>,
    {
        Self {
            generation_type: GenerationType::Amplitude,
            incoming: incoming.into_iter().map(Into::into).collect(),
            outgoing_alternatives: vec![outgoing.into_iter().map(Into::into).collect()],
            loop_count: 0..=0,
            symmetrize_initial: false,
            symmetrize_final: false,
            symmetrize_left_right: false,
            symmetrize_external_fermions: false,
        }
    }

    pub fn cross_section<I, O, PI, PO>(incoming: I, outgoing: O) -> Self
    where
        I: IntoIterator<Item = PI>,
        O: IntoIterator<Item = PO>,
        PI: Into<ParticleSelector>,
        PO: Into<ParticleSelector>,
    {
        Self {
            generation_type: GenerationType::CrossSection,
            incoming: incoming.into_iter().map(Into::into).collect(),
            outgoing_alternatives: vec![outgoing.into_iter().map(Into::into).collect()],
            loop_count: 0..=0,
            symmetrize_initial: false,
            symmetrize_final: false,
            symmetrize_left_right: false,
            symmetrize_external_fermions: false,
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
        if self.generation_type == GenerationType::Amplitude && alternatives.len() != 1 {
            return Err(ProcessError::MultipleAmplitudeFinalStates);
        }
        self.outgoing_alternatives = alternatives;
        Ok(self)
    }

    pub fn with_loop_count(mut self, minimum: usize, maximum: usize) -> Result<Self, ProcessError> {
        if minimum > maximum {
            return Err(ProcessError::InvalidLoopRange { minimum, maximum });
        }
        self.loop_count = minimum..=maximum;
        Ok(self)
    }

    pub fn symmetrize_initial(mut self, enabled: bool) -> Self {
        self.symmetrize_initial = enabled;
        self
    }

    pub fn symmetrize_final(mut self, enabled: bool) -> Self {
        self.symmetrize_final = enabled;
        self
    }

    pub fn symmetrize_left_right(mut self, enabled: bool) -> Self {
        self.symmetrize_left_right = enabled;
        self
    }

    /// Include amplitude fermions in enabled external-state symmetry classes.
    ///
    /// This is disabled by default. Cross-section symmetry is unchanged by
    /// this amplitude-specific policy switch.
    pub fn symmetrize_external_fermions(mut self, enabled: bool) -> Self {
        self.symmetrize_external_fermions = enabled;
        self
    }

    pub fn generation_type(&self) -> GenerationType {
        self.generation_type
    }

    pub fn incoming(&self) -> &[ParticleSelector] {
        &self.incoming
    }

    pub fn outgoing_alternatives(&self) -> &[Vec<ParticleSelector>] {
        &self.outgoing_alternatives
    }

    pub fn loop_count(&self) -> RangeInclusive<usize> {
        self.loop_count.clone()
    }

    pub fn symmetrizes_initial(&self) -> bool {
        self.symmetrize_initial
    }

    pub fn symmetrizes_final(&self) -> bool {
        self.symmetrize_final
    }

    pub fn symmetrizes_left_right(&self) -> bool {
        self.symmetrize_left_right
    }

    pub fn symmetrizes_external_fermions(&self) -> bool {
        self.symmetrize_external_fermions
    }

    pub fn validate(&self) -> Result<(), ProcessError> {
        if self.outgoing_alternatives.is_empty() {
            return Err(ProcessError::MissingFinalState);
        }
        if self.generation_type == GenerationType::Amplitude
            && self.outgoing_alternatives.len() != 1
        {
            return Err(ProcessError::MultipleAmplitudeFinalStates);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accepts_vacuum_processes_and_empty_cross_section_alternatives() {
        assert!(
            Process::cross_section(Vec::<i64>::new(), Vec::<i64>::new())
                .validate()
                .is_ok()
        );
        assert!(
            Process::cross_section([1_i64], [1_i64])
                .with_final_state_alternatives([Vec::<i64>::new(), vec![1]])
                .is_ok()
        );
        assert!(
            Process::amplitude(Vec::<i64>::new(), Vec::<i64>::new())
                .validate()
                .is_ok()
        );
    }
}
