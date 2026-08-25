use std::{fmt, ops::RangeInclusive};

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Particle name or PDG selector accepted by process definitions.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(untagged)]
pub enum ParticleSelector {
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
            Self::Name(name) => formatter.write_str(name),
            Self::Pdg(pdg) => pdg.fmt(formatter),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GenerationType {
    Amplitude,
    CrossSection,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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

#[derive(Debug, Clone, Error, PartialEq, Eq)]
pub enum ProcessError {
    #[error("a process must specify at least one final-state alternative")]
    MissingFinalState,
    #[error("amplitude generation accepts exactly one final-state alternative")]
    MultipleAmplitudeFinalStates,
    #[error("cross-section final-state alternative {alternative} is empty")]
    EmptyCrossSectionFinalState { alternative: usize },
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
        if self.generation_type == GenerationType::CrossSection
            && let Some(alternative) = alternatives.iter().position(Vec::is_empty)
        {
            return Err(ProcessError::EmptyCrossSectionFinalState { alternative });
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
        if self.generation_type == GenerationType::CrossSection
            && let Some(alternative) = self.outgoing_alternatives.iter().position(Vec::is_empty)
        {
            return Err(ProcessError::EmptyCrossSectionFinalState { alternative });
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_empty_cross_section_alternatives_but_accepts_vacuum_amplitudes() {
        assert_eq!(
            Process::cross_section(Vec::<i64>::new(), Vec::<i64>::new()).validate(),
            Err(ProcessError::EmptyCrossSectionFinalState { alternative: 0 })
        );
        assert_eq!(
            Process::cross_section([1_i64], [1_i64])
                .with_final_state_alternatives([Vec::<i64>::new(), vec![1]]),
            Err(ProcessError::EmptyCrossSectionFinalState { alternative: 0 })
        );
        assert!(
            Process::amplitude(Vec::<i64>::new(), Vec::<i64>::new())
                .validate()
                .is_ok()
        );
    }
}
