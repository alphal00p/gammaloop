use std::{io, path::PathBuf};

use thiserror::Error;

/// A named collection in a physics model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EntityKind {
    Order,
    Parameter,
    Particle,
    Propagator,
    LorentzStructure,
    Coupling,
    VertexRule,
    ModelFunction,
    FormFactor,
}

impl std::fmt::Display for EntityKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Order => "order",
            Self::Parameter => "parameter",
            Self::Particle => "particle",
            Self::Propagator => "propagator",
            Self::LorentzStructure => "Lorentz structure",
            Self::Coupling => "coupling",
            Self::VertexRule => "vertex rule",
            Self::ModelFunction => "model function",
            Self::FormFactor => "form factor",
        })
    }
}

/// A structural error in a model definition.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum ModelValidationError {
    #[error("model name must not be empty")]
    EmptyModelName,
    #[error("duplicate {kind} name '{name}'")]
    DuplicateName { kind: EntityKind, name: String },
    #[error("duplicate particle PDG code {pdg}")]
    DuplicatePdg { pdg: i64 },
    #[error("{kind} '{name}' references unknown {target_kind} '{target}' through '{field}'")]
    UnknownReference {
        kind: EntityKind,
        name: String,
        field: &'static str,
        target_kind: EntityKind,
        target: String,
    },
    #[error(
        "particle '{particle}' names '{antiparticle}' as its antiparticle, but '{antiparticle}' names '{actual_antiparticle}'"
    )]
    NonReciprocalAntiparticle {
        particle: String,
        antiparticle: String,
        actual_antiparticle: String,
    },
    #[error(
        "particle '{particle}' references propagator '{propagator}', which belongs to particle '{propagator_particle}'"
    )]
    PropagatorParticleMismatch {
        particle: String,
        propagator: String,
        propagator_particle: String,
    },
    #[error(
        "Lorentz structure '{lorentz_structure}' has spin signature {lorentz_spins:?}, but vertex rule '{vertex}' requires {particle_spins:?}"
    )]
    LorentzSpinMismatch {
        vertex: String,
        lorentz_structure: String,
        particle_spins: Vec<i64>,
        lorentz_spins: Vec<i64>,
    },
    #[error("external parameter '{name}' has no value")]
    ExternalParameterWithoutValue { name: String },
    #[error("internal parameter '{name}' has neither a value nor an expression")]
    UnresolvedInternalParameter { name: String },
    #[error("vertex rule '{name}' has {actual} coupling rows, but {expected} color structures")]
    CouplingRowCount {
        name: String,
        expected: usize,
        actual: usize,
    },
    #[error(
        "coupling row {row} of vertex rule '{name}' has {actual} entries, but {expected} Lorentz structures"
    )]
    CouplingColumnCount {
        name: String,
        row: usize,
        expected: usize,
        actual: usize,
    },
}

/// Errors returned by model and parameter-card APIs.
#[derive(Debug, Error)]
pub enum ModelError {
    #[error("could not read model data from '{}': {source}", path.display())]
    Read {
        path: PathBuf,
        #[source]
        source: io::Error,
    },
    #[error("could not write model data to '{}': {source}", path.display())]
    Write {
        path: PathBuf,
        #[source]
        source: io::Error,
    },
    #[error("invalid JSON: {0}")]
    Json(#[from] serde_json::Error),
    #[error(transparent)]
    Validation(#[from] ModelValidationError),
    #[error("{kind} '{key}' was not found in model '{model}'")]
    NotFound {
        model: String,
        kind: EntityKind,
        key: String,
    },
    #[error("parameter card references unknown parameter '{name}'")]
    UnknownCardParameter { name: String },
    #[error("external parameter '{name}' has no value for the default parameter card")]
    MissingCardValue { name: String },
}
