pub mod diagram_generator;

use ahash::HashMap;
use std::{fmt, str::FromStr};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FeynGenError {
    #[error("{0}")]
    GenericError(String),
}

#[derive(Debug, Clone)]
pub enum GenerationType {
    Amplitude,
    CrossSection,
}

impl fmt::Display for GenerationType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Amplitude => "Amplitude",
                Self::CrossSection => "Cross-section",
            }
        )
    }
}

impl FromStr for GenerationType {
    type Err = FeynGenError;

    fn from_str(s: &str) -> Result<Self, FeynGenError> {
        match s {
            "amplitude" => Ok(Self::Amplitude),
            "cross_section" => Ok(Self::CrossSection),
            _ => Err(FeynGenError::GenericError(format!(
                "Invalid generation type: {}",
                s
            ))),
        }
    }
}

#[derive(Debug, Clone)]
pub struct FeynGenFilters(pub Vec<FeynGenFilter>);

impl FeynGenFilters {
    pub fn get_max_bridge(&self) -> Option<usize> {
        self.0.iter().find_map(|f| match f {
            FeynGenFilter::MaxNumberOfBridges(n) => Some(*n),
            _ => None,
        })
    }

    pub fn allow_tadpoles(&self) -> bool {
        !self
            .0
            .iter()
            .any(|f| matches!(f, FeynGenFilter::NoTadpoles))
    }
}

#[derive(Debug, Clone)]
pub enum FeynGenFilter {
    No1PI,
    ParticleVeto(Vec<i64>),
    MaxNumberOfBridges(usize),
    NoTadpoles,
    CouplingOrders(HashMap<String, usize>),
}

impl fmt::Display for FeynGenFilter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::No1PI => "No1PI".into(),
                Self::ParticleVeto(pdgs) => format!(
                    "ParticleVeto({})",
                    pdgs.iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<String>>()
                        .join("|")
                ),
                Self::MaxNumberOfBridges(n) => format!("MaxNumberOfBridges({})", n),
                Self::NoTadpoles => "NoTadpoles".into(),
                Self::CouplingOrders(orders) => format!(
                    "CouplingOrders({})",
                    orders
                        .iter()
                        .map(|(k, v)| format!("{}={}", k, v))
                        .collect::<Vec<String>>()
                        .join("|")
                ),
            }
        )
    }
}

#[derive(Debug, Clone)]
pub struct FeynGenOptions {
    pub generation_type: GenerationType,
    pub initial_pdgs: Vec<i64>,
    pub final_pdgs: Vec<i64>,
    pub loop_count_range: (usize, usize),
    pub filters: FeynGenFilters,
}

impl fmt::Display for FeynGenOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Generation type: {}\nInitial PDGs: {:?}\nFinal PDGs: {:?}\nLoop count: {}\nFilters:{}{}",
            self.generation_type,
            self.initial_pdgs,
            self.final_pdgs,
            if self.loop_count_range.0 == self.loop_count_range.1 {
                format!("{}", self.loop_count_range.0)
            } else {
                format!("{:?}", self.loop_count_range)
            },
            if self.filters.0.is_empty() { " None" } else {"\n"},
            if self.filters.0.is_empty() { "".into() } else { self.filters
                .0
                .iter()
                .map(|f| format!(" > {}", f))
                .collect::<Vec<String>>()
                .join("\n")
            }
        )
    }
}
