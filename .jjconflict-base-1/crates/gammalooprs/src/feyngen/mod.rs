pub mod diagram_generator;
pub mod feynkit;

use ahash::HashMap;
use bincode_trait_derive::Decode;
use bincode_trait_derive::Encode;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::ops::RangeInclusive;
use std::{fmt, str::FromStr};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FeynGenError {
    #[error("generation interrupted by user")]
    Interrupted,
    #[error("{0}")]
    GenericError(String),
    #[error(transparent)]
    Feynkit(#[from] feynkit::FeynkitAdapterError),
    #[error(transparent)]
    Eyre(#[from] color_eyre::Report),
}

#[derive(Debug, Clone, Eq, PartialEq, Serialize, Deserialize, JsonSchema, Encode, Decode)]

pub struct GraphGroupingOptions {
    pub numerical_sample_seed: u16,
    pub number_of_numerical_samples: usize,
    pub differentiate_particle_masses_only: bool,
    pub fully_numerical_substitution_when_comparing_numerators: bool,
    pub test_canonized_numerator: bool,
    pub symmetric_polarizations: bool,
}

impl Default for GraphGroupingOptions {
    fn default() -> Self {
        Self {
            numerical_sample_seed: 3,
            number_of_numerical_samples: 5,
            differentiate_particle_masses_only: true,
            fully_numerical_substitution_when_comparing_numerators: false,
            test_canonized_numerator: false,
            symmetric_polarizations: false,
        }
    }
}

impl fmt::Display for GraphGroupingOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "differentiate_masses_only={}, test_canonized_numerator={}, #samples={}, seed={}, fully_numerical_substitution={}, symmetric_polarizations={}",
            self.differentiate_particle_masses_only,
            self.test_canonized_numerator,
            self.number_of_numerical_samples,
            self.numerical_sample_seed,
            self.fully_numerical_substitution_when_comparing_numerators,
            self.symmetric_polarizations
        )
    }
}

#[derive(Debug, Clone, Eq, PartialEq, Serialize, Deserialize, JsonSchema, Encode, Decode)]
pub enum NumeratorAwareGraphGroupingOption {
    NoGrouping,
    OnlyDetectZeroes,
    GroupIdenticalGraphUpToSign(GraphGroupingOptions),
    GroupIdenticalGraphUpToScalarRescaling(GraphGroupingOptions),
}

impl fmt::Display for NumeratorAwareGraphGroupingOption {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::NoGrouping => "no grouping",
                Self::OnlyDetectZeroes => "only detect zero numerators",
                Self::GroupIdenticalGraphUpToSign(_opts) => "up to a sign",
                Self::GroupIdenticalGraphUpToScalarRescaling(_opts) => {
                    "up to a scalar rescaling"
                }
            }
        )
    }
}

impl NumeratorAwareGraphGroupingOption {
    pub(crate) fn get_options(&self) -> Option<&GraphGroupingOptions> {
        match self {
            Self::NoGrouping => None,
            Self::OnlyDetectZeroes => None,
            Self::GroupIdenticalGraphUpToSign(opts) => Some(opts),
            Self::GroupIdenticalGraphUpToScalarRescaling(opts) => Some(opts),
        }
    }

    #[allow(dead_code)]
    pub(crate) fn description(&self) -> String {
        format!(
            "{}{}",
            self,
            self.get_options().map_or("".into(), |o| format!("({})", o))
        )
    }
}

#[derive(Debug, Clone, Eq, PartialEq, Serialize, Deserialize, JsonSchema, Encode, Decode, Copy)]
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

#[derive(Debug, Clone, Encode, Decode, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct FeynGenFilters(pub Vec<FeynGenFilter>);

impl FeynGenFilters {
    pub(crate) fn get_coupling_orders(&self) -> Option<&HashMap<String, (usize, Option<usize>)>> {
        self.0.iter().find_map(|f| {
            if let FeynGenFilter::CouplingOrders(o) = f {
                Some(o)
            } else {
                None
            }
        })
    }

    pub(crate) fn get_perturbative_orders(&self) -> Option<&HashMap<String, usize>> {
        self.0.iter().find_map(|f| {
            if let FeynGenFilter::PerturbativeOrders(o) = f {
                Some(o)
            } else {
                None
            }
        })
    }

    pub(crate) fn get_loop_count_range(&self) -> Option<(usize, usize)> {
        self.0.iter().find_map(|f| {
            if let FeynGenFilter::LoopCountRange(o) = f {
                Some(*o)
            } else {
                None
            }
        })
    }
}

#[derive(Debug, Clone, Copy, Encode, Decode, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct SelfEnergyFilterOptions {
    pub veto_self_energy_of_massive_lines: bool,
    pub veto_self_energy_of_massless_lines: bool,
    pub veto_only_scaleless_self_energy: bool,
}

impl Default for SelfEnergyFilterOptions {
    fn default() -> Self {
        Self {
            veto_self_energy_of_massive_lines: true,
            veto_self_energy_of_massless_lines: true,
            veto_only_scaleless_self_energy: false,
        }
    }
}

impl fmt::Display for SelfEnergyFilterOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut descr = vec![];
        if self.veto_self_energy_of_massive_lines && !self.veto_self_energy_of_massless_lines {
            descr.push("only massive legs")
        } else if !self.veto_self_energy_of_massive_lines && self.veto_self_energy_of_massless_lines
        {
            descr.push("only massless legs")
        };
        if self.veto_only_scaleless_self_energy {
            descr.push("only scaleless self-energies")
        };
        write!(f, "{}", descr.join(" | "))
    }
}

#[derive(Debug, Clone, Copy, Encode, Decode, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct SnailFilterOptions {
    pub veto_snails_attached_to_massive_lines: bool,
    pub veto_snails_attached_to_massless_lines: bool,
    pub veto_only_scaleless_snails: bool,
}

impl Default for SnailFilterOptions {
    fn default() -> Self {
        Self {
            veto_snails_attached_to_massive_lines: false,
            veto_snails_attached_to_massless_lines: true,
            veto_only_scaleless_snails: false,
        }
    }
}

impl fmt::Display for SnailFilterOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut descr = vec![];
        if self.veto_snails_attached_to_massive_lines
            && !self.veto_snails_attached_to_massless_lines
        {
            descr.push("only attached to massive legs")
        } else if !self.veto_snails_attached_to_massive_lines
            && self.veto_snails_attached_to_massless_lines
        {
            descr.push("only attached to massless legs")
        };
        if self.veto_only_scaleless_snails {
            descr.push("only scaleless snails")
        };
        write!(f, "{}", descr.join(" | "))
    }
}

#[derive(Debug, Clone, Copy, Encode, Decode, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct TadpolesFilterOptions {
    pub veto_tadpoles_attached_to_massive_lines: bool,
    pub veto_tadpoles_attached_to_massless_lines: bool,
    pub veto_only_scaleless_tadpoles: bool,
}

impl Default for TadpolesFilterOptions {
    fn default() -> Self {
        Self {
            veto_tadpoles_attached_to_massive_lines: true,
            veto_tadpoles_attached_to_massless_lines: true,
            veto_only_scaleless_tadpoles: false,
        }
    }
}

impl fmt::Display for TadpolesFilterOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut descr = vec![];
        if self.veto_tadpoles_attached_to_massive_lines
            && !self.veto_tadpoles_attached_to_massless_lines
        {
            descr.push("only attached to massive legs")
        } else if !self.veto_tadpoles_attached_to_massive_lines
            && self.veto_tadpoles_attached_to_massless_lines
        {
            descr.push("only attached to massless legs")
        };
        if self.veto_only_scaleless_tadpoles {
            descr.push("only scaleless tadpoles")
        };
        write!(f, "{}", descr.join(" | "))
    }
}

#[derive(Debug, Clone, Copy, Encode, Decode, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct SewedFilterOptions {
    pub filter_tadpoles: bool,
}

#[derive(Debug, Clone, Encode, Decode, Serialize, Deserialize, JsonSchema, PartialEq)]
pub enum FeynGenFilter {
    SelfEnergyFilter(SelfEnergyFilterOptions),
    TadpolesFilter(TadpolesFilterOptions),
    ZeroSnailsFilter(SnailFilterOptions),
    SewedFilter(SewedFilterOptions),
    /// A list of vetoed pdgs
    ParticleVeto(Vec<i64>),
    VertexAllow(Vec<String>),
    VertexVeto(Vec<String>),
    MaxNumberOfBridges(usize),
    /// A map between the coupling order name and a range of orders, inclusive, with an optional upper bound
    CouplingOrders(HashMap<String, (usize, Option<usize>)>),
    /// A range of loop counts, inclusive
    LoopCountRange((usize, usize)),
    /// A range of blob counts, inclusive
    BlobRange(RangeInclusive<usize>),
    SpectatorRange(RangeInclusive<usize>),
    PerturbativeOrders(HashMap<String, usize>),
    FermionLoopCountRange((usize, usize)),
    FactorizedLoopTopologiesCountRange((usize, usize)),
}

impl fmt::Display for FeynGenFilter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::SelfEnergyFilter(opts) => format!("NoExternalSelfEnergy({})", opts),
                Self::ParticleVeto(pdgs) => format!(
                    "ParticleVeto({})",
                    pdgs.iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<String>>()
                        .join("|")
                ),
                Self::VertexVeto(vetos) => format!(
                    "VertexVeto({})",
                    vetos
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<String>>()
                        .join("|")
                ),
                Self::VertexAllow(allowed) => format!(
                    "VertexAllow({})",
                    allowed
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<String>>()
                        .join("|")
                ),
                Self::SpectatorRange(r) => format!("SpectatorRange({:?})", r),
                Self::BlobRange(r) => format!("BlobRange({:?})", r),
                Self::MaxNumberOfBridges(n) => format!("MaxNumberOfBridges({})", n),
                Self::TadpolesFilter(opts) => format!("NoTadpoles({})", opts),
                Self::ZeroSnailsFilter(opts) => format!("NoZeroSnails({})", opts),
                Self::CouplingOrders(orders) => format!(
                    "CouplingOrders({})",
                    orders
                        .iter()
                        .map(|(k, (v_min, v_max_opt))| {
                            if let Some(v_max) = v_max_opt {
                                if v_min == v_max {
                                    format!("{}=={}", k, v_min)
                                } else {
                                    format!("{}=[{}..{}]", k, v_min, v_max)
                                }
                            } else {
                                format!("{}>={}", k, v_min)
                            }
                        })
                        .collect::<Vec<String>>()
                        .join("|")
                ),
                Self::PerturbativeOrders(orders) => format!(
                    "PerturbativeOrders({})",
                    orders
                        .iter()
                        .map(|(k, v)| format!("{}={}", k, v))
                        .collect::<Vec<String>>()
                        .join("|")
                ),
                Self::LoopCountRange((loop_count_min, loop_count_max)) =>
                    format!("LoopCountRange({{{},{}}})", loop_count_min, loop_count_max),
                Self::FermionLoopCountRange((loop_count_min, loop_count_max)) => format!(
                    "FermionLoopCountRange({{{},{}}})",
                    loop_count_min, loop_count_max
                ),
                Self::FactorizedLoopTopologiesCountRange((loop_count_min, loop_count_max)) =>
                    format!(
                        "NFactorizableLoopRange({{{},{}}})",
                        loop_count_min, loop_count_max
                    ),
                Self::SewedFilter(SewedFilterOptions { filter_tadpoles }) => format!(
                    "SewedCrossSectionFilter(filter_tadpoles={{{}}})",
                    filter_tadpoles
                ),
            }
        )
    }
}

#[cfg(test)]
pub mod test;
