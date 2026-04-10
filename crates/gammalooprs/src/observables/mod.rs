use crate::model::Model;
use crate::momentum::FourMomentum;
use crate::settings::RuntimeSettings;
use crate::utils::serde_utils::{
    IsDefault, is_false, is_float, is_true, is_usize, show_defaults_helper,
};
use crate::utils::{F, FloatLike};
use bincode_trait_derive::{Decode, Encode};
use eyre::{Result, eyre};
use schemars::JsonSchema;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use smallvec::{SmallVec, smallvec};
use spenso::algebra::complex::Complex;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

pub mod clustering;
pub mod events;

pub use clustering::{ClusteringResult, Jet, JetAlgorithm, JetClustering};
pub use events::{
    AdditionalWeightKey, CutInfo, Event, EventGroup, EventGroupList, GenericAdditionalWeightInfo,
    GenericEvent, GenericEventGroup, GenericEventGroupList,
};

pub type QuantitiesSettings = BTreeMap<String, QuantitySettings>;
pub type ObservablesSettings = BTreeMap<String, ObservableSettings>;
pub type SelectorsSettings = BTreeMap<String, SelectorSettings>;

#[derive(Debug, Clone, Default)]
pub struct HistogramProcessInfo {
    pub graph_names: Vec<String>,
    pub graph_to_group_id: Vec<usize>,
    pub graph_group_master_names: Vec<String>,
    pub orientation_labels_by_group: Vec<Vec<String>>,
    pub lmb_channel_labels_by_group: Vec<Vec<String>>,
}

#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum EntrySelection {
    #[default]
    All,
    LeadingOnly,
    NthOnly,
}

#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum SelectorReduction {
    #[default]
    AnyInRange,
    AllInRange,
}

#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum ObservableValueTransform {
    #[default]
    Identity,
    Log10,
}

#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum ObservablePhase {
    #[default]
    Real,
    Imag,
}

#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
pub enum ObservableFileFormat {
    #[serde(rename = "none")]
    None,
    #[serde(rename = "hwu", alias = "HwU")]
    Hwu,
    #[default]
    #[serde(rename = "json")]
    Json,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
pub enum FilterQuantity {
    #[serde(rename = "E")]
    Energy,
    #[serde(rename = "CosTheta")]
    CosThetaP,
    #[serde(rename = "PT")]
    PT,
    #[serde(rename = "y")]
    Rapidity,
    #[serde(rename = "eta")]
    PseudoRapidity,
    #[serde(rename = "Px")]
    Px,
    #[serde(rename = "Py")]
    Py,
    #[serde(rename = "Pz")]
    Pz,
    #[serde(rename = "Mass")]
    Mass,
}

impl FilterQuantity {
    pub(crate) fn setting_name(&self) -> &'static str {
        match self {
            FilterQuantity::Energy => "E",
            FilterQuantity::CosThetaP => "CosTheta",
            FilterQuantity::PT => "PT",
            FilterQuantity::Rapidity => "y",
            FilterQuantity::PseudoRapidity => "eta",
            FilterQuantity::Px => "Px",
            FilterQuantity::Py => "Py",
            FilterQuantity::Pz => "Pz",
            FilterQuantity::Mass => "Mass",
        }
    }

    fn project_momentum<T: FloatLike>(
        &self,
        momentum: &FourMomentum<F<T>>,
        incoming_beam: Option<&FourMomentum<F<T>>>,
    ) -> Option<F<T>> {
        match self {
            FilterQuantity::Energy => Some(momentum.temporal.value.clone()),
            FilterQuantity::CosThetaP => incoming_beam.map(|beam| {
                let beam_spatial = beam.spatial.clone();
                let momentum_spatial = momentum.spatial.clone();
                beam_spatial.clone() * momentum_spatial.clone()
                    / (beam_spatial.norm() * momentum_spatial.norm())
            }),
            FilterQuantity::PT => Some(momentum.pt()),
            FilterQuantity::Rapidity => Some(momentum.rapidity()),
            FilterQuantity::PseudoRapidity => Some(momentum.spatial.pseudo_rap()),
            FilterQuantity::Px => Some(momentum.spatial.px.clone()),
            FilterQuantity::Py => Some(momentum.spatial.py.clone()),
            FilterQuantity::Pz => Some(momentum.spatial.pz.clone()),
            FilterQuantity::Mass => Some(momentum.square().abs().sqrt()),
        }
    }
}

impl fmt::Display for FilterQuantity {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.setting_name())
    }
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
pub enum QuantityOrdering {
    #[serde(rename = "PT")]
    PT,
    #[serde(rename = "Energy")]
    Energy,
    #[serde(rename = "AbsRapidity")]
    AbsRapidity,
    #[serde(rename = "Quantity")]
    Quantity,
}

impl QuantityOrdering {
    pub(crate) fn setting_name(&self) -> &'static str {
        match self {
            QuantityOrdering::PT => "PT",
            QuantityOrdering::Energy => "Energy",
            QuantityOrdering::AbsRapidity => "AbsRapidity",
            QuantityOrdering::Quantity => "Quantity",
        }
    }

    fn scalar_sort_key<T: FloatLike>(
        &self,
        momentum: &FourMomentum<F<T>>,
        quantity_value: &F<T>,
    ) -> F<T> {
        match self {
            QuantityOrdering::PT => momentum.pt(),
            QuantityOrdering::Energy => momentum.temporal.value.clone(),
            QuantityOrdering::AbsRapidity => momentum.rapidity().abs(),
            QuantityOrdering::Quantity => quantity_value.clone(),
        }
    }
}

impl fmt::Display for QuantityOrdering {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.setting_name())
    }
}

#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
pub enum QuantityOrder {
    #[serde(rename = "Ascending")]
    Ascending,
    #[default]
    #[serde(rename = "Descending")]
    Descending,
}

impl QuantityOrder {
    pub(crate) fn setting_name(&self) -> &'static str {
        match self {
            QuantityOrder::Ascending => "Ascending",
            QuantityOrder::Descending => "Descending",
        }
    }
}

impl fmt::Display for QuantityOrder {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.setting_name())
    }
}

#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum QuantityComputation {
    #[default]
    Scalar,
    Count,
    Pair,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
pub enum PairQuantity {
    #[serde(rename = "DeltaR")]
    DeltaR,
}

impl PairQuantity {
    pub(crate) fn setting_name(&self) -> &'static str {
        match self {
            PairQuantity::DeltaR => "DeltaR",
        }
    }

    fn project_momenta<T: FloatLike>(
        &self,
        left: &FourMomentum<F<T>>,
        right: &FourMomentum<F<T>>,
    ) -> F<T> {
        match self {
            PairQuantity::DeltaR => left.delta_r(right),
        }
    }
}

impl fmt::Display for PairQuantity {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.setting_name())
    }
}

#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum PairingMode {
    #[default]
    AllPairs,
}

impl PairingMode {
    pub(crate) fn setting_name(&self) -> &'static str {
        match self {
            PairingMode::AllPairs => "all_pairs",
        }
    }
}

impl fmt::Display for PairingMode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.setting_name())
    }
}

fn filter_quantity_completion_example() -> FilterQuantity {
    FilterQuantity::PT
}

fn pair_quantity_completion_example() -> PairQuantity {
    PairQuantity::DeltaR
}

fn pairing_mode_completion_example() -> PairingMode {
    PairingMode::AllPairs
}

fn quantity_ordering_completion_example() -> QuantityOrdering {
    QuantityOrdering::Quantity
}

fn quantity_order_completion_example() -> QuantityOrder {
    QuantityOrder::Descending
}

#[derive(Debug, Clone, Copy)]
enum QuantitySourceKind {
    Particle,
    Jet,
}

impl QuantitySourceKind {
    fn default_scalar_ordering(self) -> QuantityOrdering {
        match self {
            QuantitySourceKind::Particle => QuantityOrdering::Quantity,
            QuantitySourceKind::Jet => QuantityOrdering::PT,
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(default, deny_unknown_fields)]
pub struct QuantityComputationSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(
        description = "Quantity computation mode: scalar projection per object, exact object count, or pairwise quantity."
    )]
    pub computation: QuantityComputation,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(
        description = "Per-object scalar quantity used when computation = \"scalar\".",
        example = filter_quantity_completion_example()
    )]
    pub quantity: Option<FilterQuantity>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(
        description = "Pairwise quantity used when computation = \"pair\".",
        example = pair_quantity_completion_example()
    )]
    pub pair_quantity: Option<PairQuantity>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(
        description = "Pairing strategy used when computation = \"pair\".",
        example = pairing_mode_completion_example()
    )]
    pub pairing: Option<PairingMode>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(
        description = "Ordering key used before entry_selection. Particle scalar quantities default to Quantity, jet scalar quantities default to PT, and pair quantities default to Quantity.",
        example = quantity_ordering_completion_example()
    )]
    pub ordering: Option<QuantityOrdering>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(
        description = "Sorting direction applied together with ordering.",
        example = quantity_order_completion_example()
    )]
    pub order: QuantityOrder,
}

impl QuantityComputationSettings {
    pub fn scalar(quantity: FilterQuantity) -> Self {
        Self {
            computation: QuantityComputation::Scalar,
            quantity: Some(quantity),
            pair_quantity: None,
            pairing: None,
            ordering: None,
            order: QuantityOrder::Descending,
        }
    }

    fn try_normalized_for_source(&self, source: QuantitySourceKind) -> Result<Self> {
        match self.computation {
            QuantityComputation::Scalar => Ok(Self {
                computation: QuantityComputation::Scalar,
                quantity: Some(self.quantity.unwrap_or(FilterQuantity::PT)),
                pair_quantity: None,
                pairing: None,
                ordering: Some(self.ordering.unwrap_or(source.default_scalar_ordering())),
                order: self.order,
            }),
            QuantityComputation::Count => Ok(Self {
                computation: QuantityComputation::Count,
                quantity: None,
                pair_quantity: None,
                pairing: None,
                ordering: None,
                order: QuantityOrder::Descending,
            }),
            QuantityComputation::Pair => {
                let ordering = self.ordering.unwrap_or(QuantityOrdering::Quantity);
                if ordering != QuantityOrdering::Quantity {
                    return Err(eyre!(
                        "Pair quantities only support ordering=Quantity, got ordering={ordering}"
                    ));
                }

                Ok(Self {
                    computation: QuantityComputation::Pair,
                    quantity: None,
                    pair_quantity: Some(self.pair_quantity.unwrap_or(PairQuantity::DeltaR)),
                    pairing: Some(self.pairing.unwrap_or_default()),
                    ordering: Some(QuantityOrdering::Quantity),
                    order: self.order,
                })
            }
        }
    }

    fn resolve_for_source(
        &self,
        source: QuantitySourceKind,
    ) -> Result<ResolvedQuantityComputation> {
        let normalized = self.try_normalized_for_source(source)?;
        match normalized.computation {
            QuantityComputation::Scalar => Ok(ResolvedQuantityComputation::Scalar {
                quantity: normalized.quantity.unwrap(),
                ordering: normalized.ordering.unwrap(),
                order: normalized.order,
            }),
            QuantityComputation::Count => Ok(ResolvedQuantityComputation::Count),
            QuantityComputation::Pair => Ok(ResolvedQuantityComputation::Pair {
                quantity: normalized.pair_quantity.unwrap(),
                pairing: normalized.pairing.unwrap(),
                order: normalized.order,
            }),
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(default, deny_unknown_fields)]
pub struct JetClusteringSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub algorithm: JetAlgorithm,
    #[serde(skip_serializing_if = "is_float::<0>")]
    pub dR: f64,
    #[serde(skip_serializing_if = "is_float::<0>")]
    pub min_jpt: f64,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    #[schemars(
        description = "Optional list of PDG IDs allowed to be clustered. Use null to derive the default from the model, or a compact list like [-1,1,21,82] in CLI key-value assignments.",
        example = clustered_pdgs_completion_example()
    )]
    pub clustered_pdgs: Option<Vec<isize>>,
}

fn clustered_pdgs_completion_example() -> Option<Vec<isize>> {
    Some(vec![-1, 1, 21, 82])
}

#[derive(Debug, Clone, PartialEq)]
struct ResolvedJetClusteringSettings {
    algorithm: JetAlgorithm,
    d_r: f64,
    min_jpt: f64,
    clustered_pdgs: Vec<isize>,
}

impl JetClusteringSettings {
    fn resolve(&self, model: Option<&Model>) -> Result<ResolvedJetClusteringSettings> {
        Ok(ResolvedJetClusteringSettings {
            algorithm: self.algorithm,
            d_r: self.dR,
            min_jpt: self.min_jpt,
            clustered_pdgs: self.resolve_clustered_pdgs(model)?,
        })
    }

    fn resolve_clustered_pdgs(&self, model: Option<&Model>) -> Result<Vec<isize>> {
        let clustered_pdgs = match &self.clustered_pdgs {
            Some(clustered_pdgs) => clustered_pdgs.clone(),
            None => Self::default_clustered_pdgs(model.ok_or_else(|| {
                eyre!(
                    "Cannot resolve default clustered_pdgs without a model. Provide a model-aware event-processing runtime or set clustered_pdgs explicitly."
                )
            })?)?,
        };
        Ok(Self::normalize_clustered_pdgs(clustered_pdgs))
    }

    fn default_clustered_pdgs(model: &Model) -> Result<Vec<isize>> {
        Ok(model
            .particles
            .iter()
            .filter(|particle| particle.is_qcd_charged())
            .map(|particle| {
                particle
                    .has_zero_resolved_mass(model)
                    .map(|has_zero_mass| has_zero_mass.then_some(particle.pdg_code))
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect())
    }

    fn normalize_clustered_pdgs(mut clustered_pdgs: Vec<isize>) -> Vec<isize> {
        clustered_pdgs.sort_unstable();
        clustered_pdgs.dedup();
        clustered_pdgs
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
pub struct ParticleQuantitySettings {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub pdgs: Vec<isize>,
    #[serde(flatten)]
    pub computation: QuantityComputationSettings,
}

impl ParticleQuantitySettings {
    fn try_normalized(&self) -> Result<Self> {
        Ok(Self {
            pdgs: self.pdgs.clone(),
            computation: self
                .computation
                .try_normalized_for_source(QuantitySourceKind::Particle)?,
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
pub struct JetQuantitySettings {
    #[serde(flatten)]
    pub clustering: JetClusteringSettings,
    #[serde(flatten)]
    pub computation: QuantityComputationSettings,
}

impl JetQuantitySettings {
    fn try_normalized(&self) -> Result<Self> {
        Ok(Self {
            clustering: self.clustering.clone(),
            computation: self
                .computation
                .try_normalized_for_source(QuantitySourceKind::Jet)?,
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[serde(tag = "type", rename_all = "snake_case")]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
pub enum QuantitySettings {
    Particle(ParticleQuantitySettings),
    Jet(JetQuantitySettings),
    AFB {},
    Integral {},
    GraphId {},
    GraphGroupId {},
    OrientationId {},
    LmbChannelId {},
}

impl QuantitySettings {
    pub fn try_normalized(&self) -> Result<Self> {
        match self {
            QuantitySettings::Particle(settings) => {
                Ok(QuantitySettings::Particle(settings.try_normalized()?))
            }
            QuantitySettings::Jet(settings) => {
                Ok(QuantitySettings::Jet(settings.try_normalized()?))
            }
            QuantitySettings::AFB {} => Ok(QuantitySettings::AFB {}),
            QuantitySettings::Integral {} => Ok(QuantitySettings::Integral {}),
            QuantitySettings::GraphId {} => Ok(QuantitySettings::GraphId {}),
            QuantitySettings::GraphGroupId {} => Ok(QuantitySettings::GraphGroupId {}),
            QuantitySettings::OrientationId {} => Ok(QuantitySettings::OrientationId {}),
            QuantitySettings::LmbChannelId {} => Ok(QuantitySettings::LmbChannelId {}),
        }
    }

    pub fn normalized(&self) -> Self {
        self.try_normalized()
            .expect("quantity settings should normalize successfully")
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
pub struct ValueRangeSelectorSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(description = "Optional lower bound. Use null to disable the lower cut.")]
    pub min: Option<f64>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(description = "Optional upper bound. Use null to disable the upper cut.")]
    pub max: Option<f64>,
    #[serde(default, skip_serializing_if = "is_default_selector_reduction")]
    pub reduction: SelectorReduction,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
pub struct DiscreteRangeSelectorSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(description = "Optional inclusive lower bound. Use null to disable.")]
    pub min: Option<isize>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(description = "Optional inclusive upper bound. Use null to disable.")]
    pub max: Option<isize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
pub struct CountRangeSelectorSettings {
    pub min_count: usize,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub max_count: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[serde(tag = "selector", rename_all = "snake_case")]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
pub enum SelectorDefinitionSettings {
    ValueRange(ValueRangeSelectorSettings),
    DiscreteRange(DiscreteRangeSelectorSettings),
    CountRange(CountRangeSelectorSettings),
}

#[derive(Debug, Clone, Encode, Decode, PartialEq)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
pub struct SelectorSettings {
    pub quantity: String,
    pub active: bool,
    pub entry_selection: EntrySelection,
    pub entry_index: usize,
    pub selector: SelectorDefinitionSettings,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
enum ValueRangeSelectorSerdeTag {
    ValueRange,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
enum DiscreteRangeSelectorSerdeTag {
    DiscreteRange,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
enum CountRangeSelectorSerdeTag {
    CountRange,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct ValueRangeSelectorSettingsSerde {
    quantity: String,
    #[serde(default = "default_true", skip_serializing_if = "is_true")]
    active: bool,
    #[serde(default, skip_serializing_if = "is_default_entry_selection")]
    entry_selection: EntrySelection,
    #[serde(default, skip_serializing_if = "is_default_entry_index")]
    entry_index: usize,
    selector: ValueRangeSelectorSerdeTag,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    min: Option<f64>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    max: Option<f64>,
    #[serde(default, skip_serializing_if = "is_default_selector_reduction")]
    reduction: SelectorReduction,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct DiscreteRangeSelectorSettingsSerde {
    quantity: String,
    #[serde(default = "default_true", skip_serializing_if = "is_true")]
    active: bool,
    #[serde(default, skip_serializing_if = "is_default_entry_selection")]
    entry_selection: EntrySelection,
    #[serde(default, skip_serializing_if = "is_default_entry_index")]
    entry_index: usize,
    selector: DiscreteRangeSelectorSerdeTag,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    min: Option<isize>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    max: Option<isize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct CountRangeSelectorSettingsSerde {
    quantity: String,
    #[serde(default = "default_true", skip_serializing_if = "is_true")]
    active: bool,
    #[serde(default, skip_serializing_if = "is_default_entry_selection")]
    entry_selection: EntrySelection,
    #[serde(default, skip_serializing_if = "is_default_entry_index")]
    entry_index: usize,
    selector: CountRangeSelectorSerdeTag,
    min_count: usize,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    max_count: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[allow(clippy::enum_variant_names)]
#[serde(untagged)]
enum SelectorSettingsSerde {
    ValueRange(ValueRangeSelectorSettingsSerde),
    DiscreteRange(DiscreteRangeSelectorSettingsSerde),
    CountRange(CountRangeSelectorSettingsSerde),
}

impl JsonSchema for SelectorSettings {
    fn schema_name() -> std::borrow::Cow<'static, str> {
        "SelectorSettings".into()
    }

    fn json_schema(generator: &mut schemars::SchemaGenerator) -> schemars::Schema {
        <SelectorSettingsSerde as JsonSchema>::json_schema(generator)
    }
}

impl Serialize for SelectorSettings {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match &self.selector {
            SelectorDefinitionSettings::ValueRange(selector) => ValueRangeSelectorSettingsSerde {
                quantity: self.quantity.clone(),
                active: self.active,
                entry_selection: self.entry_selection,
                entry_index: self.entry_index,
                selector: ValueRangeSelectorSerdeTag::ValueRange,
                min: selector.min,
                max: selector.max,
                reduction: selector.reduction,
            }
            .serialize(serializer),
            SelectorDefinitionSettings::DiscreteRange(selector) => {
                DiscreteRangeSelectorSettingsSerde {
                    quantity: self.quantity.clone(),
                    active: self.active,
                    entry_selection: self.entry_selection,
                    entry_index: self.entry_index,
                    selector: DiscreteRangeSelectorSerdeTag::DiscreteRange,
                    min: selector.min,
                    max: selector.max,
                }
                .serialize(serializer)
            }
            SelectorDefinitionSettings::CountRange(selector) => CountRangeSelectorSettingsSerde {
                quantity: self.quantity.clone(),
                active: self.active,
                entry_selection: self.entry_selection,
                entry_index: self.entry_index,
                selector: CountRangeSelectorSerdeTag::CountRange,
                min_count: selector.min_count,
                max_count: selector.max_count,
            }
            .serialize(serializer),
        }
    }
}

impl<'de> Deserialize<'de> for SelectorSettings {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let settings = SelectorSettingsSerde::deserialize(deserializer)?;
        Ok(match settings {
            SelectorSettingsSerde::ValueRange(settings) => SelectorSettings {
                quantity: settings.quantity,
                active: settings.active,
                entry_selection: settings.entry_selection,
                entry_index: settings.entry_index,
                selector: SelectorDefinitionSettings::ValueRange(ValueRangeSelectorSettings {
                    min: settings.min,
                    max: settings.max,
                    reduction: settings.reduction,
                }),
            },
            SelectorSettingsSerde::DiscreteRange(settings) => SelectorSettings {
                quantity: settings.quantity,
                active: settings.active,
                entry_selection: settings.entry_selection,
                entry_index: settings.entry_index,
                selector: SelectorDefinitionSettings::DiscreteRange(
                    DiscreteRangeSelectorSettings {
                        min: settings.min,
                        max: settings.max,
                    },
                ),
            },
            SelectorSettingsSerde::CountRange(settings) => SelectorSettings {
                quantity: settings.quantity,
                active: settings.active,
                entry_selection: settings.entry_selection,
                entry_index: settings.entry_index,
                selector: SelectorDefinitionSettings::CountRange(CountRangeSelectorSettings {
                    min_count: settings.min_count,
                    max_count: settings.max_count,
                }),
            },
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(default, deny_unknown_fields)]
pub struct ContinuousHistogramSettings {
    #[serde(skip_serializing_if = "is_float::<0>")]
    pub x_min: f64,
    #[serde(skip_serializing_if = "is_float::<0>")]
    pub x_max: f64,
    #[serde(skip_serializing_if = "is_usize::<0>")]
    pub n_bins: usize,
    #[serde(skip_serializing_if = "is_false")]
    pub log_x_axis: bool,
    #[serde(default = "default_true", skip_serializing_if = "is_true")]
    pub log_y_axis: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(
        description = "Optional histogram title. Defaults to the observable name when omitted."
    )]
    pub title: Option<String>,
    #[serde(
        default = "default_histogram_type_description",
        skip_serializing_if = "is_default_histogram_type_description"
    )]
    #[schemars(description = "HwU TYPE description written after TYPE@ in the histogram header.")]
    pub type_description: String,
}

impl Default for ContinuousHistogramSettings {
    fn default() -> Self {
        Self {
            x_min: 0.0,
            x_max: 0.0,
            n_bins: 0,
            log_x_axis: false,
            log_y_axis: true,
            title: None,
            type_description: default_histogram_type_description(),
        }
    }
}

#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
pub enum DiscreteBinOrdering {
    #[default]
    AscendingBinId,
    ValueDescending,
    AbsValueDescending,
}

impl DiscreteBinOrdering {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AscendingBinId => "ascending_bin_id",
            Self::ValueDescending => "value_descending",
            Self::AbsValueDescending => "abs_value_descending",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum DiscreteBinDomainSettings {
    ExplicitRange { min: isize, max: isize },
    SingleBin,
    GraphIds,
    GraphGroupIds,
    OrientationIds,
    LmbChannelIds,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum DiscreteBinLabelsSettings {
    Custom { labels: Vec<String> },
    BinId,
    GraphName,
    GraphGroupMasterName,
    Orientation,
    LmbChannelEdgeIds,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(default, deny_unknown_fields)]
pub struct DiscreteHistogramSettings {
    pub domain: DiscreteBinDomainSettings,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub ordering: DiscreteBinOrdering,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub labels: Option<DiscreteBinLabelsSettings>,
    #[serde(default = "default_true", skip_serializing_if = "is_true")]
    pub log_y_axis: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    #[schemars(
        description = "Optional histogram title. Defaults to the observable name when omitted."
    )]
    pub title: Option<String>,
    #[serde(
        default = "default_histogram_type_description",
        skip_serializing_if = "is_default_histogram_type_description"
    )]
    #[schemars(description = "HwU TYPE description written after TYPE@ in the histogram header.")]
    pub type_description: String,
}

impl Default for DiscreteHistogramSettings {
    fn default() -> Self {
        Self {
            domain: DiscreteBinDomainSettings::SingleBin,
            ordering: DiscreteBinOrdering::AscendingBinId,
            labels: None,
            log_y_axis: true,
            title: None,
            type_description: default_histogram_type_description(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[serde(tag = "kind", rename_all = "snake_case")]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
pub enum HistogramSettings {
    Continuous(ContinuousHistogramSettings),
    Discrete(DiscreteHistogramSettings),
}

impl Default for HistogramSettings {
    fn default() -> Self {
        Self::Continuous(ContinuousHistogramSettings::default())
    }
}

#[derive(Debug, Clone, Encode, Decode, PartialEq)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
pub struct ObservableSettings {
    pub quantity: String,
    pub selections: Vec<String>,
    pub entry_selection: EntrySelection,
    pub entry_index: usize,
    pub value_transform: ObservableValueTransform,
    pub phase: ObservablePhase,
    pub misbinning_max_normalized_distance: Option<f64>,
    pub histogram: HistogramSettings,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
enum HistogramKindSerde {
    #[default]
    Continuous,
    Discrete,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct ObservableSettingsSerde {
    quantity: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    selections: Vec<String>,
    #[serde(default, skip_serializing_if = "is_default_entry_selection")]
    entry_selection: EntrySelection,
    #[serde(default, skip_serializing_if = "is_default_entry_index")]
    entry_index: usize,
    #[serde(default, skip_serializing_if = "is_default_value_transform")]
    value_transform: ObservableValueTransform,
    #[serde(default, skip_serializing_if = "is_default_observable_phase")]
    phase: ObservablePhase,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    misbinning_max_normalized_distance: Option<f64>,
    #[serde(default)]
    kind: Option<HistogramKindSerde>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    x_min: Option<f64>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    x_max: Option<f64>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    n_bins: Option<usize>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    log_x_axis: Option<bool>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    log_y_axis: Option<bool>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    title: Option<String>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    type_description: Option<String>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    domain: Option<DiscreteBinDomainSettings>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    ordering: Option<DiscreteBinOrdering>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    labels: Option<DiscreteBinLabelsSettings>,
}

impl JsonSchema for ObservableSettings {
    fn schema_name() -> std::borrow::Cow<'static, str> {
        "ObservableSettings".into()
    }

    fn json_schema(generator: &mut schemars::SchemaGenerator) -> schemars::Schema {
        <ObservableSettingsSerde as JsonSchema>::json_schema(generator)
    }
}

impl Serialize for ObservableSettings {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut raw = ObservableSettingsSerde {
            quantity: self.quantity.clone(),
            selections: self.selections.clone(),
            entry_selection: self.entry_selection,
            entry_index: self.entry_index,
            value_transform: self.value_transform,
            phase: self.phase,
            misbinning_max_normalized_distance: self.misbinning_max_normalized_distance,
            kind: None,
            x_min: None,
            x_max: None,
            n_bins: None,
            log_x_axis: None,
            log_y_axis: None,
            title: None,
            type_description: None,
            domain: None,
            ordering: None,
            labels: None,
        };
        match &self.histogram {
            HistogramSettings::Continuous(histogram) => {
                raw.kind = Some(HistogramKindSerde::Continuous);
                raw.x_min = Some(histogram.x_min);
                raw.x_max = Some(histogram.x_max);
                raw.n_bins = Some(histogram.n_bins);
                raw.log_x_axis = Some(histogram.log_x_axis);
                raw.log_y_axis = Some(histogram.log_y_axis);
                raw.title = histogram.title.clone();
                raw.type_description = Some(histogram.type_description.clone());
            }
            HistogramSettings::Discrete(histogram) => {
                raw.kind = Some(HistogramKindSerde::Discrete);
                raw.log_y_axis = Some(histogram.log_y_axis);
                raw.title = histogram.title.clone();
                raw.type_description = Some(histogram.type_description.clone());
                raw.domain = Some(histogram.domain.clone());
                raw.ordering = Some(histogram.ordering);
                raw.labels = histogram.labels.clone();
            }
        }
        raw.serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for ObservableSettings {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let raw = ObservableSettingsSerde::deserialize(deserializer)?;
        let inferred_kind = raw.kind.unwrap_or_else(|| {
            if raw.domain.is_some() || raw.ordering.is_some() || raw.labels.is_some() {
                HistogramKindSerde::Discrete
            } else {
                HistogramKindSerde::Continuous
            }
        });
        let histogram = match inferred_kind {
            HistogramKindSerde::Continuous => {
                if raw.domain.is_some() || raw.ordering.is_some() || raw.labels.is_some() {
                    return Err(serde::de::Error::custom(
                        "continuous observables cannot specify discrete histogram fields",
                    ));
                }
                HistogramSettings::Continuous(ContinuousHistogramSettings {
                    x_min: raw.x_min.unwrap_or_default(),
                    x_max: raw.x_max.unwrap_or_default(),
                    n_bins: raw.n_bins.unwrap_or_default(),
                    log_x_axis: raw.log_x_axis.unwrap_or(false),
                    log_y_axis: raw.log_y_axis.unwrap_or(true),
                    title: raw.title,
                    type_description: raw
                        .type_description
                        .unwrap_or_else(default_histogram_type_description),
                })
            }
            HistogramKindSerde::Discrete => {
                if raw.x_min.is_some()
                    || raw.x_max.is_some()
                    || raw.n_bins.is_some()
                    || raw.log_x_axis.is_some()
                {
                    return Err(serde::de::Error::custom(
                        "discrete observables cannot specify continuous histogram fields",
                    ));
                }
                HistogramSettings::Discrete(DiscreteHistogramSettings {
                    domain: raw.domain.unwrap_or(DiscreteBinDomainSettings::SingleBin),
                    ordering: raw.ordering.unwrap_or_default(),
                    labels: raw.labels,
                    log_y_axis: raw.log_y_axis.unwrap_or(true),
                    title: raw.title,
                    type_description: raw
                        .type_description
                        .unwrap_or_else(default_histogram_type_description),
                })
            }
        };
        Ok(ObservableSettings {
            quantity: raw.quantity,
            selections: raw.selections,
            entry_selection: raw.entry_selection,
            entry_index: raw.entry_index,
            value_transform: raw.value_transform,
            phase: raw.phase,
            misbinning_max_normalized_distance: raw.misbinning_max_normalized_distance,
            histogram,
        })
    }
}

fn default_true() -> bool {
    true
}

fn default_histogram_type_description() -> String {
    "AL".to_string()
}

fn is_default_histogram_type_description(value: &String) -> bool {
    show_defaults_helper(value == &default_histogram_type_description())
}

fn is_default_entry_selection(selection: &EntrySelection) -> bool {
    show_defaults_helper(selection == &EntrySelection::All)
}

fn is_default_entry_index(index: &usize) -> bool {
    show_defaults_helper(*index == 0)
}

fn is_default_selector_reduction(reduction: &SelectorReduction) -> bool {
    show_defaults_helper(reduction == &SelectorReduction::AnyInRange)
}

fn is_default_value_transform(transform: &ObservableValueTransform) -> bool {
    show_defaults_helper(transform == &ObservableValueTransform::Identity)
}

fn is_default_observable_phase(phase: &ObservablePhase) -> bool {
    show_defaults_helper(phase == &ObservablePhase::Real)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ObservableCoordinateKind {
    Continuous,
    Discrete,
}

#[derive(Debug, Clone)]
pub enum ObservableCoordinate<T: FloatLike> {
    Continuous(F<T>),
    Discrete(isize),
}

#[derive(Debug, Clone)]
pub struct ObservableEntry<T: FloatLike> {
    pub coordinate: ObservableCoordinate<T>,
    pub weight_modifier: Complex<F<T>>,
}

impl<T: FloatLike> ObservableEntry<T> {
    fn unit_continuous(value: F<T>) -> Self {
        Self {
            weight_modifier: unit_complex(&value),
            coordinate: ObservableCoordinate::Continuous(value),
        }
    }

    fn unit_discrete(bin_id: isize, reference: &F<T>) -> Self {
        Self {
            weight_modifier: unit_complex(reference),
            coordinate: ObservableCoordinate::Discrete(bin_id),
        }
    }
}

type ObservableEntries<T> = SmallVec<[ObservableEntry<T>; 8]>;

type ClusteringHandle = usize;

#[derive(Debug, Clone)]
struct CompiledClustering {
    settings: ResolvedJetClusteringSettings,
    clustering: JetClustering,
}

#[derive(Debug, Clone, Default)]
pub(crate) struct CompiledClusteringRegistry {
    entries: Vec<CompiledClustering>,
}

impl CompiledClusteringRegistry {
    fn register(
        &mut self,
        settings: &JetClusteringSettings,
        model: Option<&Model>,
    ) -> Result<ClusteringHandle> {
        let resolved_settings = settings.resolve(model)?;
        if let Some(index) = self
            .entries
            .iter()
            .position(|entry| entry.settings == resolved_settings)
        {
            Ok(index)
        } else {
            let index = self.entries.len();
            self.entries.push(CompiledClustering {
                clustering: JetClustering::new(
                    resolved_settings.algorithm,
                    resolved_settings.d_r,
                    resolved_settings.min_jpt,
                    resolved_settings.clustered_pdgs.clone(),
                ),
                settings: resolved_settings,
            });
            Ok(index)
        }
    }

    fn get(&self, handle: ClusteringHandle) -> &JetClustering {
        &self.entries[handle].clustering
    }

    fn len(&self) -> usize {
        self.entries.len()
    }
}

#[derive(Debug, Clone)]
#[allow(clippy::enum_variant_names)]
enum SelectorCriterion {
    ValueRange {
        selection: EntrySelection,
        entry_index: usize,
        reduction: SelectorReduction,
        min: Option<f64>,
        max: Option<f64>,
    },
    DiscreteRange {
        selection: EntrySelection,
        entry_index: usize,
        min: Option<isize>,
        max: Option<isize>,
    },
    CountRange {
        selection: EntrySelection,
        entry_index: usize,
        min: usize,
        max: Option<usize>,
    },
}

impl SelectorCriterion {
    fn passes<T: FloatLike>(&self, entries: &[ObservableEntry<T>]) -> bool {
        match self {
            SelectorCriterion::ValueRange {
                selection,
                entry_index,
                reduction,
                min,
                max,
            } => {
                let selected = apply_entry_selection(entries, *selection, *entry_index);
                if selected.is_empty() {
                    return false;
                }
                match reduction {
                    SelectorReduction::AnyInRange => selected.iter().any(|entry| {
                        matches!(
                            &entry.coordinate,
                            ObservableCoordinate::Continuous(value)
                                if value_in_range(value.into_ff64().0, *min, *max)
                        )
                    }),
                    SelectorReduction::AllInRange => selected.iter().all(|entry| {
                        matches!(
                            &entry.coordinate,
                            ObservableCoordinate::Continuous(value)
                                if value_in_range(value.into_ff64().0, *min, *max)
                        )
                    }),
                }
            }
            SelectorCriterion::DiscreteRange {
                selection,
                entry_index,
                min,
                max,
            } => {
                let selected = apply_entry_selection(entries, *selection, *entry_index);
                if selected.is_empty() {
                    return false;
                }
                selected.iter().all(|entry| {
                    matches!(
                        &entry.coordinate,
                        ObservableCoordinate::Discrete(value)
                            if discrete_value_in_range(*value, *min, *max)
                    )
                })
            }
            SelectorCriterion::CountRange {
                selection,
                entry_index,
                min,
                max,
            } => {
                let count = apply_entry_selection(entries, *selection, *entry_index).len();
                count >= *min && max.is_none_or(|max| count <= max)
            }
        }
    }
}

#[derive(Debug, Clone, Default)]
struct IntegralDefinition;

impl IntegralDefinition {
    fn process_event<T: FloatLike>(&self, event: &GenericEvent<T>) -> ObservableEntries<T> {
        let reference = event.weight.re.clone();
        smallvec![ObservableEntry::unit_discrete(0, &reference)]
    }
}

#[derive(Debug, Clone)]
enum ResolvedQuantityComputation {
    Scalar {
        quantity: FilterQuantity,
        ordering: QuantityOrdering,
        order: QuantityOrder,
    },
    Count,
    Pair {
        quantity: PairQuantity,
        pairing: PairingMode,
        order: QuantityOrder,
    },
}

impl ResolvedQuantityComputation {
    fn coordinate_kind(&self) -> ObservableCoordinateKind {
        match self {
            ResolvedQuantityComputation::Scalar { .. }
            | ResolvedQuantityComputation::Pair { .. } => ObservableCoordinateKind::Continuous,
            ResolvedQuantityComputation::Count => ObservableCoordinateKind::Discrete,
        }
    }

    fn process_momenta<T: FloatLike>(
        &self,
        momenta: &[&FourMomentum<F<T>>],
        incoming_beam: Option<&FourMomentum<F<T>>>,
        event: &GenericEvent<T>,
    ) -> ObservableEntries<T> {
        match self {
            ResolvedQuantityComputation::Scalar {
                quantity,
                ordering,
                order,
            } => {
                let mut entries = momenta
                    .iter()
                    .filter_map(|momentum| {
                        quantity
                            .project_momentum(momentum, incoming_beam)
                            .map(|value| SortableObservableEntry {
                                sort_key: ordering.scalar_sort_key(momentum, &value),
                                entry: ObservableEntry::unit_continuous(value),
                            })
                    })
                    .collect::<SmallVec<[SortableObservableEntry<T>; 8]>>();
                sort_observable_entries(&mut entries, *order);
                entries.into_iter().map(|entry| entry.entry).collect()
            }
            ResolvedQuantityComputation::Count => {
                let reference = event_representative_one(event);
                smallvec![ObservableEntry::unit_discrete(
                    momenta.len() as isize,
                    &reference,
                )]
            }
            ResolvedQuantityComputation::Pair {
                quantity,
                pairing,
                order,
            } => match pairing {
                PairingMode::AllPairs => {
                    let mut entries = SmallVec::<[SortableObservableEntry<T>; 8]>::new();
                    for left_index in 0..momenta.len() {
                        for right_index in left_index + 1..momenta.len() {
                            let value =
                                quantity.project_momenta(momenta[left_index], momenta[right_index]);
                            entries.push(SortableObservableEntry {
                                sort_key: value.clone(),
                                entry: ObservableEntry::unit_continuous(value),
                            });
                        }
                    }
                    sort_observable_entries(&mut entries, *order);
                    entries.into_iter().map(|entry| entry.entry).collect()
                }
            },
        }
    }

    fn supports_misbinning_mitigation(&self) -> bool {
        !matches!(self, ResolvedQuantityComputation::Count)
    }
}

#[derive(Debug, Clone)]
enum QuantitySourceDefinition {
    Particle { pdgs: Vec<isize> },
    Jet { clustering_handle: ClusteringHandle },
}

impl QuantitySourceDefinition {
    fn required_clustering_handle(&self) -> Option<ClusteringHandle> {
        match self {
            QuantitySourceDefinition::Particle { .. } => None,
            QuantitySourceDefinition::Jet { clustering_handle } => Some(*clustering_handle),
        }
    }

    fn collect_momenta<'a, T: FloatLike>(
        &self,
        event: &'a GenericEvent<T>,
    ) -> SmallVec<[&'a FourMomentum<F<T>>; 8]> {
        match self {
            QuantitySourceDefinition::Particle { pdgs } => event
                .cut_info
                .particle_pdgs
                .1
                .iter()
                .copied()
                .zip(event.kinematic_configuration.1.iter())
                .filter(|(pdg, _)| pdgs.contains(pdg))
                .map(|(_, momentum)| momentum)
                .collect(),
            QuantitySourceDefinition::Jet { clustering_handle } => event
                .cached_clustering(*clustering_handle)
                .expect("jet quantity requires precomputed clustering")
                .jets
                .iter()
                .map(|jet| &jet.momentum)
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
struct ObjectQuantityDefinition {
    source: QuantitySourceDefinition,
    computation: ResolvedQuantityComputation,
}

impl ObjectQuantityDefinition {
    fn particle(settings: &ParticleQuantitySettings) -> Result<Self> {
        Ok(Self {
            source: QuantitySourceDefinition::Particle {
                pdgs: settings.pdgs.clone(),
            },
            computation: settings
                .computation
                .resolve_for_source(QuantitySourceKind::Particle)?,
        })
    }

    fn jet(
        settings: &JetQuantitySettings,
        clustering_registry: &mut CompiledClusteringRegistry,
        model: Option<&Model>,
    ) -> Result<Self> {
        Ok(Self {
            source: QuantitySourceDefinition::Jet {
                clustering_handle: clustering_registry.register(&settings.clustering, model)?,
            },
            computation: settings
                .computation
                .resolve_for_source(QuantitySourceKind::Jet)?,
        })
    }

    fn process_event<T: FloatLike>(&self, event: &GenericEvent<T>) -> ObservableEntries<T> {
        let incoming_beam = event.kinematic_configuration.0.get(1);
        let momenta = self.source.collect_momenta(event);
        self.computation
            .process_momenta(&momenta, incoming_beam, event)
    }

    fn required_clustering_handle(&self) -> Option<ClusteringHandle> {
        self.source.required_clustering_handle()
    }

    fn supports_misbinning_mitigation(&self) -> bool {
        self.computation.supports_misbinning_mitigation()
    }

    fn coordinate_kind(&self) -> ObservableCoordinateKind {
        self.computation.coordinate_kind()
    }
}

#[derive(Debug, Clone, Copy)]
#[allow(clippy::enum_variant_names)]
enum MetadataQuantityKind {
    GraphId,
    GraphGroupId,
    OrientationId,
    LmbChannelId,
}

#[derive(Debug, Clone, Copy)]
struct MetadataQuantityDefinition {
    kind: MetadataQuantityKind,
}

impl MetadataQuantityDefinition {
    fn process_event<T: FloatLike>(&self, event: &GenericEvent<T>) -> ObservableEntries<T> {
        let value = match self.kind {
            MetadataQuantityKind::GraphId => Some(event.cut_info.graph_id as isize),
            MetadataQuantityKind::GraphGroupId => {
                event.cut_info.graph_group_id.map(|id| id as isize)
            }
            MetadataQuantityKind::OrientationId => {
                event.cut_info.orientation_id.map(|id| id as isize)
            }
            MetadataQuantityKind::LmbChannelId => {
                event.cut_info.lmb_channel_id.map(|id| id as isize)
            }
        };
        let Some(bin_id) = value else {
            return ObservableEntries::new();
        };
        let reference = event_representative_one(event);
        smallvec![ObservableEntry::unit_discrete(bin_id, &reference)]
    }
}

#[derive(Debug, Clone, Default)]
struct ForwardBackwardDefinition;

impl ForwardBackwardDefinition {
    fn process_event<T: FloatLike>(&self, event: &GenericEvent<T>) -> ObservableEntries<T> {
        let Some(beam) = event.kinematic_configuration.0.get(1).cloned() else {
            return ObservableEntries::new();
        };

        let Some(outgoing) = event
            .cut_info
            .particle_pdgs
            .1
            .iter()
            .copied()
            .zip(event.kinematic_configuration.1.iter().cloned())
            .find(|(pdg, _)| (1..=6).contains(pdg))
            .map(|(_, momentum)| momentum)
        else {
            return ObservableEntries::new();
        };

        let beam_spatial = beam.spatial.clone();
        let outgoing_spatial = outgoing.spatial.clone();
        smallvec![ObservableEntry::unit_continuous(
            beam_spatial.clone() * outgoing_spatial.clone()
                / (beam_spatial.norm() * outgoing_spatial.norm()),
        )]
    }
}

#[derive(Debug, Clone)]
enum ObservableDefinition {
    Integral(IntegralDefinition),
    ObjectQuantity(ObjectQuantityDefinition),
    Metadata(MetadataQuantityDefinition),
    ForwardBackward(ForwardBackwardDefinition),
}

impl ObservableDefinition {
    fn from_settings(
        settings: &QuantitySettings,
        clustering_registry: &mut CompiledClusteringRegistry,
        model: Option<&Model>,
    ) -> Result<Self> {
        Ok(match settings {
            QuantitySettings::Particle(settings) => {
                ObservableDefinition::ObjectQuantity(ObjectQuantityDefinition::particle(settings)?)
            }
            QuantitySettings::Jet(settings) => ObservableDefinition::ObjectQuantity(
                ObjectQuantityDefinition::jet(settings, clustering_registry, model)?,
            ),
            QuantitySettings::AFB {} => {
                ObservableDefinition::ForwardBackward(ForwardBackwardDefinition)
            }
            QuantitySettings::Integral {} => ObservableDefinition::Integral(IntegralDefinition),
            QuantitySettings::GraphId {} => {
                ObservableDefinition::Metadata(MetadataQuantityDefinition {
                    kind: MetadataQuantityKind::GraphId,
                })
            }
            QuantitySettings::GraphGroupId {} => {
                ObservableDefinition::Metadata(MetadataQuantityDefinition {
                    kind: MetadataQuantityKind::GraphGroupId,
                })
            }
            QuantitySettings::OrientationId {} => {
                ObservableDefinition::Metadata(MetadataQuantityDefinition {
                    kind: MetadataQuantityKind::OrientationId,
                })
            }
            QuantitySettings::LmbChannelId {} => {
                ObservableDefinition::Metadata(MetadataQuantityDefinition {
                    kind: MetadataQuantityKind::LmbChannelId,
                })
            }
        })
    }

    fn process_event<T: FloatLike>(&self, event: &GenericEvent<T>) -> ObservableEntries<T> {
        match self {
            ObservableDefinition::Integral(definition) => definition.process_event(event),
            ObservableDefinition::ObjectQuantity(definition) => definition.process_event(event),
            ObservableDefinition::Metadata(definition) => definition.process_event(event),
            ObservableDefinition::ForwardBackward(definition) => definition.process_event(event),
        }
    }

    fn required_clustering_handle(&self) -> Option<ClusteringHandle> {
        match self {
            ObservableDefinition::ObjectQuantity(definition) => {
                definition.required_clustering_handle()
            }
            _ => None,
        }
    }

    fn supports_misbinning_mitigation(&self) -> bool {
        match self {
            ObservableDefinition::ObjectQuantity(definition) => {
                definition.supports_misbinning_mitigation()
            }
            ObservableDefinition::Integral(_) | ObservableDefinition::Metadata(_) => false,
            ObservableDefinition::ForwardBackward(_) => true,
        }
    }

    fn coordinate_kind(&self) -> ObservableCoordinateKind {
        match self {
            ObservableDefinition::Integral(_) | ObservableDefinition::Metadata(_) => {
                ObservableCoordinateKind::Discrete
            }
            ObservableDefinition::ObjectQuantity(definition) => definition.coordinate_kind(),
            ObservableDefinition::ForwardBackward(_) => ObservableCoordinateKind::Continuous,
        }
    }
}

fn event_representative_one<T: FloatLike>(event: &GenericEvent<T>) -> F<T> {
    event
        .kinematic_configuration
        .0
        .first()
        .map(|momentum| momentum.temporal.value.one())
        .or_else(|| {
            event
                .kinematic_configuration
                .1
                .first()
                .map(|momentum| momentum.temporal.value.one())
        })
        .unwrap_or_else(|| event.weight.re.one())
}

#[derive(Debug, Clone)]
struct SortableObservableEntry<T: FloatLike> {
    sort_key: F<T>,
    entry: ObservableEntry<T>,
}

fn sort_observable_entries<T: FloatLike>(
    entries: &mut SmallVec<[SortableObservableEntry<T>; 8]>,
    order: QuantityOrder,
) {
    entries.sort_by(|lhs, rhs| compare_sort_keys(&lhs.sort_key, &rhs.sort_key, order));
}

fn compare_sort_keys<T: FloatLike>(lhs: &F<T>, rhs: &F<T>, order: QuantityOrder) -> Ordering {
    let base = lhs.partial_cmp(rhs).unwrap_or(Ordering::Equal);
    match order {
        QuantityOrder::Ascending => base,
        QuantityOrder::Descending => base.reverse(),
    }
}

fn ensure_event_clustering<T: FloatLike>(
    event: &mut GenericEvent<T>,
    clustering_handle: ClusteringHandle,
    clustering_registry: &CompiledClusteringRegistry,
) {
    event.ensure_clustering_slots(clustering_registry.len());
    if event.cached_clustering(clustering_handle).is_some() {
        return;
    }

    let clustering = clustering_registry.get(clustering_handle).clone();
    let clustering_result = clustering.process_event(&*event);
    event.derived_observable_data.clustered_jets[clustering_handle] = Some(clustering_result);
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct ObservableBinAccumulator {
    pub sum_weights: f64,
    new_sum_weights: f64,
    pub sum_weights_squared: f64,
    new_sum_weights_squared: f64,
    pub entry_count: usize,
    new_entry_count: usize,
    pub mitigated_fill_count: usize,
    new_mitigated_fill_count: usize,
}

impl ObservableBinAccumulator {
    fn add_sample(&mut self, sample: f64, entry_count: usize, mitigated_fill_count: usize) {
        self.new_sum_weights += sample;
        self.new_sum_weights_squared += sample * sample;
        self.new_entry_count += entry_count;
        self.new_mitigated_fill_count += mitigated_fill_count;
    }

    fn merge_samples(&mut self, other: &mut ObservableBinAccumulator) {
        self.new_sum_weights += other.new_sum_weights;
        self.new_sum_weights_squared += other.new_sum_weights_squared;
        self.new_entry_count += other.new_entry_count;
        self.new_mitigated_fill_count += other.new_mitigated_fill_count;

        other.new_sum_weights = 0.0;
        other.new_sum_weights_squared = 0.0;
        other.new_entry_count = 0;
        other.new_mitigated_fill_count = 0;
    }

    fn total_sum_weights(&self) -> f64 {
        self.sum_weights + self.new_sum_weights
    }

    fn total_sum_weights_squared(&self) -> f64 {
        self.sum_weights_squared + self.new_sum_weights_squared
    }

    fn total_entry_count(&self) -> usize {
        self.entry_count + self.new_entry_count
    }

    fn total_mitigated_fill_count(&self) -> usize {
        self.mitigated_fill_count + self.new_mitigated_fill_count
    }

    fn average(&self, sample_count: usize) -> f64 {
        if sample_count == 0 {
            0.0
        } else {
            self.total_sum_weights() / sample_count as f64
        }
    }

    fn error(&self, sample_count: usize) -> f64 {
        if sample_count <= 1 {
            return 0.0;
        }

        let n = sample_count as f64;
        let sum = self.total_sum_weights();
        let sum_sq = self.total_sum_weights_squared();
        let variance_numerator = sum_sq - (sum * sum) / n;
        if !variance_numerator.is_finite() || variance_numerator <= 0.0 {
            0.0
        } else {
            (variance_numerator / (n * (n - 1.0))).sqrt()
        }
    }

    fn update_iter(&mut self) {
        self.sum_weights += self.new_sum_weights;
        self.sum_weights_squared += self.new_sum_weights_squared;
        self.entry_count += self.new_entry_count;
        self.mitigated_fill_count += self.new_mitigated_fill_count;

        self.new_sum_weights = 0.0;
        self.new_sum_weights_squared = 0.0;
        self.new_entry_count = 0;
        self.new_mitigated_fill_count = 0;
    }

    fn rescale(&mut self, factor: f64) {
        self.sum_weights *= factor;
        self.new_sum_weights *= factor;
        self.sum_weights_squared *= factor * factor;
        self.new_sum_weights_squared *= factor * factor;
    }

    fn from_snapshot(snapshot: &HistogramBinSnapshot) -> Self {
        Self {
            sum_weights: snapshot.sum_weights,
            new_sum_weights: 0.0,
            sum_weights_squared: snapshot.sum_weights_squared,
            new_sum_weights_squared: 0.0,
            entry_count: snapshot.entry_count,
            new_entry_count: 0,
            mitigated_fill_count: snapshot.mitigated_fill_count,
            new_mitigated_fill_count: 0,
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ObservableHistogramStatistics {
    pub in_range_entry_count: usize,
    pub nan_value_count: usize,
    pub mitigated_pair_count: usize,
    new_in_range_entry_count: usize,
    new_nan_value_count: usize,
    new_mitigated_pair_count: usize,
}

impl ObservableHistogramStatistics {
    fn register_in_range_entry(&mut self) {
        self.new_in_range_entry_count += 1;
    }

    fn register_nan(&mut self) {
        self.new_nan_value_count += 1;
    }

    fn register_mitigated_pair(&mut self) {
        self.new_mitigated_pair_count += 1;
    }

    fn merge_samples(&mut self, other: &mut ObservableHistogramStatistics) {
        self.new_in_range_entry_count += other.new_in_range_entry_count;
        self.new_nan_value_count += other.new_nan_value_count;
        self.new_mitigated_pair_count += other.new_mitigated_pair_count;

        other.new_in_range_entry_count = 0;
        other.new_nan_value_count = 0;
        other.new_mitigated_pair_count = 0;
    }

    fn update_iter(&mut self) {
        self.in_range_entry_count += self.new_in_range_entry_count;
        self.nan_value_count += self.new_nan_value_count;
        self.mitigated_pair_count += self.new_mitigated_pair_count;

        self.new_in_range_entry_count = 0;
        self.new_nan_value_count = 0;
        self.new_mitigated_pair_count = 0;
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, JsonSchema)]
pub struct HistogramBinSnapshot {
    pub x_min: Option<f64>,
    pub x_max: Option<f64>,
    pub bin_id: Option<isize>,
    pub label: Option<String>,
    pub entry_count: usize,
    pub sum_weights: f64,
    pub sum_weights_squared: f64,
    pub mitigated_fill_count: usize,
}

impl HistogramBinSnapshot {
    pub fn average(&self, sample_count: usize) -> f64 {
        if sample_count == 0 {
            0.0
        } else {
            self.sum_weights / sample_count as f64
        }
    }

    pub fn error(&self, sample_count: usize) -> f64 {
        if sample_count <= 1 {
            return 0.0;
        }

        let n = sample_count as f64;
        let variance_numerator =
            self.sum_weights_squared - (self.sum_weights * self.sum_weights) / n;
        if !variance_numerator.is_finite() || variance_numerator <= 0.0 {
            0.0
        } else {
            (variance_numerator / (n * (n - 1.0))).sqrt()
        }
    }

    fn merge_in_place(&mut self, other: &HistogramBinSnapshot) {
        self.entry_count += other.entry_count;
        self.sum_weights += other.sum_weights;
        self.sum_weights_squared += other.sum_weights_squared;
        self.mitigated_fill_count += other.mitigated_fill_count;
    }

    fn rescale_in_place(&mut self, factor: f64) {
        self.sum_weights *= factor;
        self.sum_weights_squared *= factor * factor;
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, JsonSchema)]
pub struct HistogramStatisticsSnapshot {
    pub in_range_entry_count: usize,
    pub nan_value_count: usize,
    pub mitigated_pair_count: usize,
}

#[derive(
    Debug, Clone, Copy, Default, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema,
)]
#[serde(rename_all = "snake_case")]
pub enum HistogramSnapshotKind {
    #[default]
    Continuous,
    Discrete,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, JsonSchema)]
pub struct HistogramSnapshot {
    pub kind: HistogramSnapshotKind,
    pub title: String,
    pub type_description: String,
    pub phase: ObservablePhase,
    pub value_transform: ObservableValueTransform,
    pub supports_misbinning_mitigation: bool,
    pub x_min: Option<f64>,
    pub x_max: Option<f64>,
    pub sample_count: usize,
    pub log_x_axis: bool,
    pub log_y_axis: bool,
    pub discrete_min_bin_id: Option<isize>,
    pub discrete_ordering: Option<DiscreteBinOrdering>,
    pub bins: Vec<HistogramBinSnapshot>,
    pub underflow_bin: HistogramBinSnapshot,
    pub overflow_bin: HistogramBinSnapshot,
    pub statistics: HistogramStatisticsSnapshot,
}

impl HistogramSnapshot {
    pub fn to_json_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let writer = BufWriter::new(File::create(path.as_ref())?);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }

    pub fn from_json_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = BufReader::new(File::open(path.as_ref())?);
        Ok(serde_json::from_reader(reader)?)
    }

    fn write_hwu_block<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        let x_axis_mode = if self.log_x_axis { "LOG" } else { "LIN" };
        let y_axis_mode = if self.log_y_axis { "LOG" } else { "LIN" };

        writeln!(writer, "##& xmin & xmax & central value & dy &\n")?;
        writeln!(
            writer,
            "<histogram> {} \"{} |X_AXIS@{} |Y_AXIS@{} |TYPE@{}\"",
            self.bins.len(),
            self.title,
            x_axis_mode,
            y_axis_mode,
            self.type_description,
        )?;

        for bin in &self.bins {
            let (x_min, x_max) = match self.kind {
                HistogramSnapshotKind::Continuous => (
                    bin.x_min.unwrap_or(self.x_min.unwrap_or_default()),
                    bin.x_max.unwrap_or(self.x_max.unwrap_or_default()),
                ),
                HistogramSnapshotKind::Discrete => {
                    let bin_id = bin.bin_id.expect("discrete histogram bin missing bin id");
                    (bin_id as f64, (bin_id + 1) as f64)
                }
            };
            writeln!(
                writer,
                "  {:.8e}   {:.8e}   {:.8e}   {:.8e}",
                x_min,
                x_max,
                bin.average(self.sample_count),
                bin.error(self.sample_count),
            )?;
        }

        writeln!(writer, "<\\histogram>")?;
        Ok(())
    }

    pub fn merge(&self, other: &HistogramSnapshot) -> Result<Self> {
        let mut merged = self.clone();
        merged.merge_in_place(other)?;
        Ok(merged)
    }

    pub fn merged(&self, other: &HistogramSnapshot) -> Result<Self> {
        self.merge(other)
    }

    pub fn merge_in_place(&mut self, other: &HistogramSnapshot) -> Result<()> {
        if self.title != other.title
            || self.type_description != other.type_description
            || self.phase != other.phase
            || self.value_transform != other.value_transform
            || self.supports_misbinning_mitigation != other.supports_misbinning_mitigation
            || self.kind != other.kind
            || self.x_min != other.x_min
            || self.x_max != other.x_max
            || self.log_x_axis != other.log_x_axis
            || self.log_y_axis != other.log_y_axis
            || self.discrete_min_bin_id != other.discrete_min_bin_id
            || self.bins.len() != other.bins.len()
        {
            return Err(eyre!(
                "Cannot merge incompatible histogram snapshots '{}'",
                self.title
            ));
        }

        self.sample_count += other.sample_count;
        if self.kind == HistogramSnapshotKind::Continuous {
            for (bin, other_bin) in self.bins.iter_mut().zip(other.bins.iter()) {
                bin.merge_in_place(other_bin);
            }
        } else {
            let other_by_id = other
                .bins
                .iter()
                .filter_map(|bin| bin.bin_id.map(|bin_id| (bin_id, bin)))
                .collect::<BTreeMap<_, _>>();
            for bin in &mut self.bins {
                let bin_id = bin.bin_id.ok_or_else(|| {
                    eyre!(
                        "Discrete histogram snapshot '{}' contains a bin without bin_id",
                        self.title
                    )
                })?;
                let other_bin = other_by_id.get(&bin_id).ok_or_else(|| {
                    eyre!(
                        "Cannot merge discrete histogram snapshot '{}' because bin_id {} is missing",
                        self.title,
                        bin_id
                    )
                })?;
                bin.merge_in_place(other_bin);
            }
        }
        self.underflow_bin.merge_in_place(&other.underflow_bin);
        self.overflow_bin.merge_in_place(&other.overflow_bin);
        self.statistics.in_range_entry_count += other.statistics.in_range_entry_count;
        self.statistics.nan_value_count += other.statistics.nan_value_count;
        self.statistics.mitigated_pair_count += other.statistics.mitigated_pair_count;
        Ok(())
    }

    pub fn rescale(&mut self, factor: f64) {
        for bin in &mut self.bins {
            bin.rescale_in_place(factor);
        }
        self.underflow_bin.rescale_in_place(factor);
        self.overflow_bin.rescale_in_place(factor);
    }

    pub fn rescaled(&self, factor: f64) -> Self {
        let mut scaled = self.clone();
        scaled.rescale(factor);
        scaled
    }

    pub fn change_bin_ordering(&mut self, ordering: DiscreteBinOrdering) -> Result<()> {
        if self.kind != HistogramSnapshotKind::Discrete {
            return Err(eyre!(
                "Cannot change bin ordering on continuous histogram '{}'",
                self.title
            ));
        }
        self.discrete_ordering = Some(ordering);
        sort_discrete_snapshot_bins(&mut self.bins, ordering, self.sample_count);
        Ok(())
    }

    pub fn changed_bin_ordering(&self, ordering: DiscreteBinOrdering) -> Result<Self> {
        let mut reordered = self.clone();
        reordered.change_bin_ordering(ordering)?;
        Ok(reordered)
    }

    pub fn rebin(&self, contiguous_bins: usize) -> Result<Self> {
        if self.kind == HistogramSnapshotKind::Discrete {
            return Ok(self.clone());
        }
        if contiguous_bins == 0 {
            return Err(eyre!("Rebinning factor must be strictly positive."));
        }
        if self.bins.is_empty() || !self.bins.len().is_multiple_of(contiguous_bins) {
            return Err(eyre!(
                "Rebinning factor {} does not divide the {} histogram bins exactly.",
                contiguous_bins,
                self.bins.len()
            ));
        }

        let bins = self
            .bins
            .chunks(contiguous_bins)
            .map(|chunk| {
                let mut rebinned = chunk[0].clone();
                rebinned.x_min = chunk.first().and_then(|bin| bin.x_min);
                rebinned.x_max = chunk.last().and_then(|bin| bin.x_max);
                for bin in &chunk[1..] {
                    rebinned.merge_in_place(bin);
                }
                rebinned
            })
            .collect();

        Ok(Self {
            title: self.title.clone(),
            type_description: self.type_description.clone(),
            phase: self.phase,
            value_transform: self.value_transform,
            supports_misbinning_mitigation: self.supports_misbinning_mitigation,
            kind: self.kind,
            x_min: self.x_min,
            x_max: self.x_max,
            sample_count: self.sample_count,
            log_x_axis: self.log_x_axis,
            log_y_axis: self.log_y_axis,
            discrete_min_bin_id: self.discrete_min_bin_id,
            discrete_ordering: self.discrete_ordering,
            bins,
            underflow_bin: self.underflow_bin.clone(),
            overflow_bin: self.overflow_bin.clone(),
            statistics: self.statistics.clone(),
        })
    }

    pub fn rebinned(&self, contiguous_bins: usize) -> Result<Self> {
        self.rebin(contiguous_bins)
    }

    pub fn into_accumulator_state(self) -> HistogramAccumulatorState {
        HistogramAccumulatorState::from_snapshot(&self)
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct HistogramAccumulatorState {
    pub kind: HistogramSnapshotKind,
    pub title: String,
    pub type_description: String,
    pub phase: ObservablePhase,
    pub value_transform: ObservableValueTransform,
    pub supports_misbinning_mitigation: bool,
    pub x_min: f64,
    pub x_max: f64,
    pub sample_count: usize,
    new_sample_count: usize,
    pub log_x_axis: bool,
    pub log_y_axis: bool,
    pub discrete_min_bin_id: Option<isize>,
    pub discrete_ordering: Option<DiscreteBinOrdering>,
    pub bin_labels: Vec<Option<String>>,
    pub bins: Vec<ObservableBinAccumulator>,
    pub underflow_bin: ObservableBinAccumulator,
    pub overflow_bin: ObservableBinAccumulator,
    pub statistics: ObservableHistogramStatistics,
}

impl HistogramAccumulatorState {
    #[allow(clippy::too_many_arguments)]
    pub fn continuous(
        title: String,
        type_description: String,
        phase: ObservablePhase,
        value_transform: ObservableValueTransform,
        x_min: f64,
        x_max: f64,
        log_x_axis: bool,
        log_y_axis: bool,
        n_bins: usize,
    ) -> Self {
        Self::new_continuous(
            title,
            type_description,
            phase,
            value_transform,
            true,
            x_min,
            x_max,
            log_x_axis,
            log_y_axis,
            n_bins,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn discrete(
        title: String,
        type_description: String,
        phase: ObservablePhase,
        min_bin_id: isize,
        max_bin_id: isize,
        ordering: DiscreteBinOrdering,
        log_y_axis: bool,
        bin_labels: Vec<Option<String>>,
    ) -> Result<Self> {
        Self::new_discrete(
            title,
            type_description,
            phase,
            ObservableValueTransform::Identity,
            min_bin_id,
            max_bin_id,
            ordering,
            log_y_axis,
            bin_labels,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_continuous(
        title: String,
        type_description: String,
        phase: ObservablePhase,
        value_transform: ObservableValueTransform,
        supports_misbinning_mitigation: bool,
        x_min: f64,
        x_max: f64,
        log_x_axis: bool,
        log_y_axis: bool,
        n_bins: usize,
    ) -> Self {
        Self {
            kind: HistogramSnapshotKind::Continuous,
            title,
            type_description,
            phase,
            value_transform,
            supports_misbinning_mitigation,
            x_min,
            x_max,
            sample_count: 0,
            new_sample_count: 0,
            log_x_axis,
            log_y_axis,
            discrete_min_bin_id: None,
            discrete_ordering: None,
            bin_labels: Vec::new(),
            bins: vec![ObservableBinAccumulator::default(); n_bins],
            underflow_bin: ObservableBinAccumulator::default(),
            overflow_bin: ObservableBinAccumulator::default(),
            statistics: ObservableHistogramStatistics::default(),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn new_discrete(
        title: String,
        type_description: String,
        phase: ObservablePhase,
        value_transform: ObservableValueTransform,
        min_bin_id: isize,
        max_bin_id: isize,
        ordering: DiscreteBinOrdering,
        log_y_axis: bool,
        bin_labels: Vec<Option<String>>,
    ) -> Result<Self> {
        if max_bin_id < min_bin_id {
            return Err(eyre!(
                "invalid discrete histogram range: max_bin_id ({max_bin_id}) must be >= min_bin_id ({min_bin_id})"
            ));
        }
        let n_bins = max_bin_id.saturating_sub(min_bin_id) as usize + 1;
        Ok(Self {
            kind: HistogramSnapshotKind::Discrete,
            title,
            type_description,
            phase,
            value_transform,
            supports_misbinning_mitigation: false,
            x_min: min_bin_id as f64,
            x_max: (max_bin_id + 1) as f64,
            sample_count: 0,
            new_sample_count: 0,
            log_x_axis: false,
            log_y_axis,
            discrete_min_bin_id: Some(min_bin_id),
            discrete_ordering: Some(ordering),
            bin_labels,
            bins: vec![ObservableBinAccumulator::default(); n_bins],
            underflow_bin: ObservableBinAccumulator::default(),
            overflow_bin: ObservableBinAccumulator::default(),
            statistics: ObservableHistogramStatistics::default(),
        })
    }

    fn cleared_clone(&self) -> Self {
        match self.kind {
            HistogramSnapshotKind::Continuous => Self::new_continuous(
                self.title.clone(),
                self.type_description.clone(),
                self.phase,
                self.value_transform,
                self.supports_misbinning_mitigation,
                self.x_min,
                self.x_max,
                self.log_x_axis,
                self.log_y_axis,
                self.bins.len(),
            ),
            HistogramSnapshotKind::Discrete => Self::new_discrete(
                self.title.clone(),
                self.type_description.clone(),
                self.phase,
                self.value_transform,
                self.discrete_min_bin_id
                    .expect("discrete histogram missing min bin id"),
                self.discrete_min_bin_id
                    .expect("discrete histogram missing min bin id")
                    + self.bins.len() as isize
                    - 1,
                self.discrete_ordering
                    .expect("discrete histogram missing ordering"),
                self.log_y_axis,
                self.bin_labels.clone(),
            )
            .expect("cleared_clone should preserve a valid discrete histogram range"),
        }
    }

    fn discrete_bin_id(&self, index: usize) -> isize {
        self.discrete_min_bin_id
            .expect("discrete histogram missing min bin id")
            + index as isize
    }

    fn discrete_display_indices(&self) -> Vec<usize> {
        let mut indices: Vec<_> = (0..self.bins.len()).collect();
        let ordering = self
            .discrete_ordering
            .expect("discrete histogram missing ordering");
        let sample_count = self.sample_count + self.new_sample_count;
        indices.sort_by(|lhs, rhs| {
            let lhs_bin = &self.bins[*lhs];
            let rhs_bin = &self.bins[*rhs];
            let base = match ordering {
                DiscreteBinOrdering::AscendingBinId => Ordering::Equal,
                DiscreteBinOrdering::ValueDescending => rhs_bin
                    .average(sample_count)
                    .partial_cmp(&lhs_bin.average(sample_count))
                    .unwrap_or(Ordering::Equal),
                DiscreteBinOrdering::AbsValueDescending => rhs_bin
                    .average(sample_count)
                    .abs()
                    .partial_cmp(&lhs_bin.average(sample_count).abs())
                    .unwrap_or(Ordering::Equal),
            };
            if base == Ordering::Equal {
                self.discrete_bin_id(*lhs).cmp(&self.discrete_bin_id(*rhs))
            } else {
                base
            }
        });
        indices
    }

    pub fn snapshot(&self) -> HistogramSnapshot {
        let bins = match self.kind {
            HistogramSnapshotKind::Continuous => self
                .bins
                .iter()
                .enumerate()
                .map(|(index, bin)| {
                    let x_min = (self.x_max - self.x_min) * index as f64 / self.bins.len() as f64
                        + self.x_min;
                    let x_max = (self.x_max - self.x_min) * (index + 1) as f64
                        / self.bins.len() as f64
                        + self.x_min;
                    HistogramBinSnapshot {
                        x_min: Some(x_min),
                        x_max: Some(x_max),
                        bin_id: None,
                        label: None,
                        entry_count: bin.total_entry_count(),
                        sum_weights: bin.total_sum_weights(),
                        sum_weights_squared: bin.total_sum_weights_squared(),
                        mitigated_fill_count: bin.total_mitigated_fill_count(),
                    }
                })
                .collect(),
            HistogramSnapshotKind::Discrete => self
                .discrete_display_indices()
                .into_iter()
                .map(|index| {
                    let bin = &self.bins[index];
                    HistogramBinSnapshot {
                        x_min: None,
                        x_max: None,
                        bin_id: Some(self.discrete_bin_id(index)),
                        label: self.bin_labels.get(index).cloned().flatten(),
                        entry_count: bin.total_entry_count(),
                        sum_weights: bin.total_sum_weights(),
                        sum_weights_squared: bin.total_sum_weights_squared(),
                        mitigated_fill_count: bin.total_mitigated_fill_count(),
                    }
                })
                .collect(),
        };

        HistogramSnapshot {
            kind: self.kind,
            title: self.title.clone(),
            type_description: self.type_description.clone(),
            phase: self.phase,
            value_transform: self.value_transform,
            supports_misbinning_mitigation: self.supports_misbinning_mitigation,
            x_min: (self.kind == HistogramSnapshotKind::Continuous).then_some(self.x_min),
            x_max: (self.kind == HistogramSnapshotKind::Continuous).then_some(self.x_max),
            sample_count: self.sample_count + self.new_sample_count,
            log_x_axis: self.log_x_axis,
            log_y_axis: self.log_y_axis,
            discrete_min_bin_id: self.discrete_min_bin_id,
            discrete_ordering: self.discrete_ordering,
            bins,
            underflow_bin: HistogramBinSnapshot {
                x_min: None,
                x_max: (self.kind == HistogramSnapshotKind::Continuous).then_some(self.x_min),
                bin_id: None,
                label: None,
                entry_count: self.underflow_bin.total_entry_count(),
                sum_weights: self.underflow_bin.total_sum_weights(),
                sum_weights_squared: self.underflow_bin.total_sum_weights_squared(),
                mitigated_fill_count: self.underflow_bin.total_mitigated_fill_count(),
            },
            overflow_bin: HistogramBinSnapshot {
                x_min: (self.kind == HistogramSnapshotKind::Continuous).then_some(self.x_max),
                x_max: None,
                bin_id: None,
                label: None,
                entry_count: self.overflow_bin.total_entry_count(),
                sum_weights: self.overflow_bin.total_sum_weights(),
                sum_weights_squared: self.overflow_bin.total_sum_weights_squared(),
                mitigated_fill_count: self.overflow_bin.total_mitigated_fill_count(),
            },
            statistics: HistogramStatisticsSnapshot {
                in_range_entry_count: self.statistics.in_range_entry_count
                    + self.statistics.new_in_range_entry_count,
                nan_value_count: self.statistics.nan_value_count
                    + self.statistics.new_nan_value_count,
                mitigated_pair_count: self.statistics.mitigated_pair_count
                    + self.statistics.new_mitigated_pair_count,
            },
        }
    }

    fn merge_samples(&mut self, other: &mut HistogramAccumulatorState) -> Result<()> {
        self.ensure_compatible(other)?;
        self.new_sample_count += other.new_sample_count;
        other.new_sample_count = 0;
        for (bin, other_bin) in self.bins.iter_mut().zip(other.bins.iter_mut()) {
            bin.merge_samples(other_bin);
        }
        self.underflow_bin.merge_samples(&mut other.underflow_bin);
        self.overflow_bin.merge_samples(&mut other.overflow_bin);
        self.statistics.merge_samples(&mut other.statistics);
        Ok(())
    }

    pub fn merge_in_place(&mut self, other: &mut HistogramAccumulatorState) -> Result<()> {
        self.merge_samples(other)
    }

    pub fn merged_with(&self, other: &HistogramAccumulatorState) -> Result<Self> {
        let mut merged = self.clone();
        let mut other_clone = other.clone();
        merged.merge_in_place(&mut other_clone)?;
        Ok(merged)
    }

    fn update_result(&mut self) {
        self.sample_count += self.new_sample_count;
        self.new_sample_count = 0;
        for bin in &mut self.bins {
            bin.update_iter();
        }
        self.underflow_bin.update_iter();
        self.overflow_bin.update_iter();
        self.statistics.update_iter();
    }

    pub fn update_results(&mut self) {
        self.update_result();
    }

    pub fn fill_continuous_sample(&mut self, entries: &[(f64, f64)]) -> Result<()> {
        if self.kind != HistogramSnapshotKind::Continuous {
            return Err(eyre!(
                "Cannot fill continuous entries into discrete histogram '{}'",
                self.title
            ));
        }
        self.new_sample_count += 1;
        for (value, projected_weight) in entries {
            if !value.is_finite() {
                self.statistics.register_nan();
                continue;
            }
            let Some(bin_position) =
                histogram_bin_position(*value, self.x_min, self.x_max, self.bins.len())
            else {
                self.statistics.register_nan();
                continue;
            };
            match bin_position {
                HistogramBinPosition::Underflow => {
                    self.underflow_bin.add_sample(*projected_weight, 1, 0);
                }
                HistogramBinPosition::Overflow => {
                    self.overflow_bin.add_sample(*projected_weight, 1, 0);
                }
                HistogramBinPosition::InRange(bin_index) => {
                    self.statistics.register_in_range_entry();
                    self.bins[bin_index].add_sample(*projected_weight, 1, 0);
                }
            }
        }
        Ok(())
    }

    pub fn fill_discrete_sample(&mut self, entries: &[(isize, f64)]) -> Result<()> {
        if self.kind != HistogramSnapshotKind::Discrete {
            return Err(eyre!(
                "Cannot fill discrete entries into continuous histogram '{}'",
                self.title
            ));
        }
        self.new_sample_count += 1;
        let min_bin_id = self
            .discrete_min_bin_id
            .expect("discrete histogram missing min bin id");
        for (bin_id, projected_weight) in entries {
            match histogram_bin_position_discrete(*bin_id, min_bin_id, self.bins.len()) {
                HistogramBinPosition::Underflow => {
                    self.underflow_bin.add_sample(*projected_weight, 1, 0);
                }
                HistogramBinPosition::Overflow => {
                    self.overflow_bin.add_sample(*projected_weight, 1, 0);
                }
                HistogramBinPosition::InRange(bin_index) => {
                    self.statistics.register_in_range_entry();
                    self.bins[bin_index].add_sample(*projected_weight, 1, 0);
                }
            }
        }
        Ok(())
    }

    fn write_hwu_block<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        let x_axis_mode = if self.log_x_axis { "LOG" } else { "LIN" };
        let y_axis_mode = if self.log_y_axis { "LOG" } else { "LIN" };

        writeln!(writer, "##& xmin & xmax & central value & dy &\n")?;
        writeln!(
            writer,
            "<histogram> {} \"{} |X_AXIS@{} |Y_AXIS@{} |TYPE@{}\"",
            self.bins.len(),
            self.title,
            x_axis_mode,
            y_axis_mode,
            self.type_description,
        )?;

        let display_indices = if self.kind == HistogramSnapshotKind::Discrete {
            self.discrete_display_indices()
        } else {
            (0..self.bins.len()).collect()
        };

        for index in display_indices {
            let bin = &self.bins[index];
            let (c1, c2) = if self.kind == HistogramSnapshotKind::Discrete {
                let bin_id = self.discrete_bin_id(index);
                (bin_id as f64, (bin_id + 1) as f64)
            } else {
                let c1 =
                    (self.x_max - self.x_min) * index as f64 / self.bins.len() as f64 + self.x_min;
                let c2 = (self.x_max - self.x_min) * (index + 1) as f64 / self.bins.len() as f64
                    + self.x_min;
                (c1, c2)
            };
            writeln!(
                writer,
                "  {:.8e}   {:.8e}   {:.8e}   {:.8e}",
                c1,
                c2,
                bin.average(self.sample_count + self.new_sample_count),
                bin.error(self.sample_count + self.new_sample_count),
            )?;
        }

        writeln!(writer, "<\\histogram>")?;
        Ok(())
    }

    fn from_snapshot(snapshot: &HistogramSnapshot) -> Self {
        let mut bins = snapshot.bins.clone();
        if snapshot.kind == HistogramSnapshotKind::Discrete {
            bins.sort_by_key(|bin| bin.bin_id.unwrap_or_default());
        }
        let discrete_min_bin_id = if snapshot.kind == HistogramSnapshotKind::Discrete {
            Some(
                bins.first()
                    .and_then(|bin| bin.bin_id)
                    .or(snapshot.discrete_min_bin_id)
                    .unwrap_or(0),
            )
        } else {
            None
        };
        let x_min = snapshot
            .x_min
            .unwrap_or_else(|| discrete_min_bin_id.unwrap_or_default() as f64);
        let x_max = snapshot.x_max.unwrap_or_else(|| {
            let min = discrete_min_bin_id.unwrap_or_default();
            (min + bins.len() as isize) as f64
        });
        Self {
            kind: snapshot.kind,
            title: snapshot.title.clone(),
            type_description: snapshot.type_description.clone(),
            phase: snapshot.phase,
            value_transform: snapshot.value_transform,
            supports_misbinning_mitigation: snapshot.supports_misbinning_mitigation,
            x_min,
            x_max,
            sample_count: snapshot.sample_count,
            new_sample_count: 0,
            log_x_axis: snapshot.log_x_axis,
            log_y_axis: snapshot.log_y_axis,
            discrete_min_bin_id,
            discrete_ordering: snapshot.discrete_ordering,
            bin_labels: bins.iter().map(|bin| bin.label.clone()).collect(),
            bins: bins
                .iter()
                .map(ObservableBinAccumulator::from_snapshot)
                .collect(),
            underflow_bin: ObservableBinAccumulator::from_snapshot(&snapshot.underflow_bin),
            overflow_bin: ObservableBinAccumulator::from_snapshot(&snapshot.overflow_bin),
            statistics: ObservableHistogramStatistics {
                in_range_entry_count: snapshot.statistics.in_range_entry_count,
                nan_value_count: snapshot.statistics.nan_value_count,
                mitigated_pair_count: snapshot.statistics.mitigated_pair_count,
                new_in_range_entry_count: 0,
                new_nan_value_count: 0,
                new_mitigated_pair_count: 0,
            },
        }
    }

    pub fn rebin(&self, contiguous_bins: usize) -> Result<Self> {
        self.snapshot()
            .rebin(contiguous_bins)
            .map(|s| s.into_accumulator_state())
    }

    pub fn rescale(&mut self, factor: f64) {
        for bin in &mut self.bins {
            bin.rescale(factor);
        }
        self.underflow_bin.rescale(factor);
        self.overflow_bin.rescale(factor);
    }

    pub fn rescaled(&self, factor: f64) -> Self {
        let mut scaled = self.clone();
        scaled.rescale(factor);
        scaled
    }

    pub fn change_bin_ordering(&mut self, ordering: DiscreteBinOrdering) -> Result<()> {
        if self.kind != HistogramSnapshotKind::Discrete {
            return Err(eyre!(
                "Cannot change bin ordering on continuous histogram '{}'",
                self.title
            ));
        }
        self.discrete_ordering = Some(ordering);
        Ok(())
    }

    pub fn changed_bin_ordering(&self, ordering: DiscreteBinOrdering) -> Result<Self> {
        let mut reordered = self.clone();
        reordered.change_bin_ordering(ordering)?;
        Ok(reordered)
    }

    fn ensure_compatible(&self, other: &HistogramAccumulatorState) -> Result<()> {
        if self.title != other.title
            || self.type_description != other.type_description
            || self.phase != other.phase
            || self.value_transform != other.value_transform
            || self.supports_misbinning_mitigation != other.supports_misbinning_mitigation
            || self.kind != other.kind
            || self.x_min != other.x_min
            || self.x_max != other.x_max
            || self.log_x_axis != other.log_x_axis
            || self.log_y_axis != other.log_y_axis
            || self.discrete_min_bin_id != other.discrete_min_bin_id
            || self.bin_labels != other.bin_labels
            || self.bins.len() != other.bins.len()
        {
            return Err(eyre!(
                "Cannot merge incompatible histogram accumulators '{}'",
                self.title
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ObservableAccumulatorBundle {
    pub histograms: BTreeMap<String, HistogramAccumulatorState>,
}

impl ObservableAccumulatorBundle {
    pub fn merge_samples(&mut self, other: &mut ObservableAccumulatorBundle) -> Result<()> {
        for (name, other_histogram) in other.histograms.iter_mut() {
            let histogram = self
                .histograms
                .get_mut(name)
                .ok_or_else(|| eyre!("Cannot merge unknown observable accumulator '{}'", name))?;
            histogram.merge_samples(other_histogram)?;
        }
        Ok(())
    }

    pub fn update_results(&mut self) {
        for histogram in self.histograms.values_mut() {
            histogram.update_result();
        }
    }

    pub fn snapshot_bundle(&self) -> ObservableSnapshotBundle {
        ObservableSnapshotBundle {
            histograms: self
                .histograms
                .iter()
                .map(|(name, histogram)| (name.clone(), histogram.snapshot()))
                .collect(),
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, JsonSchema)]
pub struct ObservableSnapshotBundle {
    pub histograms: BTreeMap<String, HistogramSnapshot>,
}

impl ObservableSnapshotBundle {
    pub fn to_json_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let writer = BufWriter::new(File::create(path.as_ref())?);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }

    pub fn from_json_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = BufReader::new(File::open(path.as_ref())?);
        Ok(serde_json::from_reader(reader)?)
    }

    pub fn write_hwu_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let mut writer = BufWriter::new(File::create(path.as_ref())?);
        for histogram in self.histograms.values() {
            histogram.write_hwu_block(&mut writer)?;
        }
        Ok(())
    }

    pub fn into_accumulator_bundle(self) -> ObservableAccumulatorBundle {
        ObservableAccumulatorBundle {
            histograms: self
                .histograms
                .into_iter()
                .map(|(name, histogram)| (name, histogram.into_accumulator_state()))
                .collect(),
        }
    }

    pub fn merge_in_place(&mut self, other: &ObservableSnapshotBundle) -> Result<()> {
        for (name, other_histogram) in other.histograms.iter() {
            let histogram = self
                .histograms
                .get_mut(name)
                .ok_or_else(|| eyre!("Cannot merge unknown histogram snapshot '{}'", name))?;
            histogram.merge_in_place(other_histogram)?;
        }
        Ok(())
    }
}

trait EventSelector {
    fn process_event<T: FloatLike>(
        &mut self,
        event: &mut GenericEvent<T>,
        clustering_registry: &CompiledClusteringRegistry,
    ) -> bool;
}

#[derive(Debug, Clone)]
struct ConfiguredSelector {
    definition: ObservableDefinition,
    criteria: Vec<SelectorCriterion>,
}

impl ConfiguredSelector {
    pub(crate) fn new(
        settings: &SelectorSettings,
        quantities: &QuantitiesSettings,
        clustering_registry: &mut CompiledClusteringRegistry,
        model: Option<&Model>,
    ) -> Result<Self> {
        let quantity = quantities.get(&settings.quantity).ok_or_else(|| {
            eyre!(
                "Selector references unknown quantity '{}'",
                settings.quantity
            )
        })?;

        let definition = ObservableDefinition::from_settings(quantity, clustering_registry, model)?;
        let coordinate_kind = definition.coordinate_kind();

        let criterion = match &settings.selector {
            SelectorDefinitionSettings::ValueRange(selector) => {
                if coordinate_kind != ObservableCoordinateKind::Continuous {
                    return Err(eyre!(
                        "Selector '{}' uses value_range with a discrete quantity '{}'",
                        settings.quantity,
                        settings.quantity
                    ));
                }
                SelectorCriterion::ValueRange {
                    selection: settings.entry_selection,
                    entry_index: settings.entry_index,
                    reduction: selector.reduction,
                    min: selector.min,
                    max: selector.max,
                }
            }
            SelectorDefinitionSettings::DiscreteRange(selector) => {
                if coordinate_kind != ObservableCoordinateKind::Discrete {
                    return Err(eyre!(
                        "Selector '{}' uses discrete_range with a continuous quantity '{}'",
                        settings.quantity,
                        settings.quantity
                    ));
                }
                SelectorCriterion::DiscreteRange {
                    selection: settings.entry_selection,
                    entry_index: settings.entry_index,
                    min: selector.min,
                    max: selector.max,
                }
            }
            SelectorDefinitionSettings::CountRange(selector) => SelectorCriterion::CountRange {
                selection: settings.entry_selection,
                entry_index: settings.entry_index,
                min: selector.min_count,
                max: selector.max_count,
            },
        };

        Ok(Self {
            definition,
            criteria: vec![criterion],
        })
    }

    fn required_clustering_handle(&self) -> Option<ClusteringHandle> {
        self.definition.required_clustering_handle()
    }

    fn passes_event<T: FloatLike>(
        &self,
        event: &mut GenericEvent<T>,
        clustering_registry: &CompiledClusteringRegistry,
    ) -> bool {
        if let Some(clustering_handle) = self.definition.required_clustering_handle() {
            ensure_event_clustering(event, clustering_handle, clustering_registry);
        }
        let entries = self.definition.process_event(event);
        self.criteria
            .iter()
            .all(|criterion| criterion.passes(&entries))
    }

    fn passes_prepared_event<T: FloatLike>(&self, event: &GenericEvent<T>) -> bool {
        let entries = self.definition.process_event(event);
        self.criteria
            .iter()
            .all(|criterion| criterion.passes(&entries))
    }
}

impl EventSelector for ConfiguredSelector {
    fn process_event<T: FloatLike>(
        &mut self,
        event: &mut GenericEvent<T>,
        clustering_registry: &CompiledClusteringRegistry,
    ) -> bool {
        self.passes_event(event, clustering_registry)
    }
}

#[derive(Clone)]
enum Selectors {
    Configured(ConfiguredSelector),
}

impl Selectors {
    fn process_event<T: FloatLike>(
        &mut self,
        event: &mut GenericEvent<T>,
        clustering_registry: &CompiledClusteringRegistry,
    ) -> bool {
        match self {
            Selectors::Configured(selector) => selector.process_event(event, clustering_registry),
        }
    }
}

#[derive(Clone, Default)]
pub struct EventProcessingRuntime {
    clustering_registry: CompiledClusteringRegistry,
    observable_clustering_handles: Vec<ClusteringHandle>,
    selectors: Vec<Selectors>,
    observables: BTreeMap<String, Observables>,
}

impl EventProcessingRuntime {
    pub fn from_settings(settings: &RuntimeSettings) -> Result<Self> {
        Self::from_settings_inner(settings, None, None)
    }

    pub fn from_settings_with_model(settings: &RuntimeSettings, model: &Model) -> Result<Self> {
        Self::from_settings_inner(settings, Some(model), None)
    }

    pub fn from_settings_with_model_and_process_info(
        settings: &RuntimeSettings,
        model: &Model,
        process_info: &HistogramProcessInfo,
    ) -> Result<Self> {
        Self::from_settings_inner(settings, Some(model), Some(process_info))
    }

    fn from_settings_inner(
        settings: &RuntimeSettings,
        model: Option<&Model>,
        process_info: Option<&HistogramProcessInfo>,
    ) -> Result<Self> {
        let mut clustering_registry = CompiledClusteringRegistry::default();
        let selector_registry = settings
            .selectors
            .iter()
            .map(|(name, selector)| {
                ConfiguredSelector::new(
                    selector,
                    &settings.quantities,
                    &mut clustering_registry,
                    model,
                )
                .map(|configured| (name.clone(), configured))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;
        let selectors = settings
            .selectors
            .iter()
            .filter(|(_, selector)| selector.active)
            .map(|(name, _selector)| Ok(Selectors::Configured(selector_registry[name].clone())))
            .collect::<Result<Vec<_>>>()?;

        let mut observable_clustering_handles = Vec::new();
        let observables = settings
            .observables
            .iter()
            .map(|(name, observable)| {
                Observables::from_settings(
                    name,
                    observable,
                    &settings.quantities,
                    &settings.selectors,
                    &selector_registry,
                    process_info,
                    &mut clustering_registry,
                    model,
                )
                .map(|observable| {
                    for handle in observable.required_clustering_handles() {
                        if !observable_clustering_handles.contains(&handle) {
                            observable_clustering_handles.push(handle);
                        }
                    }
                    (name.clone(), observable)
                })
            })
            .collect::<Result<BTreeMap<_, _>>>()?;

        Ok(Self {
            clustering_registry,
            observable_clustering_handles,
            selectors,
            observables,
        })
    }

    pub fn has_selectors(&self) -> bool {
        !self.selectors.is_empty()
    }

    pub fn has_observables(&self) -> bool {
        !self.observables.is_empty()
    }

    fn process_event_internal<T: FloatLike>(
        &mut self,
        event: &mut GenericEvent<T>,
        prepare_observables: bool,
    ) -> bool {
        for selector in self.selectors.iter_mut() {
            if !selector.process_event(event, &self.clustering_registry) {
                return false;
            }
        }

        if prepare_observables && self.has_observables() {
            for &handle in &self.observable_clustering_handles {
                ensure_event_clustering(event, handle, &self.clustering_registry);
            }
        }
        true
    }

    pub fn process_event<T: FloatLike>(&mut self, event: &mut GenericEvent<T>) -> bool {
        self.process_event_internal(event, true)
    }

    pub(crate) fn process_event_for_selectors<T: FloatLike>(
        &mut self,
        event: &mut GenericEvent<T>,
    ) -> bool {
        self.process_event_internal(event, false)
    }

    pub fn process_event_groups<T: FloatLike>(&mut self, event_groups: &GenericEventGroupList<T>) {
        for observable in self.observables.values_mut() {
            observable.process_event_groups(event_groups, &self.clustering_registry);
        }
    }

    pub fn merge_samples(&mut self, other: &mut EventProcessingRuntime) -> Result<()> {
        for (name, other_observable) in other.observables.iter_mut() {
            let observable = self
                .observables
                .get_mut(name)
                .ok_or_else(|| eyre!("Cannot merge unknown observable '{}'", name))?;
            observable.merge_samples(other_observable)?;
        }
        Ok(())
    }

    pub fn update_results(&mut self, iter: usize) {
        for observable in self.observables.values_mut() {
            observable.update_result(iter);
        }
    }

    pub fn cleared_observable_clone(&self) -> Self {
        Self {
            clustering_registry: self.clustering_registry.clone(),
            observable_clustering_handles: self.observable_clustering_handles.clone(),
            selectors: self.selectors.clone(),
            observables: self
                .observables
                .iter()
                .map(|(name, observable)| (name.clone(), observable.cleared_clone()))
                .collect(),
        }
    }

    pub fn snapshot_bundle(&self) -> ObservableSnapshotBundle {
        ObservableSnapshotBundle {
            histograms: self
                .observables
                .iter()
                .map(|(name, observable)| (name.clone(), observable.snapshot()))
                .collect(),
        }
    }

    pub fn restore_snapshot_bundle(&mut self, bundle: &ObservableSnapshotBundle) -> Result<()> {
        if self.observables.len() != bundle.histograms.len() {
            return Err(eyre!(
                "Cannot restore observables from a snapshot with a different histogram count"
            ));
        }

        for (name, observable) in self.observables.iter_mut() {
            let snapshot = bundle
                .histograms
                .get(name)
                .ok_or_else(|| eyre!("Cannot restore unknown histogram snapshot '{}'", name))?;
            observable.restore_snapshot(snapshot)?;
        }

        Ok(())
    }

    pub fn accumulator_bundle(&self) -> ObservableAccumulatorBundle {
        ObservableAccumulatorBundle {
            histograms: self
                .observables
                .iter()
                .map(|(name, observable)| (name.clone(), observable.accumulator_state()))
                .collect(),
        }
    }

    pub fn merge_accumulator_bundle(
        &mut self,
        other: &mut ObservableAccumulatorBundle,
    ) -> Result<()> {
        for (name, other_histogram) in other.histograms.iter_mut() {
            let observable = self
                .observables
                .get_mut(name)
                .ok_or_else(|| eyre!("Cannot merge unknown observable '{}'", name))?;
            observable.merge_accumulator_state(other_histogram)?;
        }
        Ok(())
    }
}

pub(crate) trait Observable {
    fn process_event_groups<T: FloatLike>(
        &mut self,
        event_groups: &GenericEventGroupList<T>,
        clustering_registry: &CompiledClusteringRegistry,
    );

    fn merge_samples(&mut self, other: &mut Self) -> Result<()>
    where
        Self: Sized;

    fn update_result(&mut self, iter: usize);
}

#[derive(Debug, Clone)]
struct PendingBinContribution {
    value: f64,
    bin_index: usize,
    projected_weight: f64,
    mitigated: bool,
}

#[derive(Debug, Clone, Default)]
struct GroupBinContribution {
    projected_weight: f64,
    entry_count: usize,
    mitigated_fill_count: usize,
}

#[derive(Debug, Clone)]
pub struct HistogramObservable {
    definition: ObservableDefinition,
    selections: Vec<ConfiguredSelector>,
    entry_selection: EntrySelection,
    entry_index: usize,
    misbinning_max_normalized_distance: Option<f64>,
    state: HistogramAccumulatorState,
    pending_contributions: Vec<PendingBinContribution>,
    candidate_pairs: Vec<(f64, usize, usize)>,
    used_contributions: Vec<bool>,
    grouped_contributions: Vec<GroupBinContribution>,
    grouped_underflow: GroupBinContribution,
    grouped_overflow: GroupBinContribution,
    touched_bins: Vec<usize>,
}

impl HistogramObservable {
    fn new(
        definition: ObservableDefinition,
        selections: Vec<ConfiguredSelector>,
        quantities: &QuantitiesSettings,
        quantity_settings: &QuantitySettings,
        observable_name: &str,
        settings: &ObservableSettings,
        selector_settings: &SelectorsSettings,
        process_info: Option<&HistogramProcessInfo>,
    ) -> Result<Self> {
        let state = match (&settings.histogram, definition.coordinate_kind()) {
            (HistogramSettings::Continuous(histogram), ObservableCoordinateKind::Continuous) => {
                let supports_misbinning_mitigation = definition.supports_misbinning_mitigation();
                if settings.misbinning_max_normalized_distance.is_some()
                    && !supports_misbinning_mitigation
                {
                    return Err(eyre!(
                        "Observable '{}' does not support misbinning mitigation.",
                        observable_name
                    ));
                }
                HistogramAccumulatorState::new_continuous(
                    histogram
                        .title
                        .clone()
                        .unwrap_or_else(|| observable_name.to_string()),
                    histogram.type_description.clone(),
                    settings.phase,
                    settings.value_transform,
                    supports_misbinning_mitigation,
                    histogram.x_min,
                    histogram.x_max,
                    histogram.log_x_axis,
                    histogram.log_y_axis,
                    histogram.n_bins,
                )
            }
            (HistogramSettings::Discrete(histogram), ObservableCoordinateKind::Discrete) => {
                if settings.misbinning_max_normalized_distance.is_some() {
                    return Err(eyre!(
                        "Observable '{}' uses a discrete histogram, so misbinning mitigation is not supported.",
                        observable_name
                    ));
                }
                if settings.value_transform != ObservableValueTransform::Identity {
                    return Err(eyre!(
                        "Observable '{}' uses a discrete histogram, so value_transform must stay identity.",
                        observable_name
                    ));
                }
                let (min_bin_id, max_bin_id, bin_labels) = resolve_discrete_histogram_layout(
                    quantities,
                    quantity_settings,
                    observable_name,
                    histogram,
                    &settings.selections,
                    selector_settings,
                    process_info,
                )?;
                HistogramAccumulatorState::new_discrete(
                    histogram
                        .title
                        .clone()
                        .unwrap_or_else(|| observable_name.to_string()),
                    histogram.type_description.clone(),
                    settings.phase,
                    settings.value_transform,
                    min_bin_id,
                    max_bin_id,
                    histogram.ordering,
                    histogram.log_y_axis,
                    bin_labels,
                )?
            }
            (HistogramSettings::Continuous(_), ObservableCoordinateKind::Discrete) => {
                return Err(eyre!(
                    "Observable '{}' uses a discrete quantity and therefore requires a discrete histogram definition.",
                    observable_name
                ));
            }
            (HistogramSettings::Discrete(_), ObservableCoordinateKind::Continuous) => {
                return Err(eyre!(
                    "Observable '{}' uses a continuous quantity and therefore requires a continuous histogram definition.",
                    observable_name
                ));
            }
        };
        let n_bins = state.bins.len();

        Ok(Self {
            definition,
            selections,
            entry_selection: settings.entry_selection,
            entry_index: settings.entry_index,
            misbinning_max_normalized_distance: settings.misbinning_max_normalized_distance,
            state,
            pending_contributions: Vec::with_capacity(16),
            candidate_pairs: Vec::with_capacity(16),
            used_contributions: Vec::with_capacity(16),
            grouped_contributions: vec![GroupBinContribution::default(); n_bins],
            grouped_underflow: GroupBinContribution::default(),
            grouped_overflow: GroupBinContribution::default(),
            touched_bins: Vec::with_capacity(n_bins.min(16)),
        })
    }

    fn build_group_contributions<T: FloatLike>(
        &mut self,
        event_group: &GenericEventGroup<T>,
        _clustering_registry: &CompiledClusteringRegistry,
    ) {
        self.pending_contributions.clear();
        self.pending_contributions
            .reserve(event_group.len().saturating_mul(4));

        for event in event_group.iter() {
            if self
                .selections
                .iter()
                .any(|selector| !selector.passes_prepared_event(event))
            {
                continue;
            }
            let entries = self.definition.process_event(event);
            for entry in apply_entry_selection(&entries, self.entry_selection, self.entry_index) {
                let (value, bin_position) = match &entry.coordinate {
                    ObservableCoordinate::Continuous(entry_value) => {
                        let transformed_value =
                            transform_value(entry_value, self.state.value_transform);
                        let value = transformed_value.into_ff64().0;
                        if !value.is_finite() {
                            self.state.statistics.register_nan();
                            continue;
                        }

                        let Some(bin_position) = histogram_bin_position(
                            value,
                            self.state.x_min,
                            self.state.x_max,
                            self.state.bins.len(),
                        ) else {
                            self.state.statistics.register_nan();
                            continue;
                        };
                        (value, bin_position)
                    }
                    ObservableCoordinate::Discrete(bin_id) => (
                        *bin_id as f64,
                        histogram_bin_position_discrete(
                            *bin_id,
                            self.state
                                .discrete_min_bin_id
                                .expect("discrete histogram missing min bin id"),
                            self.state.bins.len(),
                        ),
                    ),
                };

                let projected_weight = self
                    .state
                    .phase
                    .project(&combined_entry_weight(event, &entry))
                    .into_ff64()
                    .0;
                if !projected_weight.is_finite() {
                    self.state.statistics.register_nan();
                    continue;
                }

                match bin_position {
                    HistogramBinPosition::Underflow => {
                        self.grouped_underflow.projected_weight += projected_weight;
                        self.grouped_underflow.entry_count += 1;
                    }
                    HistogramBinPosition::Overflow => {
                        self.grouped_overflow.projected_weight += projected_weight;
                        self.grouped_overflow.entry_count += 1;
                    }
                    HistogramBinPosition::InRange(bin_index) => {
                        self.state.statistics.register_in_range_entry();
                        self.pending_contributions.push(PendingBinContribution {
                            value,
                            bin_index,
                            projected_weight,
                            mitigated: false,
                        });
                    }
                }
            }
        }
    }

    fn mitigate_group_misbinning(&mut self) {
        if self.state.kind == HistogramSnapshotKind::Discrete {
            return;
        }
        let Some(max_distance) = self.misbinning_max_normalized_distance else {
            return;
        };

        if self.pending_contributions.len() < 2 || self.state.bins.len() < 2 || max_distance <= 0.0
        {
            return;
        }

        let bin_width = (self.state.x_max - self.state.x_min) / self.state.bins.len() as f64;
        if !bin_width.is_finite() || bin_width <= 0.0 {
            return;
        }

        self.candidate_pairs.clear();
        for i in 0..self.pending_contributions.len() {
            for j in (i + 1)..self.pending_contributions.len() {
                if self.pending_contributions[i].projected_weight == 0.0
                    || self.pending_contributions[j].projected_weight == 0.0
                {
                    continue;
                }
                if self.pending_contributions[i].projected_weight.signum()
                    == self.pending_contributions[j].projected_weight.signum()
                {
                    continue;
                }
                if self.pending_contributions[i]
                    .bin_index
                    .abs_diff(self.pending_contributions[j].bin_index)
                    != 1
                {
                    continue;
                }

                let normalized_distance = (self.pending_contributions[i].value
                    - self.pending_contributions[j].value)
                    .abs()
                    / bin_width;
                if normalized_distance <= max_distance {
                    self.candidate_pairs.push((normalized_distance, i, j));
                }
            }
        }

        self.candidate_pairs
            .sort_by(|lhs, rhs| match lhs.0.partial_cmp(&rhs.0) {
                Some(ordering) => ordering,
                None => Ordering::Greater,
            });

        self.used_contributions.clear();
        self.used_contributions
            .resize(self.pending_contributions.len(), false);
        for (_, i, j) in self.candidate_pairs.iter().copied() {
            if self.used_contributions[i] || self.used_contributions[j] {
                continue;
            }

            let averaged_weight = 0.5
                * (self.pending_contributions[i].projected_weight
                    + self.pending_contributions[j].projected_weight);
            self.pending_contributions[i].projected_weight = averaged_weight;
            self.pending_contributions[j].projected_weight = averaged_weight;
            self.pending_contributions[i].mitigated = true;
            self.pending_contributions[j].mitigated = true;
            self.used_contributions[i] = true;
            self.used_contributions[j] = true;
            self.state.statistics.register_mitigated_pair();
        }
    }

    pub fn snapshot(&self) -> HistogramSnapshot {
        self.state.snapshot()
    }

    pub fn write_to_file<P: AsRef<Path>>(
        &self,
        path: P,
        format: ObservableFileFormat,
    ) -> Result<()> {
        match format {
            ObservableFileFormat::None => Ok(()),
            ObservableFileFormat::Hwu => {
                self.write_hwu_file(path)?;
                Ok(())
            }
            ObservableFileFormat::Json => self.snapshot().to_json_file(path),
        }
    }

    fn write_hwu_file<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        let mut writer = BufWriter::new(File::create(path.as_ref())?);
        self.state.write_hwu_block(&mut writer)
    }

    fn flush_sample_contributions(&mut self) {
        if self.touched_bins.is_empty()
            && self.grouped_underflow.projected_weight == 0.0
            && self.grouped_underflow.entry_count == 0
            && self.grouped_underflow.mitigated_fill_count == 0
            && self.grouped_overflow.projected_weight == 0.0
            && self.grouped_overflow.entry_count == 0
            && self.grouped_overflow.mitigated_fill_count == 0
        {
            return;
        }

        for &bin_index in &self.touched_bins {
            let grouped = &self.grouped_contributions[bin_index];
            let bin_accumulator = &mut self.state.bins[bin_index];
            bin_accumulator.add_sample(
                grouped.projected_weight,
                grouped.entry_count,
                grouped.mitigated_fill_count,
            );
        }
        if self.grouped_underflow.projected_weight != 0.0
            || self.grouped_underflow.entry_count != 0
            || self.grouped_underflow.mitigated_fill_count != 0
        {
            self.state.underflow_bin.add_sample(
                self.grouped_underflow.projected_weight,
                self.grouped_underflow.entry_count,
                self.grouped_underflow.mitigated_fill_count,
            );
        }
        if self.grouped_overflow.projected_weight != 0.0
            || self.grouped_overflow.entry_count != 0
            || self.grouped_overflow.mitigated_fill_count != 0
        {
            self.state.overflow_bin.add_sample(
                self.grouped_overflow.projected_weight,
                self.grouped_overflow.entry_count,
                self.grouped_overflow.mitigated_fill_count,
            );
        }
    }

    fn clear_sample_contributions(&mut self) {
        self.grouped_underflow = GroupBinContribution::default();
        self.grouped_overflow = GroupBinContribution::default();
        for bin_index in self.touched_bins.drain(..) {
            self.grouped_contributions[bin_index] = GroupBinContribution::default();
        }
    }

    fn required_clustering_handles(&self) -> Vec<ClusteringHandle> {
        let mut handles = Vec::new();
        if let Some(handle) = self.definition.required_clustering_handle() {
            handles.push(handle);
        }
        for selector in &self.selections {
            if let Some(handle) = selector.required_clustering_handle()
                && !handles.contains(&handle)
            {
                handles.push(handle);
            }
        }
        handles
    }
}

fn resolve_discrete_histogram_layout(
    quantities: &QuantitiesSettings,
    quantity_settings: &QuantitySettings,
    observable_name: &str,
    histogram: &DiscreteHistogramSettings,
    selection_names: &[String],
    selector_settings: &SelectorsSettings,
    process_info: Option<&HistogramProcessInfo>,
) -> Result<(isize, isize, Vec<Option<String>>)> {
    let graph_group_context = resolve_graph_group_context(
        observable_name,
        selection_names,
        selector_settings,
        quantities,
        process_info,
    )?;
    let (min_bin_id, max_bin_id) = match histogram.domain {
        DiscreteBinDomainSettings::ExplicitRange { min, max } => (min, max),
        DiscreteBinDomainSettings::SingleBin => (0, 0),
        DiscreteBinDomainSettings::GraphIds => {
            let info = process_info.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses GraphIds but no process information is available.",
                    observable_name
                )
            })?;
            if info.graph_names.is_empty() {
                return Err(eyre!(
                    "Observable '{}' uses GraphIds but the process exposes no graphs.",
                    observable_name
                ));
            }
            (0, info.graph_names.len() as isize - 1)
        }
        DiscreteBinDomainSettings::GraphGroupIds => {
            let info = process_info.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses GraphGroupIds but no process information is available.",
                    observable_name
                )
            })?;
            if info.graph_group_master_names.is_empty() {
                return Err(eyre!(
                    "Observable '{}' uses GraphGroupIds but the process exposes no graph groups.",
                    observable_name
                ));
            }
            (0, info.graph_group_master_names.len() as isize - 1)
        }
        DiscreteBinDomainSettings::OrientationIds => {
            let info = process_info.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses OrientationIds but no process information is available.",
                    observable_name
                )
            })?;
            let group_id = graph_group_context.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses OrientationIds but no singleton graph-group context could be resolved from its selections.",
                    observable_name
                )
            })?;
            let labels = info.orientation_labels_by_group.get(group_id).ok_or_else(|| {
                eyre!(
                    "Observable '{}' resolved graph-group {} for OrientationIds, but no orientation metadata is available.",
                    observable_name,
                    group_id
                )
            })?;
            (0, labels.len() as isize - 1)
        }
        DiscreteBinDomainSettings::LmbChannelIds => {
            let info = process_info.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses LmbChannelIds but no process information is available.",
                    observable_name
                )
            })?;
            let group_id = graph_group_context.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses LmbChannelIds but no singleton graph-group context could be resolved from its selections.",
                    observable_name
                )
            })?;
            let labels = info.lmb_channel_labels_by_group.get(group_id).ok_or_else(|| {
                eyre!(
                    "Observable '{}' resolved graph-group {} for LmbChannelIds, but no LMB-channel metadata is available.",
                    observable_name,
                    group_id
                )
            })?;
            (0, labels.len() as isize - 1)
        }
    };
    if max_bin_id < min_bin_id {
        return Err(eyre!(
            "Observable '{}' has an invalid discrete histogram range [{}, {}].",
            observable_name,
            min_bin_id,
            max_bin_id
        ));
    }
    let n_bins = max_bin_id.saturating_sub(min_bin_id) as usize + 1;
    let bin_labels = match histogram.labels.as_ref() {
        None => {
            if matches!(quantity_settings, QuantitySettings::Integral {}) && n_bins == 1 {
                vec![Some("total integral".to_string())]
            } else {
                vec![None; n_bins]
            }
        }
        Some(DiscreteBinLabelsSettings::Custom { labels }) => {
            if labels.len() != n_bins {
                return Err(eyre!(
                    "Observable '{}' defines {} discrete labels for {} bins.",
                    observable_name,
                    labels.len(),
                    n_bins
                ));
            }
            labels.iter().cloned().map(Some).collect()
        }
        Some(DiscreteBinLabelsSettings::BinId) => (0..n_bins)
            .map(|offset| Some(format!("#{}", min_bin_id + offset as isize)))
            .collect(),
        Some(DiscreteBinLabelsSettings::GraphName) => {
            if !matches!(quantity_settings, QuantitySettings::GraphId {}) {
                return Err(eyre!(
                    "Observable '{}' uses GraphName labels, but its quantity is not graph_id.",
                    observable_name
                ));
            }
            let info = process_info.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses GraphName labels but no process information is available.",
                    observable_name
                )
            })?;
            info.graph_names.iter().cloned().map(Some).collect()
        }
        Some(DiscreteBinLabelsSettings::GraphGroupMasterName) => {
            if !matches!(quantity_settings, QuantitySettings::GraphGroupId {}) {
                return Err(eyre!(
                    "Observable '{}' uses GraphGroupMasterName labels, but its quantity is not graph_group_id.",
                    observable_name
                ));
            }
            let info = process_info.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses GraphGroupMasterName labels but no process information is available.",
                    observable_name
                )
            })?;
            info.graph_group_master_names
                .iter()
                .cloned()
                .map(Some)
                .collect()
        }
        Some(DiscreteBinLabelsSettings::Orientation) => {
            if !matches!(quantity_settings, QuantitySettings::OrientationId {}) {
                return Err(eyre!(
                    "Observable '{}' uses Orientation labels, but its quantity is not orientation_id.",
                    observable_name
                ));
            }
            let info = process_info.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses Orientation labels but no process information is available.",
                    observable_name
                )
            })?;
            let group_id = graph_group_context.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses Orientation labels but no singleton graph-group context could be resolved from its selections.",
                    observable_name
                )
            })?;
            info.orientation_labels_by_group
                .get(group_id)
                .ok_or_else(|| {
                    eyre!(
                        "Observable '{}' resolved graph group {} for Orientation labels, but the process only exposes {} graph groups.",
                        observable_name,
                        group_id,
                        info.orientation_labels_by_group.len()
                    )
                })?
                .iter()
                .cloned()
                .map(Some)
                .collect()
        }
        Some(DiscreteBinLabelsSettings::LmbChannelEdgeIds) => {
            if !matches!(quantity_settings, QuantitySettings::LmbChannelId {}) {
                return Err(eyre!(
                    "Observable '{}' uses LmbChannelEdgeIds labels, but its quantity is not lmb_channel_id.",
                    observable_name
                ));
            }
            let info = process_info.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses LmbChannelEdgeIds labels but no process information is available.",
                    observable_name
                )
            })?;
            let group_id = graph_group_context.ok_or_else(|| {
                eyre!(
                    "Observable '{}' uses LmbChannelEdgeIds labels but no singleton graph-group context could be resolved from its selections.",
                    observable_name
                )
            })?;
            info.lmb_channel_labels_by_group
                .get(group_id)
                .ok_or_else(|| {
                    eyre!(
                        "Observable '{}' resolved graph group {} for LmbChannelEdgeIds labels, but the process only exposes {} graph groups.",
                        observable_name,
                        group_id,
                        info.lmb_channel_labels_by_group.len()
                    )
                })?
                .iter()
                .cloned()
                .map(Some)
                .collect()
        }
    };
    if bin_labels.len() != n_bins {
        return Err(eyre!(
            "Observable '{}' resolved {} labels for {} discrete bins.",
            observable_name,
            bin_labels.len(),
            n_bins
        ));
    }
    Ok((min_bin_id, max_bin_id, bin_labels))
}

fn resolve_graph_group_context(
    observable_name: &str,
    selection_names: &[String],
    selector_settings: &SelectorsSettings,
    quantities: &QuantitiesSettings,
    process_info: Option<&HistogramProcessInfo>,
) -> Result<Option<usize>> {
    let mut resolved = None;
    for selection_name in selection_names {
        let Some(selector) = selector_settings.get(selection_name) else {
            return Err(eyre!(
                "Observable '{}' references unknown selector '{}'",
                observable_name,
                selection_name
            ));
        };
        let quantity = quantities.get(&selector.quantity).ok_or_else(|| {
            eyre!(
                "Observable '{}' selection '{}' references unknown quantity '{}'",
                observable_name,
                selection_name,
                selector.quantity
            )
        })?;
        let candidate = graph_group_context_candidate(
            observable_name,
            selection_name,
            selector,
            quantity,
            process_info,
        )?;
        if let Some(candidate) = candidate {
            if let Some(existing) = resolved
                && existing != candidate
            {
                return Err(eyre!(
                    "Observable '{}' resolves conflicting graph-group contexts {} and {} from its selections.",
                    observable_name,
                    existing,
                    candidate
                ));
            }
            resolved = Some(candidate);
        }
    }

    if resolved.is_none()
        && let Some(info) = process_info
        && info.graph_group_master_names.len() == 1
    {
        resolved = Some(0);
    }

    Ok(resolved)
}

fn graph_group_context_candidate(
    observable_name: &str,
    selection_name: &str,
    selector: &SelectorSettings,
    quantity: &QuantitySettings,
    process_info: Option<&HistogramProcessInfo>,
) -> Result<Option<usize>> {
    let SelectorDefinitionSettings::DiscreteRange(DiscreteRangeSelectorSettings {
        min: Some(min),
        max: Some(max),
    }) = &selector.selector
    else {
        return Ok(None);
    };
    if min != max {
        return Ok(None);
    }

    match quantity {
        QuantitySettings::GraphGroupId {} => {
            let group_id = singleton_selector_index(
                observable_name,
                selection_name,
                "graph_group_id",
                *min,
                process_info.map(|info| info.graph_group_master_names.len()),
            )?;
            Ok(Some(group_id))
        }
        QuantitySettings::GraphId {} => {
            let info = process_info.ok_or_else(|| {
                eyre!(
                    "Observable '{}' needs process information to map selector '{}' graph_id={} to a graph group.",
                    observable_name,
                    selection_name,
                    min
                )
            })?;
            let graph_id = singleton_selector_index(
                observable_name,
                selection_name,
                "graph_id",
                *min,
                Some(info.graph_to_group_id.len()),
            )?;
            Ok(info.graph_to_group_id.get(graph_id).copied())
        }
        _ => Ok(None),
    }
}

fn singleton_selector_index(
    observable_name: &str,
    selection_name: &str,
    quantity_label: &str,
    value: isize,
    upper_bound: Option<usize>,
) -> Result<usize> {
    let index = usize::try_from(value).map_err(|_| {
        eyre!(
            "Observable '{}' resolves selector '{}' {}={} as graph-group context, but the value must be non-negative.",
            observable_name,
            selection_name,
            quantity_label,
            value
        )
    })?;
    if let Some(upper_bound) = upper_bound
        && index >= upper_bound
    {
        return Err(eyre!(
            "Observable '{}' resolves selector '{}' {}={} as graph-group context, but the process only exposes {} valid values.",
            observable_name,
            selection_name,
            quantity_label,
            value,
            upper_bound
        ));
    }
    Ok(index)
}

impl Observable for HistogramObservable {
    fn process_event_groups<T: FloatLike>(
        &mut self,
        event_groups: &GenericEventGroupList<T>,
        clustering_registry: &CompiledClusteringRegistry,
    ) {
        self.state.new_sample_count += 1;
        self.touched_bins.clear();

        for event_group in event_groups.iter() {
            self.build_group_contributions(event_group, clustering_registry);
            self.mitigate_group_misbinning();

            for contribution in &self.pending_contributions {
                let grouped = &mut self.grouped_contributions[contribution.bin_index];
                if grouped.projected_weight == 0.0
                    && grouped.entry_count == 0
                    && grouped.mitigated_fill_count == 0
                {
                    self.touched_bins.push(contribution.bin_index);
                }
                grouped.projected_weight += contribution.projected_weight;
                grouped.entry_count += 1;
                if contribution.mitigated {
                    grouped.mitigated_fill_count += 1;
                }
            }
        }

        self.flush_sample_contributions();
        self.clear_sample_contributions();
    }

    fn merge_samples(&mut self, other: &mut HistogramObservable) -> Result<()> {
        self.state.merge_samples(&mut other.state)
    }

    fn update_result(&mut self, _iter: usize) {
        self.state.update_result();
    }
}

#[derive(Debug, Clone)]
pub enum Observables {
    Histogram(HistogramObservable),
}

impl Observables {
    #[allow(clippy::too_many_arguments)]
    fn from_settings(
        name: &str,
        settings: &ObservableSettings,
        quantities: &QuantitiesSettings,
        selector_settings: &SelectorsSettings,
        selector_registry: &BTreeMap<String, ConfiguredSelector>,
        process_info: Option<&HistogramProcessInfo>,
        clustering_registry: &mut CompiledClusteringRegistry,
        model: Option<&Model>,
    ) -> Result<Self> {
        let quantity = quantities.get(&settings.quantity).ok_or_else(|| {
            eyre!(
                "Observable '{}' references unknown quantity '{}'",
                name,
                settings.quantity
            )
        })?;
        let selections = settings
            .selections
            .iter()
            .map(|selector_name| {
                selector_registry
                    .get(selector_name)
                    .cloned()
                    .ok_or_else(|| {
                        eyre!(
                            "Observable '{}' references unknown selector '{}'",
                            name,
                            selector_name
                        )
                    })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Observables::Histogram(HistogramObservable::new(
            ObservableDefinition::from_settings(quantity, clustering_registry, model)?,
            selections,
            quantities,
            quantity,
            name,
            settings,
            selector_settings,
            process_info,
        )?))
    }

    pub(crate) fn process_event_groups<T: FloatLike>(
        &mut self,
        event_groups: &GenericEventGroupList<T>,
        clustering_registry: &CompiledClusteringRegistry,
    ) {
        match self {
            Observables::Histogram(observable) => {
                observable.process_event_groups(event_groups, clustering_registry)
            }
        }
    }

    pub(crate) fn merge_samples(&mut self, other: &mut Observables) -> Result<()> {
        match (self, other) {
            (Observables::Histogram(lhs), Observables::Histogram(rhs)) => lhs.merge_samples(rhs),
        }
    }

    pub(crate) fn required_clustering_handles(&self) -> Vec<ClusteringHandle> {
        match self {
            Observables::Histogram(observable) => observable.required_clustering_handles(),
        }
    }

    pub(crate) fn update_result(&mut self, iter: usize) {
        match self {
            Observables::Histogram(observable) => observable.update_result(iter),
        }
    }

    pub(crate) fn restore_snapshot(&mut self, snapshot: &HistogramSnapshot) -> Result<()> {
        match self {
            Observables::Histogram(observable) => observable.restore_snapshot(snapshot),
        }
    }

    pub(crate) fn cleared_clone(&self) -> Self {
        match self {
            Observables::Histogram(observable) => Observables::Histogram(HistogramObservable {
                definition: observable.definition.clone(),
                selections: observable.selections.clone(),
                entry_selection: observable.entry_selection,
                entry_index: observable.entry_index,
                misbinning_max_normalized_distance: observable.misbinning_max_normalized_distance,
                state: observable.state.cleared_clone(),
                pending_contributions: Vec::with_capacity(
                    observable.pending_contributions.capacity(),
                ),
                candidate_pairs: Vec::with_capacity(observable.candidate_pairs.capacity()),
                used_contributions: Vec::with_capacity(observable.used_contributions.capacity()),
                grouped_contributions: vec![
                    GroupBinContribution::default();
                    observable.grouped_contributions.len()
                ],
                grouped_underflow: GroupBinContribution::default(),
                grouped_overflow: GroupBinContribution::default(),
                touched_bins: Vec::with_capacity(observable.touched_bins.capacity()),
            }),
        }
    }

    pub(crate) fn snapshot(&self) -> HistogramSnapshot {
        match self {
            Observables::Histogram(observable) => observable.snapshot(),
        }
    }

    pub(crate) fn accumulator_state(&self) -> HistogramAccumulatorState {
        match self {
            Observables::Histogram(observable) => observable.state.clone(),
        }
    }

    pub(crate) fn merge_accumulator_state(
        &mut self,
        other: &mut HistogramAccumulatorState,
    ) -> Result<()> {
        match self {
            Observables::Histogram(observable) => observable.state.merge_samples(other),
        }
    }
}

impl HistogramObservable {
    fn restore_snapshot(&mut self, snapshot: &HistogramSnapshot) -> Result<()> {
        let restored_state = snapshot.clone().into_accumulator_state();
        self.state.ensure_compatible(&restored_state)?;
        self.state = restored_state;
        self.pending_contributions.clear();
        self.candidate_pairs.clear();
        self.used_contributions.clear();
        self.grouped_contributions
            .fill(GroupBinContribution::default());
        self.grouped_underflow = GroupBinContribution::default();
        self.grouped_overflow = GroupBinContribution::default();
        self.touched_bins.clear();
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
enum HistogramBinPosition {
    Underflow,
    Overflow,
    InRange(usize),
}

fn histogram_bin_position(
    value: f64,
    x_min: f64,
    x_max: f64,
    num_bins: usize,
) -> Option<HistogramBinPosition> {
    if num_bins == 0 || x_max <= x_min || !value.is_finite() {
        return None;
    }

    if value < x_min {
        return Some(HistogramBinPosition::Underflow);
    }
    if value >= x_max {
        return Some(HistogramBinPosition::Overflow);
    }

    let index = ((value - x_min) / (x_max - x_min) * num_bins as f64) as usize;
    if index < num_bins {
        Some(HistogramBinPosition::InRange(index))
    } else {
        Some(HistogramBinPosition::Overflow)
    }
}

fn histogram_bin_position_discrete(
    value: isize,
    min_bin_id: isize,
    num_bins: usize,
) -> HistogramBinPosition {
    let max_bin_id = min_bin_id + num_bins as isize - 1;
    if value < min_bin_id {
        HistogramBinPosition::Underflow
    } else if value > max_bin_id {
        HistogramBinPosition::Overflow
    } else {
        HistogramBinPosition::InRange((value - min_bin_id) as usize)
    }
}

fn sort_discrete_snapshot_bins(
    bins: &mut [HistogramBinSnapshot],
    ordering: DiscreteBinOrdering,
    sample_count: usize,
) {
    bins.sort_by(|lhs, rhs| {
        let base = match ordering {
            DiscreteBinOrdering::AscendingBinId => Ordering::Equal,
            DiscreteBinOrdering::ValueDescending => rhs
                .average(sample_count)
                .partial_cmp(&lhs.average(sample_count))
                .unwrap_or(Ordering::Equal),
            DiscreteBinOrdering::AbsValueDescending => rhs
                .average(sample_count)
                .abs()
                .partial_cmp(&lhs.average(sample_count).abs())
                .unwrap_or(Ordering::Equal),
        };
        if base == Ordering::Equal {
            lhs.bin_id
                .unwrap_or_default()
                .cmp(&rhs.bin_id.unwrap_or_default())
        } else {
            base
        }
    });
}

fn unit_complex<T: FloatLike>(reference: &F<T>) -> Complex<F<T>> {
    Complex::new(reference.one(), reference.zero())
}

fn apply_entry_selection<T: FloatLike>(
    entries: &[ObservableEntry<T>],
    selection: EntrySelection,
    entry_index: usize,
) -> ObservableEntries<T> {
    match selection {
        EntrySelection::All => entries.iter().cloned().collect(),
        EntrySelection::LeadingOnly => entries.first().into_iter().cloned().collect(),
        EntrySelection::NthOnly => entries.get(entry_index).into_iter().cloned().collect(),
    }
}

fn transform_value<T: FloatLike>(value: &F<T>, value_transform: ObservableValueTransform) -> F<T> {
    match value_transform {
        ObservableValueTransform::Identity => value.clone(),
        ObservableValueTransform::Log10 => value.log10(),
    }
}

fn value_in_range(value: f64, min: Option<f64>, max: Option<f64>) -> bool {
    if let Some(min) = min
        && value < min
    {
        return false;
    }

    if let Some(max) = max {
        value <= max
    } else {
        true
    }
}

fn discrete_value_in_range(value: isize, min: Option<isize>, max: Option<isize>) -> bool {
    if let Some(min) = min
        && value < min
    {
        return false;
    }

    if let Some(max) = max {
        value <= max
    } else {
        true
    }
}

fn combined_entry_weight<T: FloatLike>(
    event: &GenericEvent<T>,
    entry: &ObservableEntry<T>,
) -> Complex<F<T>> {
    event.weight.clone() * entry.weight_modifier.clone()
}

impl ObservablePhase {
    fn project<T: FloatLike>(&self, weight: &Complex<F<T>>) -> F<T> {
        match self {
            ObservablePhase::Real => weight.re.clone(),
            ObservablePhase::Imag => weight.im.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        ContinuousHistogramSettings, DiscreteBinDomainSettings, DiscreteBinOrdering,
        DiscreteRangeSelectorSettings, EntrySelection, HistogramAccumulatorState,
        HistogramBinSnapshot, HistogramProcessInfo, HistogramSettings, HistogramSnapshot,
        HistogramSnapshotKind, HistogramStatisticsSnapshot, JetClusteringSettings, ObservablePhase,
        ObservableSettings, ObservableValueTransform, QuantitySettings, SelectorDefinitionSettings,
        SelectorSettings, resolve_graph_group_context,
    };
    use schemars::schema_for;
    use serde_json::Value as JsonValue;
    use std::collections::BTreeMap;
    use std::fs;

    fn sample_histogram_snapshot() -> HistogramSnapshot {
        HistogramSnapshot {
            kind: HistogramSnapshotKind::Continuous,
            title: "top_pt".to_string(),
            type_description: "AL".to_string(),
            phase: ObservablePhase::Real,
            value_transform: ObservableValueTransform::Identity,
            supports_misbinning_mitigation: true,
            x_min: Some(0.0),
            x_max: Some(10.0),
            sample_count: 3,
            log_x_axis: false,
            log_y_axis: true,
            discrete_min_bin_id: None,
            discrete_ordering: None,
            bins: vec![HistogramBinSnapshot {
                x_min: Some(0.0),
                x_max: Some(10.0),
                bin_id: None,
                label: None,
                entry_count: 3,
                sum_weights: 3.75,
                sum_weights_squared: 5.625,
                mitigated_fill_count: 2,
            }],
            underflow_bin: HistogramBinSnapshot {
                x_min: None,
                x_max: Some(0.0),
                bin_id: None,
                label: None,
                entry_count: 1,
                sum_weights: 1.5,
                sum_weights_squared: 2.25,
                mitigated_fill_count: 0,
            },
            overflow_bin: HistogramBinSnapshot {
                x_min: Some(10.0),
                x_max: None,
                bin_id: None,
                label: None,
                entry_count: 2,
                sum_weights: 0.5,
                sum_weights_squared: 0.25,
                mitigated_fill_count: 1,
            },
            statistics: HistogramStatisticsSnapshot {
                in_range_entry_count: 3,
                nan_value_count: 3,
                mitigated_pair_count: 4,
            },
        }
    }

    #[test]
    fn jet_clustering_settings_normalize_explicit_clustered_pdgs() {
        let settings = JetClusteringSettings {
            clustered_pdgs: Some(vec![21, -1, 21, 1, -1]),
            ..JetClusteringSettings::default()
        };
        let resolved = settings
            .resolve(None)
            .expect("explicit clustered_pdgs should not require a model");
        assert_eq!(resolved.clustered_pdgs, vec![-1, 1, 21]);
    }

    #[test]
    fn jet_clustering_settings_require_model_for_default_clustered_pdgs() {
        let error = JetClusteringSettings::default()
            .resolve(None)
            .expect_err("default clustered_pdgs should require model context");
        assert!(
            error
                .to_string()
                .contains("Cannot resolve default clustered_pdgs without a model")
        );
    }

    #[test]
    fn histogram_snapshot_json_round_trip() {
        let snapshot = sample_histogram_snapshot();

        let file_path = std::env::temp_dir().join(format!(
            "gammaloop_histogram_snapshot_{}_{}.json",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));

        snapshot.to_json_file(&file_path).unwrap();
        let loaded = HistogramSnapshot::from_json_file(&file_path).unwrap();
        let _ = fs::remove_file(&file_path);

        assert_eq!(loaded, snapshot);
    }

    #[test]
    fn histogram_snapshot_rescale_scales_weight_sums_only() {
        let snapshot = sample_histogram_snapshot();
        let scaled = snapshot.rescaled(2.0);

        assert_eq!(scaled.sample_count, snapshot.sample_count);
        assert_eq!(scaled.statistics, snapshot.statistics);
        assert_eq!(scaled.bins[0].entry_count, snapshot.bins[0].entry_count);
        assert_eq!(
            scaled.bins[0].mitigated_fill_count,
            snapshot.bins[0].mitigated_fill_count
        );
        assert_eq!(scaled.bins[0].sum_weights, 7.5);
        assert_eq!(scaled.bins[0].sum_weights_squared, 22.5);
        assert_eq!(scaled.underflow_bin.sum_weights, 3.0);
        assert_eq!(scaled.underflow_bin.sum_weights_squared, 9.0);
        assert_eq!(scaled.overflow_bin.sum_weights, 1.0);
        assert_eq!(scaled.overflow_bin.sum_weights_squared, 1.0);
        assert_eq!(
            scaled.bins[0].average(scaled.sample_count),
            snapshot.bins[0].average(snapshot.sample_count) * 2.0
        );
        assert_eq!(
            scaled.bins[0].error(scaled.sample_count),
            snapshot.bins[0].error(snapshot.sample_count) * 2.0
        );
    }

    #[test]
    fn histogram_snapshot_hwu_header_uses_type_description() {
        let mut snapshot = sample_histogram_snapshot();
        snapshot.type_description = "SB".to_string();

        let mut output = Vec::new();
        snapshot.write_hwu_block(&mut output).unwrap();
        let output = String::from_utf8(output).unwrap();

        assert!(output.contains("|TYPE@SB\""), "{output}");
    }

    #[test]
    fn observable_settings_infer_discrete_histogram_from_domain() {
        let settings: ObservableSettings = toml::from_str(
            r#"
quantity = "integral"
domain = { type = "single_bin" }
"#,
        )
        .unwrap();

        let HistogramSettings::Discrete(histogram) = settings.histogram else {
            panic!("expected discrete histogram");
        };
        assert_eq!(histogram.domain, DiscreteBinDomainSettings::SingleBin);
    }

    #[test]
    fn observable_settings_serialize_flat_histogram_fields() {
        let settings = ObservableSettings {
            quantity: "integral".to_string(),
            selections: vec!["graph_cut".to_string()],
            entry_selection: super::EntrySelection::All,
            entry_index: 0,
            value_transform: ObservableValueTransform::Identity,
            phase: ObservablePhase::Real,
            misbinning_max_normalized_distance: None,
            histogram: HistogramSettings::Continuous(ContinuousHistogramSettings {
                x_min: 0.0,
                x_max: 1.0,
                n_bins: 1,
                log_x_axis: false,
                log_y_axis: true,
                title: Some("xs".to_string()),
                type_description: "AL".to_string(),
            }),
        };

        let serialized = toml::to_string(&settings).unwrap();
        assert!(serialized.contains("kind = \"continuous\""));
        assert!(serialized.contains("quantity = \"integral\""));
        assert!(serialized.contains("n_bins = 1"));
        assert!(!serialized.contains("histogram"));
    }

    #[test]
    fn observable_settings_schema_is_flat() {
        let schema = serde_json::to_value(schema_for!(ObservableSettings)).unwrap();
        let properties = schema
            .get("properties")
            .and_then(JsonValue::as_object)
            .expect("observable schema properties");
        assert!(properties.contains_key("kind"));
        assert!(properties.contains_key("quantity"));
        assert!(properties.contains_key("domain"));
        assert!(!properties.contains_key("histogram"));
    }

    #[test]
    fn selector_settings_schema_is_flat() {
        let schema = serde_json::to_string(&schema_for!(SelectorSettings)).unwrap();
        assert!(schema.contains("\"selector\""));
        assert!(schema.contains("\"quantity\""));
        assert!(!schema.contains("\"histogram\""));
    }

    #[test]
    fn selector_settings_serialize_flat_discrete_range() {
        let settings = SelectorSettings {
            quantity: "graph_id".to_string(),
            active: false,
            entry_selection: super::EntrySelection::All,
            entry_index: 0,
            selector: SelectorDefinitionSettings::DiscreteRange(DiscreteRangeSelectorSettings {
                min: Some(3),
                max: Some(3),
            }),
        };

        let serialized = toml::to_string(&settings).unwrap();
        assert!(serialized.contains("selector = \"discrete_range\""));
        assert!(serialized.contains("active = false"));
        assert!(serialized.contains("min = 3"));
        assert!(!serialized.contains("selector = {"));
    }

    #[test]
    fn selector_settings_deserialize_flat_discrete_range() {
        let settings: SelectorSettings = toml::from_str(
            r#"
quantity = "graph_id"
active = false
selector = "discrete_range"
min = 0
max = 0
"#,
        )
        .unwrap();

        assert_eq!(settings.quantity, "graph_id");
        assert!(!settings.active);
        match settings.selector {
            SelectorDefinitionSettings::DiscreteRange(selector) => {
                assert_eq!(selector.min, Some(0));
                assert_eq!(selector.max, Some(0));
            }
            other => panic!("expected discrete-range selector, got {other:?}"),
        }
    }

    #[test]
    fn discrete_histogram_ordering_changes_presentation_only() {
        let mut histogram = HistogramAccumulatorState::discrete(
            "jet_count".to_string(),
            "AL".to_string(),
            ObservablePhase::Real,
            0,
            1,
            DiscreteBinOrdering::AscendingBinId,
            true,
            vec![Some("zero".to_string()), Some("one".to_string())],
        )
        .unwrap();
        histogram
            .fill_discrete_sample(&[(0, 1.0), (1, 3.0)])
            .unwrap();
        histogram.update_results();

        let ascending = histogram.snapshot();
        assert_eq!(ascending.bins[0].bin_id, Some(0));
        assert_eq!(ascending.bins[1].bin_id, Some(1));
        assert_eq!(ascending.discrete_min_bin_id, Some(0));

        histogram
            .change_bin_ordering(DiscreteBinOrdering::ValueDescending)
            .unwrap();
        let reordered = histogram.snapshot();
        assert_eq!(reordered.bins[0].bin_id, Some(1));
        assert_eq!(reordered.bins[1].bin_id, Some(0));
        assert_eq!(reordered.bins[0].label.as_deref(), Some("one"));
        assert_eq!(reordered.bins[1].label.as_deref(), Some("zero"));
        assert_eq!(reordered.rebin(99).unwrap(), reordered);
    }

    #[test]
    fn discrete_histogram_constructor_rejects_invalid_range() {
        let err = HistogramAccumulatorState::discrete(
            "bad_range".to_string(),
            "AL".to_string(),
            ObservablePhase::Real,
            3,
            1,
            DiscreteBinOrdering::AscendingBinId,
            true,
            Vec::new(),
        )
        .unwrap_err();

        assert!(err.to_string().contains("invalid discrete histogram range"));
    }

    #[test]
    fn resolve_graph_group_context_uses_selector_quantity_type() {
        let quantities = BTreeMap::from([("my_graph".to_string(), QuantitySettings::GraphId {})]);
        let selectors = BTreeMap::from([(
            "graph_only".to_string(),
            SelectorSettings {
                quantity: "my_graph".to_string(),
                active: false,
                entry_selection: EntrySelection::All,
                entry_index: 0,
                selector: SelectorDefinitionSettings::DiscreteRange(
                    DiscreteRangeSelectorSettings {
                        min: Some(1),
                        max: Some(1),
                    },
                ),
            },
        )]);
        let process_info = HistogramProcessInfo {
            graph_names: vec!["g0".to_string(), "g1".to_string()],
            graph_to_group_id: vec![0, 2],
            graph_group_master_names: vec![
                "group0".to_string(),
                "group1".to_string(),
                "group2".to_string(),
            ],
            orientation_labels_by_group: vec![Vec::new(), Vec::new(), Vec::new()],
            lmb_channel_labels_by_group: vec![Vec::new(), Vec::new(), Vec::new()],
        };

        let resolved = resolve_graph_group_context(
            "obs",
            &["graph_only".to_string()],
            &selectors,
            &quantities,
            Some(&process_info),
        )
        .unwrap();

        assert_eq!(resolved, Some(2));
    }

    #[test]
    fn resolve_graph_group_context_rejects_negative_singleton_values() {
        let quantities =
            BTreeMap::from([("my_group".to_string(), QuantitySettings::GraphGroupId {})]);
        let selectors = BTreeMap::from([(
            "negative_group".to_string(),
            SelectorSettings {
                quantity: "my_group".to_string(),
                active: false,
                entry_selection: EntrySelection::All,
                entry_index: 0,
                selector: SelectorDefinitionSettings::DiscreteRange(
                    DiscreteRangeSelectorSettings {
                        min: Some(-1),
                        max: Some(-1),
                    },
                ),
            },
        )]);
        let process_info = HistogramProcessInfo {
            graph_names: Vec::new(),
            graph_to_group_id: Vec::new(),
            graph_group_master_names: vec!["group0".to_string()],
            orientation_labels_by_group: vec![Vec::new()],
            lmb_channel_labels_by_group: vec![Vec::new()],
        };

        let err = resolve_graph_group_context(
            "obs",
            &["negative_group".to_string()],
            &selectors,
            &quantities,
            Some(&process_info),
        )
        .unwrap_err();

        assert!(err.to_string().contains("must be non-negative"));
    }
}
