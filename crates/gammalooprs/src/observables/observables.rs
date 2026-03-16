use super::clustering::{JetAlgorithm, JetClustering};
use super::events::{GenericEvent, GenericEventGroup, GenericEventGroupList};
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
use tracing::info;

pub type QuantitiesSettings = BTreeMap<String, QuantitySettings>;
pub type ObservablesSettings = BTreeMap<String, ObservableSettings>;
pub type SelectorsSettings = BTreeMap<String, SelectorSettings>;

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
}

impl fmt::Display for FilterQuantity {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            FilterQuantity::Energy => write!(f, "Energy"),
            FilterQuantity::CosThetaP => write!(f, "CosTheta"),
            FilterQuantity::PT => write!(f, "Pt"),
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
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
pub struct ParticleScalarQuantitySettings {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub pdgs: Vec<isize>,
    pub quantity: FilterQuantity,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
pub struct JetPtQuantitySettings {
    #[serde(flatten)]
    pub clustering: JetClusteringSettings,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[serde(tag = "type", rename_all = "snake_case")]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
pub enum QuantitySettings {
    ParticleScalar(ParticleScalarQuantitySettings),
    JetPt(JetPtQuantitySettings),
    AFB {},
    CrossSection {},
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
pub struct ValueRangeSelectorSettings {
    pub min: f64,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub max: Option<f64>,
    #[serde(default, skip_serializing_if = "is_default_selector_reduction")]
    pub reduction: SelectorReduction,
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
    CountRange(CountRangeSelectorSettings),
}

#[derive(Debug, Clone, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
pub struct SelectorSettings {
    pub quantity: String,
    pub entry_selection: EntrySelection,
    pub entry_index: usize,
    pub selector: SelectorDefinitionSettings,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum SelectorSerdeTag {
    ValueRange,
    CountRange,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ValueRangeSelectorSettingsSerde {
    quantity: String,
    #[serde(default, skip_serializing_if = "is_default_entry_selection")]
    entry_selection: EntrySelection,
    #[serde(default, skip_serializing_if = "is_default_entry_index")]
    entry_index: usize,
    selector: SelectorSerdeTag,
    min: f64,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    max: Option<f64>,
    #[serde(default, skip_serializing_if = "is_default_selector_reduction")]
    reduction: SelectorReduction,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct CountRangeSelectorSettingsSerde {
    quantity: String,
    #[serde(default, skip_serializing_if = "is_default_entry_selection")]
    entry_selection: EntrySelection,
    #[serde(default, skip_serializing_if = "is_default_entry_index")]
    entry_index: usize,
    selector: SelectorSerdeTag,
    min_count: usize,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    max_count: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
enum SelectorSettingsSerde {
    ValueRange(ValueRangeSelectorSettingsSerde),
    CountRange(CountRangeSelectorSettingsSerde),
}

impl Serialize for SelectorSettings {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match &self.selector {
            SelectorDefinitionSettings::ValueRange(selector) => ValueRangeSelectorSettingsSerde {
                quantity: self.quantity.clone(),
                entry_selection: self.entry_selection,
                entry_index: self.entry_index,
                selector: SelectorSerdeTag::ValueRange,
                min: selector.min,
                max: selector.max,
                reduction: selector.reduction,
            }
            .serialize(serializer),
            SelectorDefinitionSettings::CountRange(selector) => CountRangeSelectorSettingsSerde {
                quantity: self.quantity.clone(),
                entry_selection: self.entry_selection,
                entry_index: self.entry_index,
                selector: SelectorSerdeTag::CountRange,
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
            SelectorSettingsSerde::ValueRange(settings) => {
                let SelectorSerdeTag::ValueRange = settings.selector else {
                    return Err(serde::de::Error::custom(
                        "selector = \"value_range\" is required for a value-range selector",
                    ));
                };
                SelectorSettings {
                    quantity: settings.quantity,
                    entry_selection: settings.entry_selection,
                    entry_index: settings.entry_index,
                    selector: SelectorDefinitionSettings::ValueRange(ValueRangeSelectorSettings {
                        min: settings.min,
                        max: settings.max,
                        reduction: settings.reduction,
                    }),
                }
            }
            SelectorSettingsSerde::CountRange(settings) => {
                let SelectorSerdeTag::CountRange = settings.selector else {
                    return Err(serde::de::Error::custom(
                        "selector = \"count_range\" is required for a count-range selector",
                    ));
                };
                SelectorSettings {
                    quantity: settings.quantity,
                    entry_selection: settings.entry_selection,
                    entry_index: settings.entry_index,
                    selector: SelectorDefinitionSettings::CountRange(CountRangeSelectorSettings {
                        min_count: settings.min_count,
                        max_count: settings.max_count,
                    }),
                }
            }
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(default, deny_unknown_fields)]
pub struct HistogramSettings {
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
}

impl Default for HistogramSettings {
    fn default() -> Self {
        Self {
            x_min: 0.0,
            x_max: 0.0,
            n_bins: 0,
            log_x_axis: false,
            log_y_axis: true,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
pub struct ObservableSettings {
    pub quantity: String,
    #[serde(default, skip_serializing_if = "is_default_entry_selection")]
    pub entry_selection: EntrySelection,
    #[serde(default, skip_serializing_if = "is_default_entry_index")]
    pub entry_index: usize,
    #[serde(default, skip_serializing_if = "is_default_value_transform")]
    pub value_transform: ObservableValueTransform,
    #[serde(default, skip_serializing_if = "is_default_observable_phase")]
    pub phase: ObservablePhase,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub misbinning_max_normalized_distance: Option<f64>,
    #[serde(flatten)]
    pub histogram: HistogramSettings,
}

fn default_true() -> bool {
    true
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

#[derive(Debug, Clone)]
pub struct ObservableEntry<T: FloatLike> {
    pub value: F<T>,
    pub weight_modifier: Complex<F<T>>,
}

impl<T: FloatLike> ObservableEntry<T> {
    fn unit(value: F<T>) -> Self {
        Self {
            weight_modifier: unit_complex(&value),
            value,
        }
    }
}

type ObservableEntries<T> = SmallVec<[ObservableEntry<T>; 8]>;

#[derive(Debug, Clone)]
enum SelectorCriterion {
    ValueRange {
        selection: EntrySelection,
        entry_index: usize,
        reduction: SelectorReduction,
        min: f64,
        max: Option<f64>,
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
                    SelectorReduction::AnyInRange => selected
                        .iter()
                        .any(|entry| value_in_range(entry.value.into_ff64().0, *min, *max)),
                    SelectorReduction::AllInRange => selected
                        .iter()
                        .all(|entry| value_in_range(entry.value.into_ff64().0, *min, *max)),
                }
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
struct CrossSectionDefinition;

impl CrossSectionDefinition {
    fn process_event<T: FloatLike>(&mut self, event: &GenericEvent<T>) -> ObservableEntries<T> {
        let reference = event.weight.re.clone();
        let half = reference.one() / reference.from_usize(2);
        smallvec![ObservableEntry::unit(half)]
    }
}

#[derive(Debug, Clone)]
struct ParticleScalarDefinition {
    pdgs: Vec<isize>,
    quantity: FilterQuantity,
}

impl ParticleScalarDefinition {
    fn new(quantity: FilterQuantity, pdgs: Vec<isize>) -> Self {
        Self { pdgs, quantity }
    }

    fn process_event<T: FloatLike>(&mut self, event: &GenericEvent<T>) -> ObservableEntries<T> {
        let incoming_beam = event.kinematic_configuration.0.get(1).cloned();
        let mut entries = ObservableEntries::new();

        for (pdg, momentum) in event
            .cut_info
            .particle_pdgs
            .1
            .iter()
            .copied()
            .zip(event.kinematic_configuration.1.iter().cloned())
        {
            if !self.pdgs.contains(&pdg) {
                continue;
            }

            let Some(value) = (match self.quantity {
                FilterQuantity::Energy => Some(momentum.temporal.value),
                FilterQuantity::PT => Some(momentum.pt()),
                FilterQuantity::CosThetaP => incoming_beam.clone().map(|beam| {
                    let beam_spatial = beam.spatial.clone();
                    let momentum_spatial = momentum.spatial.clone();
                    beam_spatial.clone() * momentum_spatial.clone()
                        / (beam_spatial.norm() * momentum_spatial.norm())
                }),
            }) else {
                continue;
            };

            entries.push(ObservableEntry::unit(value));
        }

        entries
    }
}

#[derive(Debug, Clone, Default)]
struct ForwardBackwardDefinition;

impl ForwardBackwardDefinition {
    fn process_event<T: FloatLike>(&mut self, event: &GenericEvent<T>) -> ObservableEntries<T> {
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
        smallvec![ObservableEntry::unit(
            beam_spatial.clone() * outgoing_spatial.clone()
                / (beam_spatial.norm() * outgoing_spatial.norm()),
        )]
    }
}

#[derive(Debug, Clone)]
struct JetPtDefinition {
    clustering: JetClustering,
}

impl JetPtDefinition {
    fn new(settings: &JetClusteringSettings) -> Self {
        Self {
            clustering: JetClustering::new(settings.algorithm, settings.dR, settings.min_jpt),
        }
    }

    fn process_event<T: FloatLike>(&mut self, event: &GenericEvent<T>) -> ObservableEntries<T> {
        let jets = self.clustering.process_event(event);
        let mut entries = ObservableEntries::new();
        for jet in jets.jets {
            entries.push(ObservableEntry::unit(jet.pt()));
        }
        entries
    }
}

#[derive(Debug, Clone)]
enum ObservableDefinition {
    CrossSection(CrossSectionDefinition),
    ParticleScalar(ParticleScalarDefinition),
    ForwardBackward(ForwardBackwardDefinition),
    JetPt(JetPtDefinition),
}

impl ObservableDefinition {
    fn from_settings(settings: &QuantitySettings) -> Self {
        match settings {
            QuantitySettings::ParticleScalar(settings) => ObservableDefinition::ParticleScalar(
                ParticleScalarDefinition::new(settings.quantity, settings.pdgs.clone()),
            ),
            QuantitySettings::JetPt(settings) => {
                ObservableDefinition::JetPt(JetPtDefinition::new(&settings.clustering))
            }
            QuantitySettings::AFB {} => {
                ObservableDefinition::ForwardBackward(ForwardBackwardDefinition)
            }
            QuantitySettings::CrossSection {} => {
                ObservableDefinition::CrossSection(CrossSectionDefinition)
            }
        }
    }

    fn process_event<T: FloatLike>(&mut self, event: &GenericEvent<T>) -> ObservableEntries<T> {
        match self {
            ObservableDefinition::CrossSection(definition) => definition.process_event(event),
            ObservableDefinition::ParticleScalar(definition) => definition.process_event(event),
            ObservableDefinition::ForwardBackward(definition) => definition.process_event(event),
            ObservableDefinition::JetPt(definition) => definition.process_event(event),
        }
    }
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
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, JsonSchema)]
pub struct HistogramStatisticsSnapshot {
    pub in_range_entry_count: usize,
    pub nan_value_count: usize,
    pub mitigated_pair_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, JsonSchema)]
pub struct HistogramSnapshot {
    pub title: String,
    pub phase: ObservablePhase,
    pub value_transform: ObservableValueTransform,
    pub x_min: f64,
    pub x_max: f64,
    pub sample_count: usize,
    pub log_x_axis: bool,
    pub log_y_axis: bool,
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
            "<histogram> {} \"{} |X_AXIS@{} |Y_AXIS@{} |TYPE@AL\"",
            self.bins.len(),
            self.title,
            x_axis_mode,
            y_axis_mode,
        )?;

        for bin in &self.bins {
            writeln!(
                writer,
                "  {:.8e}   {:.8e}   {:.8e}   {:.8e}",
                bin.x_min.unwrap_or(self.x_min),
                bin.x_max.unwrap_or(self.x_max),
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
            || self.phase != other.phase
            || self.value_transform != other.value_transform
            || self.x_min != other.x_min
            || self.x_max != other.x_max
            || self.log_x_axis != other.log_x_axis
            || self.log_y_axis != other.log_y_axis
            || self.bins.len() != other.bins.len()
        {
            return Err(eyre!(
                "Cannot merge incompatible histogram snapshots '{}'",
                self.title
            ));
        }

        self.sample_count += other.sample_count;
        for (bin, other_bin) in self.bins.iter_mut().zip(other.bins.iter()) {
            bin.merge_in_place(other_bin);
        }
        self.underflow_bin.merge_in_place(&other.underflow_bin);
        self.overflow_bin.merge_in_place(&other.overflow_bin);
        self.statistics.in_range_entry_count += other.statistics.in_range_entry_count;
        self.statistics.nan_value_count += other.statistics.nan_value_count;
        self.statistics.mitigated_pair_count += other.statistics.mitigated_pair_count;
        Ok(())
    }

    pub fn rebin(&self, contiguous_bins: usize) -> Result<Self> {
        if contiguous_bins == 0 {
            return Err(eyre!("Rebinning factor must be strictly positive."));
        }
        if self.bins.is_empty() || self.bins.len() % contiguous_bins != 0 {
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
            phase: self.phase,
            value_transform: self.value_transform,
            x_min: self.x_min,
            x_max: self.x_max,
            sample_count: self.sample_count,
            log_x_axis: self.log_x_axis,
            log_y_axis: self.log_y_axis,
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
    pub title: String,
    pub phase: ObservablePhase,
    pub value_transform: ObservableValueTransform,
    pub x_min: f64,
    pub x_max: f64,
    pub sample_count: usize,
    new_sample_count: usize,
    pub log_x_axis: bool,
    pub log_y_axis: bool,
    pub bins: Vec<ObservableBinAccumulator>,
    pub underflow_bin: ObservableBinAccumulator,
    pub overflow_bin: ObservableBinAccumulator,
    pub statistics: ObservableHistogramStatistics,
}

impl HistogramAccumulatorState {
    fn new(
        title: String,
        phase: ObservablePhase,
        value_transform: ObservableValueTransform,
        x_min: f64,
        x_max: f64,
        log_x_axis: bool,
        log_y_axis: bool,
        n_bins: usize,
    ) -> Self {
        Self {
            title,
            phase,
            value_transform,
            x_min,
            x_max,
            sample_count: 0,
            new_sample_count: 0,
            log_x_axis,
            log_y_axis,
            bins: vec![ObservableBinAccumulator::default(); n_bins],
            underflow_bin: ObservableBinAccumulator::default(),
            overflow_bin: ObservableBinAccumulator::default(),
            statistics: ObservableHistogramStatistics::default(),
        }
    }

    fn cleared_clone(&self) -> Self {
        Self::new(
            self.title.clone(),
            self.phase,
            self.value_transform,
            self.x_min,
            self.x_max,
            self.log_x_axis,
            self.log_y_axis,
            self.bins.len(),
        )
    }

    fn snapshot(&self) -> HistogramSnapshot {
        let bins = self
            .bins
            .iter()
            .enumerate()
            .map(|(index, bin)| {
                let x_min =
                    (self.x_max - self.x_min) * index as f64 / self.bins.len() as f64 + self.x_min;
                let x_max = (self.x_max - self.x_min) * (index + 1) as f64 / self.bins.len() as f64
                    + self.x_min;
                HistogramBinSnapshot {
                    x_min: Some(x_min),
                    x_max: Some(x_max),
                    entry_count: bin.total_entry_count(),
                    sum_weights: bin.total_sum_weights(),
                    sum_weights_squared: bin.total_sum_weights_squared(),
                    mitigated_fill_count: bin.total_mitigated_fill_count(),
                }
            })
            .collect();

        HistogramSnapshot {
            title: self.title.clone(),
            phase: self.phase,
            value_transform: self.value_transform,
            x_min: self.x_min,
            x_max: self.x_max,
            sample_count: self.sample_count + self.new_sample_count,
            log_x_axis: self.log_x_axis,
            log_y_axis: self.log_y_axis,
            bins,
            underflow_bin: HistogramBinSnapshot {
                x_min: None,
                x_max: Some(self.x_min),
                entry_count: self.underflow_bin.total_entry_count(),
                sum_weights: self.underflow_bin.total_sum_weights(),
                sum_weights_squared: self.underflow_bin.total_sum_weights_squared(),
                mitigated_fill_count: self.underflow_bin.total_mitigated_fill_count(),
            },
            overflow_bin: HistogramBinSnapshot {
                x_min: Some(self.x_max),
                x_max: None,
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

    fn update_result(&mut self) {
        self.sample_count += self.new_sample_count;
        self.new_sample_count = 0;
        for bin in &mut self.bins {
            bin.update_iter();
        }
        self.underflow_bin.update_iter();
        self.overflow_bin.update_iter();
        self.statistics.update_iter();

        for (i, bin) in self.bins.iter().enumerate() {
            let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
            let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64;
            info!(
                "{}={}: {} +/- {}",
                c1,
                c2,
                bin.average(self.sample_count),
                bin.error(self.sample_count)
            );
        }

        info!(
            "{} stats: entries={}, underflow={}, overflow={}, nan_values={}, mitigated_pairs={}",
            self.title,
            self.statistics.in_range_entry_count,
            self.underflow_bin.total_entry_count(),
            self.overflow_bin.total_entry_count(),
            self.statistics.nan_value_count,
            self.statistics.mitigated_pair_count,
        );
    }

    fn write_hwu_block<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        let x_axis_mode = if self.log_x_axis { "LOG" } else { "LIN" };
        let y_axis_mode = if self.log_y_axis { "LOG" } else { "LIN" };

        writeln!(writer, "##& xmin & xmax & central value & dy &\n")?;
        writeln!(
            writer,
            "<histogram> {} \"{} |X_AXIS@{} |Y_AXIS@{} |TYPE@AL\"",
            self.bins.len(),
            self.title,
            x_axis_mode,
            y_axis_mode,
        )?;

        for (i, bin) in self.bins.iter().enumerate() {
            let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
            let c2 =
                (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64 + self.x_min;
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
        Self {
            title: snapshot.title.clone(),
            phase: snapshot.phase,
            value_transform: snapshot.value_transform,
            x_min: snapshot.x_min,
            x_max: snapshot.x_max,
            sample_count: snapshot.sample_count,
            new_sample_count: 0,
            log_x_axis: snapshot.log_x_axis,
            log_y_axis: snapshot.log_y_axis,
            bins: snapshot
                .bins
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

    fn ensure_compatible(&self, other: &HistogramAccumulatorState) -> Result<()> {
        if self.title != other.title
            || self.phase != other.phase
            || self.value_transform != other.value_transform
            || self.x_min != other.x_min
            || self.x_max != other.x_max
            || self.log_x_axis != other.log_x_axis
            || self.log_y_axis != other.log_y_axis
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

pub trait EventSelector {
    fn process_event<T: FloatLike>(&mut self, event: &GenericEvent<T>) -> bool;
}

#[derive(Clone, Default)]
pub struct NoEventSelector;

impl EventSelector for NoEventSelector {
    fn process_event<T: FloatLike>(&mut self, _event: &GenericEvent<T>) -> bool {
        true
    }
}

#[derive(Debug, Clone)]
pub struct ConfiguredSelector {
    definition: ObservableDefinition,
    criteria: Vec<SelectorCriterion>,
}

impl ConfiguredSelector {
    pub(crate) fn new(
        settings: &SelectorSettings,
        quantities: &QuantitiesSettings,
    ) -> Result<Self> {
        let quantity = quantities.get(&settings.quantity).ok_or_else(|| {
            eyre!(
                "Selector references unknown quantity '{}'",
                settings.quantity
            )
        })?;

        let criterion = match &settings.selector {
            SelectorDefinitionSettings::ValueRange(selector) => SelectorCriterion::ValueRange {
                selection: settings.entry_selection,
                entry_index: settings.entry_index,
                reduction: selector.reduction,
                min: selector.min,
                max: selector.max,
            },
            SelectorDefinitionSettings::CountRange(selector) => SelectorCriterion::CountRange {
                selection: settings.entry_selection,
                entry_index: settings.entry_index,
                min: selector.min_count,
                max: selector.max_count,
            },
        };

        Ok(Self {
            definition: ObservableDefinition::from_settings(quantity),
            criteria: vec![criterion],
        })
    }
}

impl EventSelector for ConfiguredSelector {
    fn process_event<T: FloatLike>(&mut self, event: &GenericEvent<T>) -> bool {
        let entries = self.definition.process_event(event);
        self.criteria
            .iter()
            .all(|criterion| criterion.passes(&entries))
    }
}

#[derive(Clone)]
pub enum Selectors {
    All(NoEventSelector),
    Configured(ConfiguredSelector),
}

impl Selectors {
    pub(crate) fn from_settings(
        settings: &SelectorSettings,
        quantities: &QuantitiesSettings,
    ) -> Result<Self> {
        Ok(Selectors::Configured(ConfiguredSelector::new(
            settings, quantities,
        )?))
    }

    pub(crate) fn process_event<T: FloatLike>(&mut self, event: &GenericEvent<T>) -> bool {
        match self {
            Selectors::All(selector) => selector.process_event(event),
            Selectors::Configured(selector) => selector.process_event(event),
        }
    }
}

#[derive(Clone, Default)]
pub struct EventProcessingRuntime {
    selectors: Vec<Selectors>,
    observables: BTreeMap<String, Observables>,
}

impl EventProcessingRuntime {
    pub fn from_settings(settings: &RuntimeSettings) -> Result<Self> {
        let selectors = settings
            .selectors
            .values()
            .map(|selector| Selectors::from_settings(selector, &settings.quantities))
            .collect::<Result<Vec<_>>>()?;

        let observables = settings
            .observables
            .iter()
            .map(|(name, observable)| {
                Observables::from_settings(name, observable, &settings.quantities)
                    .map(|observable| (name.clone(), observable))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;

        Ok(Self {
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

    pub fn process_event<T: FloatLike>(&mut self, event: &GenericEvent<T>) -> bool {
        for selector in self.selectors.iter_mut() {
            if !selector.process_event(event) {
                return false;
            }
        }
        true
    }

    pub fn process_event_groups<T: FloatLike>(&mut self, event_groups: &GenericEventGroupList<T>) {
        for observable in self.observables.values_mut() {
            observable.process_event_groups(event_groups);
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

pub trait Observable {
    fn process_event_groups<T: FloatLike>(&mut self, event_groups: &GenericEventGroupList<T>);

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
        quantity_settings: &QuantitySettings,
        observable_name: &str,
        settings: &ObservableSettings,
    ) -> Self {
        let mut n_bins = settings.histogram.n_bins;
        let mut x_min = settings.histogram.x_min;
        let mut x_max = settings.histogram.x_max;
        if matches!(quantity_settings, QuantitySettings::CrossSection {}) {
            n_bins = 1;
            if x_max <= x_min {
                x_min = 0.0;
                x_max = 1.0;
            }
        }

        Self {
            definition,
            entry_selection: settings.entry_selection,
            entry_index: settings.entry_index,
            misbinning_max_normalized_distance: settings.misbinning_max_normalized_distance,
            state: HistogramAccumulatorState::new(
                observable_name.to_string(),
                settings.phase,
                settings.value_transform,
                x_min,
                x_max,
                settings.histogram.log_x_axis,
                settings.histogram.log_y_axis,
                n_bins,
            ),
            pending_contributions: Vec::with_capacity(16),
            candidate_pairs: Vec::with_capacity(16),
            used_contributions: Vec::with_capacity(16),
            grouped_contributions: vec![GroupBinContribution::default(); n_bins],
            grouped_underflow: GroupBinContribution::default(),
            grouped_overflow: GroupBinContribution::default(),
            touched_bins: Vec::with_capacity(n_bins.min(16)),
        }
    }

    fn build_group_contributions<T: FloatLike>(&mut self, event_group: &GenericEventGroup<T>) {
        self.pending_contributions.clear();
        self.pending_contributions
            .reserve(event_group.len().saturating_mul(4));

        for event in event_group.iter() {
            let entries = self.definition.process_event(event);
            for entry in apply_entry_selection(&entries, self.entry_selection, self.entry_index) {
                let transformed_value = transform_value(&entry.value, self.state.value_transform);
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

                let projected_weight = self
                    .state
                    .phase
                    .project(&combined_entry_weight(event, &entry))
                    .into_ff64()
                    .0;

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
}

impl Observable for HistogramObservable {
    fn process_event_groups<T: FloatLike>(&mut self, event_groups: &GenericEventGroupList<T>) {
        self.state.new_sample_count += 1;
        self.touched_bins.clear();

        for event_group in event_groups.iter() {
            self.build_group_contributions(event_group);
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
    pub(crate) fn from_settings(
        name: &str,
        settings: &ObservableSettings,
        quantities: &QuantitiesSettings,
    ) -> Result<Self> {
        let quantity = quantities.get(&settings.quantity).ok_or_else(|| {
            eyre!(
                "Observable '{}' references unknown quantity '{}'",
                name,
                settings.quantity
            )
        })?;

        Ok(Observables::Histogram(HistogramObservable::new(
            ObservableDefinition::from_settings(quantity),
            quantity,
            name,
            settings,
        )))
    }

    pub(crate) fn process_event_groups<T: FloatLike>(
        &mut self,
        event_groups: &GenericEventGroupList<T>,
    ) {
        match self {
            Observables::Histogram(observable) => observable.process_event_groups(event_groups),
        }
    }

    pub(crate) fn merge_samples(&mut self, other: &mut Observables) -> Result<()> {
        match (self, other) {
            (Observables::Histogram(lhs), Observables::Histogram(rhs)) => lhs.merge_samples(rhs),
        }
    }

    pub(crate) fn update_result(&mut self, iter: usize) {
        match self {
            Observables::Histogram(observable) => observable.update_result(iter),
        }
    }

    pub(crate) fn cleared_clone(&self) -> Self {
        match self {
            Observables::Histogram(observable) => Observables::Histogram(HistogramObservable {
                definition: observable.definition.clone(),
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

fn value_in_range(value: f64, min: f64, max: Option<f64>) -> bool {
    if value < min {
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
        HistogramBinSnapshot, HistogramSnapshot, HistogramStatisticsSnapshot, ObservablePhase,
        ObservableValueTransform,
    };
    use std::fs;

    #[test]
    fn histogram_snapshot_json_round_trip() {
        let snapshot = HistogramSnapshot {
            title: "top_pt".to_string(),
            phase: ObservablePhase::Real,
            value_transform: ObservableValueTransform::Identity,
            x_min: 0.0,
            x_max: 10.0,
            log_x_axis: false,
            log_y_axis: true,
            bins: vec![HistogramBinSnapshot {
                x_min: 0.0,
                x_max: 10.0,
                average: 1.25,
                error: 0.5,
                mitigated_fill_count: 2,
            }],
            statistics: HistogramStatisticsSnapshot {
                underflow_count: 1,
                overflow_count: 2,
                nan_value_count: 3,
                mitigated_pair_count: 4,
            },
        };

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
}
