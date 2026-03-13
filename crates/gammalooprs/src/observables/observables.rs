use super::clustering::{JetAlgorithm, JetClustering};
use super::events::{GenericEvent, GenericEventGroup};
use crate::settings::RuntimeSettings;
use crate::utils::serde_utils::{
    IsDefault, is_false, is_float, is_true, is_usize, show_defaults_helper,
};
use crate::utils::{F, FloatLike};
use bincode_trait_derive::{Decode, Encode};
use eyre::{Result, eyre};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smallvec::{SmallVec, smallvec};
use spenso::algebra::complex::Complex;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use symbolica::numerical_integration::StatisticsAccumulator;
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
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
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
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
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
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[serde(deny_unknown_fields)]
pub struct ParticleScalarQuantitySettings {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub pdgs: Vec<isize>,
    pub quantity: FilterQuantity,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[serde(deny_unknown_fields)]
pub struct JetPtQuantitySettings {
    #[serde(flatten)]
    pub clustering: JetClusteringSettings,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[serde(tag = "type", rename_all = "snake_case")]
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
pub enum QuantitySettings {
    ParticleScalar(ParticleScalarQuantitySettings),
    JetPt(JetPtQuantitySettings),
    AFB {},
    CrossSection {},
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
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
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[serde(deny_unknown_fields)]
pub struct CountRangeSelectorSettings {
    pub min_count: usize,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub max_count: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[serde(tag = "selector", rename_all = "snake_case")]
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
pub enum SelectorDefinitionSettings {
    ValueRange(ValueRangeSelectorSettings),
    CountRange(CountRangeSelectorSettings),
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[serde(deny_unknown_fields)]
pub struct SelectorSettings {
    pub quantity: String,
    #[serde(default, skip_serializing_if = "is_default_entry_selection")]
    pub entry_selection: EntrySelection,
    #[serde(default, skip_serializing_if = "is_default_entry_index")]
    pub entry_index: usize,
    #[serde(flatten)]
    pub selector: SelectorDefinitionSettings,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[allow(non_snake_case)]
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
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
#[cfg_attr(feature = "python_api", pyo3::pyclass)]
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
        let reference = event.integrand_weight.re.clone();
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

#[derive(Debug, Clone, Default)]
pub struct ObservableBinAccumulator {
    pub stats: StatisticsAccumulator<F<f64>>,
    pub mitigated_fill_count: usize,
    new_mitigated_fill_count: usize,
}

impl ObservableBinAccumulator {
    fn add_sample(&mut self, sample: F<f64>, mitigated_fill_count: usize) {
        self.stats.add_sample(sample, None);
        self.new_mitigated_fill_count += mitigated_fill_count;
    }

    fn merge_samples(&mut self, other: &mut ObservableBinAccumulator) {
        self.stats.merge_samples(&mut other.stats);
        self.new_mitigated_fill_count += other.new_mitigated_fill_count;
        other.new_mitigated_fill_count = 0;
    }

    fn update_iter(&mut self) {
        self.stats.update_iter(false);
        self.mitigated_fill_count += self.new_mitigated_fill_count;
        self.new_mitigated_fill_count = 0;
    }
}

#[derive(Debug, Clone, Default)]
pub struct ObservableHistogramStatistics {
    pub underflow_count: usize,
    pub overflow_count: usize,
    pub nan_value_count: usize,
    pub mitigated_pair_count: usize,
    new_underflow_count: usize,
    new_overflow_count: usize,
    new_nan_value_count: usize,
    new_mitigated_pair_count: usize,
}

impl ObservableHistogramStatistics {
    fn register_underflow(&mut self) {
        self.new_underflow_count += 1;
    }

    fn register_overflow(&mut self) {
        self.new_overflow_count += 1;
    }

    fn register_nan(&mut self) {
        self.new_nan_value_count += 1;
    }

    fn register_mitigated_pair(&mut self) {
        self.new_mitigated_pair_count += 1;
    }

    fn merge_samples(&mut self, other: &mut ObservableHistogramStatistics) {
        self.new_underflow_count += other.new_underflow_count;
        self.new_overflow_count += other.new_overflow_count;
        self.new_nan_value_count += other.new_nan_value_count;
        self.new_mitigated_pair_count += other.new_mitigated_pair_count;

        other.new_underflow_count = 0;
        other.new_overflow_count = 0;
        other.new_nan_value_count = 0;
        other.new_mitigated_pair_count = 0;
    }

    fn update_iter(&mut self) {
        self.underflow_count += self.new_underflow_count;
        self.overflow_count += self.new_overflow_count;
        self.nan_value_count += self.new_nan_value_count;
        self.mitigated_pair_count += self.new_mitigated_pair_count;

        self.new_underflow_count = 0;
        self.new_overflow_count = 0;
        self.new_nan_value_count = 0;
        self.new_mitigated_pair_count = 0;
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, JsonSchema)]
pub struct HistogramBinSnapshot {
    pub x_min: f64,
    pub x_max: f64,
    pub average: f64,
    pub error: f64,
    pub mitigated_fill_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, JsonSchema)]
pub struct HistogramStatisticsSnapshot {
    pub underflow_count: usize,
    pub overflow_count: usize,
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
    pub log_x_axis: bool,
    pub log_y_axis: bool,
    pub bins: Vec<HistogramBinSnapshot>,
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
}

impl EventProcessingRuntime {
    pub(crate) fn from_settings(settings: &RuntimeSettings) -> Result<Self> {
        let selectors = settings
            .selectors
            .values()
            .map(|selector| Selectors::from_settings(selector, &settings.quantities))
            .collect::<Result<Vec<_>>>()?;

        Ok(Self { selectors })
    }

    pub(crate) fn has_selectors(&self) -> bool {
        !self.selectors.is_empty()
    }

    pub(crate) fn process_event<T: FloatLike>(&mut self, event: &GenericEvent<T>) -> bool {
        for selector in self.selectors.iter_mut() {
            if !selector.process_event(event) {
                return false;
            }
        }
        true
    }
}

pub trait Observable {
    fn process_event<T: FloatLike>(
        &mut self,
        event_group: &GenericEventGroup<T>,
        integrator_weight: F<T>,
    );

    fn merge_samples(&mut self, other: &mut Self)
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
    mitigated_fill_count: usize,
}

#[derive(Debug, Clone)]
pub struct HistogramObservable {
    definition: ObservableDefinition,
    entry_selection: EntrySelection,
    entry_index: usize,
    value_transform: ObservableValueTransform,
    phase: ObservablePhase,
    misbinning_max_normalized_distance: Option<f64>,
    x_min: f64,
    x_max: f64,
    title: String,
    bins: Vec<ObservableBinAccumulator>,
    histogram_stats: ObservableHistogramStatistics,
    log_x_axis: bool,
    log_y_axis: bool,
    pending_contributions: Vec<PendingBinContribution>,
    candidate_pairs: Vec<(f64, usize, usize)>,
    used_contributions: Vec<bool>,
    grouped_contributions: Vec<GroupBinContribution>,
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
            value_transform: settings.value_transform,
            phase: settings.phase,
            misbinning_max_normalized_distance: settings.misbinning_max_normalized_distance,
            x_min,
            x_max,
            title: observable_name.to_string(),
            bins: vec![ObservableBinAccumulator::default(); n_bins],
            histogram_stats: ObservableHistogramStatistics::default(),
            log_x_axis: settings.histogram.log_x_axis,
            log_y_axis: settings.histogram.log_y_axis,
            pending_contributions: Vec::with_capacity(16),
            candidate_pairs: Vec::with_capacity(16),
            used_contributions: Vec::with_capacity(16),
            grouped_contributions: vec![GroupBinContribution::default(); n_bins],
            touched_bins: Vec::with_capacity(n_bins.min(16)),
        }
    }

    fn build_group_contributions<T: FloatLike>(
        &mut self,
        event_group: &GenericEventGroup<T>,
        integrator_weight: F<T>,
    ) {
        self.pending_contributions.clear();
        self.pending_contributions
            .reserve(event_group.len().saturating_mul(4));

        for event in event_group.iter() {
            let entries = self.definition.process_event(event);
            for entry in apply_entry_selection(&entries, self.entry_selection, self.entry_index) {
                let transformed_value = transform_value(&entry.value, self.value_transform);
                let value = transformed_value.into_ff64().0;
                if !value.is_finite() {
                    self.histogram_stats.register_nan();
                    continue;
                }

                let Some(bin_position) =
                    histogram_bin_position(value, self.x_min, self.x_max, self.bins.len())
                else {
                    self.histogram_stats.register_nan();
                    continue;
                };

                let projected_weight = self
                    .phase
                    .project(&combined_entry_weight(event, &entry, &integrator_weight))
                    .into_ff64()
                    .0;

                match bin_position {
                    HistogramBinPosition::Underflow => self.histogram_stats.register_underflow(),
                    HistogramBinPosition::Overflow => self.histogram_stats.register_overflow(),
                    HistogramBinPosition::InRange(bin_index) => {
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

        if self.pending_contributions.len() < 2 || self.bins.len() < 2 || max_distance <= 0.0 {
            return;
        }

        let bin_width = (self.x_max - self.x_min) / self.bins.len() as f64;
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
            self.histogram_stats.register_mitigated_pair();
        }
    }

    pub fn snapshot(&self) -> HistogramSnapshot {
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
                    x_min,
                    x_max,
                    average: bin.stats.avg.0,
                    error: bin.stats.err.0,
                    mitigated_fill_count: bin.mitigated_fill_count,
                }
            })
            .collect();

        HistogramSnapshot {
            title: self.title.clone(),
            phase: self.phase,
            value_transform: self.value_transform,
            x_min: self.x_min,
            x_max: self.x_max,
            log_x_axis: self.log_x_axis,
            log_y_axis: self.log_y_axis,
            bins,
            statistics: HistogramStatisticsSnapshot {
                underflow_count: self.histogram_stats.underflow_count,
                overflow_count: self.histogram_stats.overflow_count,
                nan_value_count: self.histogram_stats.nan_value_count,
                mitigated_pair_count: self.histogram_stats.mitigated_pair_count,
            },
        }
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
                c1, c2, bin.stats.avg.0, bin.stats.err.0,
            )?;
        }

        writeln!(writer, "<\\histogram>")?;
        Ok(())
    }
}

impl Observable for HistogramObservable {
    fn process_event<T: FloatLike>(
        &mut self,
        event_group: &GenericEventGroup<T>,
        integrator_weight: F<T>,
    ) {
        self.build_group_contributions(event_group, integrator_weight);
        self.mitigate_group_misbinning();

        self.touched_bins.clear();
        for contribution in &self.pending_contributions {
            let grouped = &mut self.grouped_contributions[contribution.bin_index];
            if grouped.projected_weight == 0.0 && grouped.mitigated_fill_count == 0 {
                self.touched_bins.push(contribution.bin_index);
            }
            grouped.projected_weight += contribution.projected_weight;
            if contribution.mitigated {
                grouped.mitigated_fill_count += 1;
            }
        }

        for (bin_index, bin_accumulator) in self.bins.iter_mut().enumerate() {
            let grouped = &self.grouped_contributions[bin_index];
            if grouped.projected_weight != 0.0 || grouped.mitigated_fill_count != 0 {
                bin_accumulator
                    .add_sample(F(grouped.projected_weight), grouped.mitigated_fill_count);
            } else {
                bin_accumulator.add_sample(F(0.0), 0);
            }
        }

        for bin_index in self.touched_bins.drain(..) {
            self.grouped_contributions[bin_index] = GroupBinContribution::default();
        }
    }

    fn merge_samples(&mut self, other: &mut HistogramObservable) {
        for (bin, other_bin) in self.bins.iter_mut().zip(other.bins.iter_mut()) {
            bin.merge_samples(other_bin);
        }
        self.histogram_stats
            .merge_samples(&mut other.histogram_stats);
    }

    fn update_result(&mut self, _iter: usize) {
        for bin in &mut self.bins {
            bin.update_iter();
        }
        self.histogram_stats.update_iter();

        for (i, bin) in self.bins.iter().enumerate() {
            let c1 = (self.x_max - self.x_min) * i as f64 / self.bins.len() as f64 + self.x_min;
            let c2 = (self.x_max - self.x_min) * (i + 1) as f64 / self.bins.len() as f64;
            info!("{}={}: {} +/- {}", c1, c2, bin.stats.avg, bin.stats.err);
        }

        info!(
            "{} stats: underflow={}, overflow={}, nan_values={}, mitigated_pairs={}",
            self.title,
            self.histogram_stats.underflow_count,
            self.histogram_stats.overflow_count,
            self.histogram_stats.nan_value_count,
            self.histogram_stats.mitigated_pair_count,
        );
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

    pub(crate) fn process_event<T: FloatLike>(
        &mut self,
        event_group: &GenericEventGroup<T>,
        integrator_weight: F<T>,
    ) {
        match self {
            Observables::Histogram(observable) => {
                observable.process_event(event_group, integrator_weight)
            }
        }
    }

    pub(crate) fn merge_samples(&mut self, other: &mut Observables) {
        match (self, other) {
            (Observables::Histogram(lhs), Observables::Histogram(rhs)) => lhs.merge_samples(rhs),
        }
    }

    pub(crate) fn update_result(&mut self, iter: usize) {
        match self {
            Observables::Histogram(observable) => observable.update_result(iter),
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
    integrator_weight: &F<T>,
) -> Complex<F<T>> {
    let combined = event.integrand_weight.clone()
        * event.observable_weight.clone()
        * entry.weight_modifier.clone();
    Complex::new(
        combined.re * integrator_weight.clone(),
        combined.im * integrator_weight.clone(),
    )
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
