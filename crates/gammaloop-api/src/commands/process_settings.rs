use std::{collections::BTreeMap, sync::OnceLock};

use color_eyre::Result;
use gammalooprs::{
    observables::{
        CountRangeSelectorSettings, EntrySelection, FilterQuantity, HistogramSettings,
        JetClusteringSettings, JetQuantitySettings, ObservablePhase, ObservableSettings,
        ObservableValueTransform, PairQuantity, ParticleQuantitySettings, QuantityComputation,
        QuantityComputationSettings, QuantityOrdering, QuantitySettings,
        SelectorDefinitionSettings, SelectorReduction, SelectorSettings,
        ValueRangeSelectorSettings,
    },
    settings::RuntimeSettings,
};
use serde::Serialize;
use serde_json::Value as JsonValue;

use crate::settings_tree::{serialize_schema, serialize_settings_with_defaults};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NamedProcessSettingKind {
    Quantity,
    Observable,
    Selector,
}

impl NamedProcessSettingKind {
    pub(crate) fn singular(self) -> &'static str {
        match self {
            NamedProcessSettingKind::Quantity => "quantity",
            NamedProcessSettingKind::Observable => "observable",
            NamedProcessSettingKind::Selector => "selector",
        }
    }

    pub(crate) fn plural(self) -> &'static str {
        match self {
            NamedProcessSettingKind::Quantity => "quantities",
            NamedProcessSettingKind::Observable => "observables",
            NamedProcessSettingKind::Selector => "selectors",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ProcessSettingsCompletionEntry {
    pub process_id: usize,
    pub process_name: String,
    pub integrand_name: String,
    pub quantities: BTreeMap<String, JsonValue>,
    pub observables: BTreeMap<String, JsonValue>,
    pub selectors: BTreeMap<String, JsonValue>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SerializedRuntimeNamedSettings {
    pub quantities: BTreeMap<String, JsonValue>,
    pub observables: BTreeMap<String, JsonValue>,
    pub selectors: BTreeMap<String, JsonValue>,
}

#[derive(Debug, Clone)]
struct KindedTemplate<T> {
    kind: String,
    settings: T,
    root: JsonValue,
}

fn build_kinded_template<T>(settings: T, tag_key: &str, context: &str) -> KindedTemplate<T>
where
    T: Clone + Serialize,
{
    let root = serialize_settings_with_defaults(&settings, context)
        .expect("template settings must serialize with defaults");
    let kind = root
        .get(tag_key)
        .and_then(JsonValue::as_str)
        .unwrap_or_else(|| panic!("missing tag '{tag_key}' in {context}"))
        .to_string();
    KindedTemplate {
        kind,
        settings,
        root,
    }
}

fn quantity_templates() -> &'static [KindedTemplate<QuantitySettings>] {
    static TEMPLATES: OnceLock<Vec<KindedTemplate<QuantitySettings>>> = OnceLock::new();
    TEMPLATES.get_or_init(|| {
        vec![
            build_kinded_template(
                QuantitySettings::Particle(ParticleQuantitySettings {
                    pdgs: Vec::new(),
                    computation: QuantityComputationSettings::scalar(FilterQuantity::PT),
                }),
                "type",
                "particle quantity template",
            ),
            build_kinded_template(
                QuantitySettings::Jet(JetQuantitySettings {
                    clustering: JetClusteringSettings::default(),
                    computation: QuantityComputationSettings::scalar(FilterQuantity::PT),
                }),
                "type",
                "jet quantity template",
            ),
            build_kinded_template(QuantitySettings::AFB {}, "type", "afb quantity template"),
            build_kinded_template(
                QuantitySettings::CrossSection {},
                "type",
                "cross-section quantity template",
            ),
        ]
    })
}

fn selector_templates() -> &'static [KindedTemplate<SelectorSettings>] {
    static TEMPLATES: OnceLock<Vec<KindedTemplate<SelectorSettings>>> = OnceLock::new();
    TEMPLATES.get_or_init(|| {
        vec![
            build_kinded_template(
                SelectorSettings {
                    quantity: String::new(),
                    entry_selection: EntrySelection::All,
                    entry_index: 0,
                    selector: SelectorDefinitionSettings::ValueRange(ValueRangeSelectorSettings {
                        min: None,
                        max: None,
                        reduction: SelectorReduction::AnyInRange,
                    }),
                },
                "selector",
                "value-range selector template",
            ),
            build_kinded_template(
                SelectorSettings {
                    quantity: String::new(),
                    entry_selection: EntrySelection::All,
                    entry_index: 0,
                    selector: SelectorDefinitionSettings::CountRange(CountRangeSelectorSettings {
                        min_count: 0,
                        max_count: None,
                    }),
                },
                "selector",
                "count-range selector template",
            ),
        ]
    })
}

pub(crate) fn observable_template() -> ObservableSettings {
    ObservableSettings {
        quantity: String::new(),
        entry_selection: EntrySelection::All,
        entry_index: 0,
        value_transform: ObservableValueTransform::Identity,
        phase: ObservablePhase::Real,
        misbinning_max_normalized_distance: None,
        histogram: HistogramSettings::default(),
    }
}

pub(crate) fn quantity_kind_names() -> Vec<String> {
    quantity_templates()
        .iter()
        .map(|entry| entry.kind.clone())
        .collect()
}

pub(crate) fn selector_kind_names() -> Vec<String> {
    selector_templates()
        .iter()
        .map(|entry| entry.kind.clone())
        .collect()
}

pub(crate) fn parse_quantity_kind(raw: &str) -> std::result::Result<String, String> {
    if quantity_templates().iter().any(|entry| entry.kind == raw) {
        Ok(raw.to_string())
    } else {
        Err(format!(
            "Unknown quantity kind '{raw}', expected one of: {}",
            quantity_kind_names().join(", ")
        ))
    }
}

pub(crate) fn parse_selector_kind(raw: &str) -> std::result::Result<String, String> {
    if selector_templates().iter().any(|entry| entry.kind == raw) {
        Ok(raw.to_string())
    } else {
        Err(format!(
            "Unknown selector kind '{raw}', expected one of: {}",
            selector_kind_names().join(", ")
        ))
    }
}

pub(crate) fn quantity_template(kind: &str) -> Option<QuantitySettings> {
    quantity_templates()
        .iter()
        .find(|entry| entry.kind == kind)
        .map(|entry| entry.settings.clone())
}

pub(crate) fn selector_template(kind: &str) -> Option<SelectorSettings> {
    selector_templates()
        .iter()
        .find(|entry| entry.kind == kind)
        .map(|entry| entry.settings.clone())
}

pub(crate) fn quantity_completion_root_for_kind(kind: &str) -> Option<&'static JsonValue> {
    quantity_templates()
        .iter()
        .find(|entry| entry.kind == kind)
        .map(|entry| &entry.root)
}

pub(crate) fn selector_completion_root_for_kind(kind: &str) -> Option<&'static JsonValue> {
    selector_templates()
        .iter()
        .find(|entry| entry.kind == kind)
        .map(|entry| &entry.root)
}

pub(crate) fn observable_completion_root() -> &'static JsonValue {
    static ROOT: OnceLock<JsonValue> = OnceLock::new();
    ROOT.get_or_init(|| {
        serialize_settings_with_defaults(&observable_template(), "observable settings template")
            .expect("observable settings template must serialize with defaults")
    })
}

pub(crate) fn quantity_schema() -> &'static JsonValue {
    static ROOT: OnceLock<JsonValue> = OnceLock::new();
    ROOT.get_or_init(|| {
        serialize_schema::<QuantitySettings>("quantity settings completion")
            .expect("quantity settings schema must serialize")
    })
}

pub(crate) fn selector_schema() -> &'static JsonValue {
    static ROOT: OnceLock<JsonValue> = OnceLock::new();
    ROOT.get_or_init(|| {
        serialize_schema::<SelectorSettings>("selector settings completion")
            .expect("selector settings schema must serialize")
    })
}

pub(crate) fn observable_schema() -> &'static JsonValue {
    static ROOT: OnceLock<JsonValue> = OnceLock::new();
    ROOT.get_or_init(|| {
        serialize_schema::<ObservableSettings>("observable settings completion")
            .expect("observable settings schema must serialize")
    })
}

pub(crate) fn serialize_runtime_named_settings(
    settings: &RuntimeSettings,
) -> Result<SerializedRuntimeNamedSettings> {
    Ok(SerializedRuntimeNamedSettings {
        quantities: serialize_named_settings_map(&settings.quantities, "quantity settings")?,
        observables: serialize_named_settings_map(&settings.observables, "observable settings")?,
        selectors: serialize_named_settings_map(&settings.selectors, "selector settings")?,
    })
}

fn serialize_named_settings_map<T>(
    map: &BTreeMap<String, T>,
    context: &str,
) -> Result<BTreeMap<String, JsonValue>>
where
    T: Serialize,
{
    map.iter()
        .map(|(name, settings)| {
            serialize_settings_with_defaults(settings, context).map(|root| (name.clone(), root))
        })
        .collect()
}

pub(crate) fn quantity_kind(settings: &QuantitySettings) -> &'static str {
    match settings {
        QuantitySettings::Particle(_) => "particle",
        QuantitySettings::Jet(_) => "jet",
        QuantitySettings::AFB {} => "afb",
        QuantitySettings::CrossSection {} => "cross_section",
    }
}

pub(crate) fn selector_kind(settings: &SelectorSettings) -> &'static str {
    match &settings.selector {
        SelectorDefinitionSettings::ValueRange(_) => "value_range",
        SelectorDefinitionSettings::CountRange(_) => "count_range",
    }
}

pub(crate) fn observable_kind(_settings: &ObservableSettings) -> &'static str {
    "histogram"
}

pub(crate) fn summarize_quantity(settings: &QuantitySettings) -> String {
    match settings {
        QuantitySettings::Particle(settings) => summarize_particle_quantity(settings),
        QuantitySettings::Jet(settings) => summarize_jet_quantity(settings),
        QuantitySettings::AFB {} => "forward/backward asymmetry".to_string(),
        QuantitySettings::CrossSection {} => "cross section weight".to_string(),
    }
}

fn summarize_particle_quantity(settings: &ParticleQuantitySettings) -> String {
    let QuantitySettings::Particle(settings) =
        QuantitySettings::Particle(settings.clone()).normalized()
    else {
        unreachable!("particle quantity normalization must preserve its variant");
    };
    format!(
        "{} pdgs={}",
        summarize_quantity_computation(&settings.computation),
        format_integer_list(&settings.pdgs)
    )
}

fn summarize_jet_quantity(settings: &JetQuantitySettings) -> String {
    let QuantitySettings::Jet(settings) = QuantitySettings::Jet(settings.clone()).normalized()
    else {
        unreachable!("jet quantity normalization must preserve its variant");
    };
    format!(
        "{} algorithm={:?} dR={} min_jpt={} clustered_pdgs={}",
        summarize_quantity_computation(&settings.computation),
        settings.clustering.algorithm,
        settings.clustering.dR,
        settings.clustering.min_jpt,
        settings
            .clustering
            .clustered_pdgs
            .as_ref()
            .map(|pdgs| format!("{pdgs:?}"))
            .unwrap_or_else(|| "None (model default)".to_string())
    )
}

fn summarize_quantity_computation(settings: &QuantityComputationSettings) -> String {
    match settings.computation {
        QuantityComputation::Scalar => format!(
            "computation=scalar quantity={} ordering={} order={}",
            settings.quantity.unwrap_or(FilterQuantity::PT),
            settings.ordering.unwrap_or(QuantityOrdering::Quantity),
            settings.order
        ),
        QuantityComputation::Count => "computation=count".to_string(),
        QuantityComputation::Pair => format!(
            "computation=pair pair_quantity={} pairing={} ordering={} order={}",
            settings.pair_quantity.unwrap_or(PairQuantity::DeltaR),
            settings.pairing.unwrap_or_default(),
            settings.ordering.unwrap_or(QuantityOrdering::Quantity),
            settings.order
        ),
    }
}

pub(crate) fn summarize_observable(settings: &ObservableSettings) -> String {
    format!(
        "quantity={} selection={} transform={:?} phase={:?} bins={} range=[{}, {}] title={} type_description={}",
        format_reference(&settings.quantity),
        format_entry_selection(settings.entry_selection, settings.entry_index),
        settings.value_transform,
        settings.phase,
        settings.histogram.n_bins,
        settings.histogram.x_min,
        settings.histogram.x_max,
        settings
            .histogram
            .title
            .clone()
            .unwrap_or_else(|| "<observable name>".to_string()),
        settings.histogram.type_description,
    )
}

pub(crate) fn summarize_selector(settings: &SelectorSettings) -> String {
    let selector_summary = match &settings.selector {
        SelectorDefinitionSettings::ValueRange(selector) => format!(
            "range=[{}, {}] reduction={:?}",
            selector
                .min
                .map(|value| value.to_string())
                .unwrap_or_else(|| "-inf".to_string()),
            selector
                .max
                .map(|value| value.to_string())
                .unwrap_or_else(|| "inf".to_string()),
            selector.reduction
        ),
        SelectorDefinitionSettings::CountRange(selector) => format!(
            "count=[{}, {}]",
            selector.min_count,
            selector
                .max_count
                .map(|value| value.to_string())
                .unwrap_or_else(|| "inf".to_string())
        ),
    };

    format!(
        "quantity={} selection={} {selector_summary}",
        format_reference(&settings.quantity),
        format_entry_selection(settings.entry_selection, settings.entry_index),
    )
}

fn format_entry_selection(selection: EntrySelection, entry_index: usize) -> String {
    match selection {
        EntrySelection::All => "all".to_string(),
        EntrySelection::LeadingOnly => "leading_only".to_string(),
        EntrySelection::NthOnly => format!("nth_only[{entry_index}]"),
    }
}

fn format_reference(name: &str) -> String {
    if name.is_empty() {
        "<unset>".to_string()
    } else {
        name.to_string()
    }
}

fn format_integer_list(values: &[isize]) -> String {
    if values.is_empty() {
        "[]".to_string()
    } else {
        format!(
            "[{}]",
            values
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}
