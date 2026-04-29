use super::clustering::ClusteringResult;
use crate::momentum::{FourMomentum, Rotation};
use crate::utils::{F, FloatLike, into_complex_ff64};
use colored::Colorize;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use spenso::algebra::complex::Complex;
use std::collections::BTreeMap;
use std::fmt;
use std::ops::{Deref, DerefMut};
use tabled::{
    Table, Tabled,
    settings::{Style, style::HorizontalLine},
};

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct CutInfo {
    pub particle_pdgs: (SmallVec<[isize; 2]>, SmallVec<[isize; 4]>),
    pub cut_id: usize,
    pub graph_id: usize,
    pub graph_group_id: Option<usize>,
    pub orientation_id: Option<usize>,
    pub lmb_channel_id: Option<usize>,
    pub lmb_channel_edge_ids: Option<SmallVec<[usize; 4]>>,
}

pub type Event = GenericEvent<f64>;
pub type EventGroup = GenericEventGroup<f64>;
pub type EventGroupList = GenericEventGroupList<f64>;

#[derive(Default, Debug, Clone)]
pub struct GenericDerivedEventData<T: FloatLike> {
    pub clustered_jets: Vec<Option<ClusteringResult<T>>>,
}

impl<T: FloatLike> GenericDerivedEventData<T> {
    pub fn to_f64(&self) -> GenericDerivedEventData<f64> {
        GenericDerivedEventData {
            clustered_jets: self
                .clustered_jets
                .iter()
                .map(|result| result.as_ref().map(ClusteringResult::to_f64))
                .collect(),
        }
    }

    pub fn from_f64(data: &GenericDerivedEventData<f64>) -> Self {
        GenericDerivedEventData {
            clustered_jets: data
                .clustered_jets
                .iter()
                .map(|result| result.as_ref().map(ClusteringResult::from_f64))
                .collect(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum AdditionalWeightKey {
    FullMultiplicativeFactor,
    Original,
    ThresholdCounterterm { subset_index: usize },
    AmplitudeThresholdCounterterm {
        esurface_id: usize,
        overlap_group: usize,
    },
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct GenericEventGroup<T: FloatLike>(pub Vec<GenericEvent<T>>);

impl<T: FloatLike> Deref for GenericEventGroup<T> {
    type Target = Vec<GenericEvent<T>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: FloatLike> DerefMut for GenericEventGroup<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T: FloatLike> GenericEventGroup<T> {
    pub fn to_f64(&self) -> EventGroup {
        GenericEventGroup(self.0.iter().map(GenericEvent::to_f64).collect())
    }

    pub fn from_f64(event_group: &EventGroup) -> Self {
        GenericEventGroup(event_group.0.iter().map(GenericEvent::from_f64).collect())
    }
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct GenericEventGroupList<T: FloatLike>(pub Vec<GenericEventGroup<T>>);

impl<T: FloatLike> Deref for GenericEventGroupList<T> {
    type Target = Vec<GenericEventGroup<T>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: FloatLike> DerefMut for GenericEventGroupList<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T: FloatLike> GenericEventGroupList<T> {
    pub fn to_f64(&self) -> EventGroupList {
        GenericEventGroupList(self.0.iter().map(GenericEventGroup::to_f64).collect())
    }

    pub fn from_f64(event_group_list: &EventGroupList) -> Self {
        GenericEventGroupList(
            event_group_list
                .0
                .iter()
                .map(GenericEventGroup::from_f64)
                .collect(),
        )
    }

    pub fn push_singleton(&mut self, event: GenericEvent<T>) {
        self.0.push(GenericEventGroup(vec![event]));
    }
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct GenericAdditionalWeightInfo<T: FloatLike> {
    pub weights: BTreeMap<AdditionalWeightKey, Complex<F<T>>>,
}

impl<T: FloatLike> GenericAdditionalWeightInfo<T> {
    pub fn to_f64(&self) -> GenericAdditionalWeightInfo<f64> {
        GenericAdditionalWeightInfo {
            weights: self
                .weights
                .iter()
                .map(|(key, weight)| (*key, into_complex_ff64(weight)))
                .collect(),
        }
    }

    pub fn from_f64(additional_weights: &GenericAdditionalWeightInfo<f64>) -> Self {
        GenericAdditionalWeightInfo {
            weights: additional_weights
                .weights
                .iter()
                .map(|(key, weight)| {
                    (
                        *key,
                        Complex::new(F::from_ff64(weight.re), F::from_ff64(weight.im)),
                    )
                })
                .collect(),
        }
    }
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct GenericEvent<T: FloatLike> {
    #[allow(clippy::type_complexity)]
    pub kinematic_configuration: (
        SmallVec<[FourMomentum<F<T>>; 2]>,
        SmallVec<[FourMomentum<F<T>>; 4]>,
    ),
    pub cut_info: CutInfo,
    // Stores the fully normalized contribution associated with this event.
    pub weight: Complex<F<T>>,
    // Contains additional partial weights for future use to do post-processing reweighting.
    pub additional_weights: GenericAdditionalWeightInfo<T>,
    #[serde(skip)]
    pub derived_observable_data: GenericDerivedEventData<T>,
}

impl<T: FloatLike> GenericEvent<T> {
    pub(crate) fn inverse_rotate(&mut self, rotation: &Rotation) {
        if rotation.is_identity() {
            return;
        }

        for momentum in self.kinematic_configuration.0.iter_mut() {
            *momentum = rotation.inverse_rotate_four(momentum);
        }
        for momentum in self.kinematic_configuration.1.iter_mut() {
            *momentum = rotation.inverse_rotate_four(momentum);
        }

        // Any cached clustering belongs to the rotated frame and must be recomputed.
        self.derived_observable_data = GenericDerivedEventData::default();
    }

    pub fn ensure_clustering_slots(&mut self, n_clusterings: usize) {
        if self.derived_observable_data.clustered_jets.len() < n_clusterings {
            self.derived_observable_data
                .clustered_jets
                .resize_with(n_clusterings, || None);
        }
    }

    pub fn cached_clustering(&self, handle: usize) -> Option<&ClusteringResult<T>> {
        self.derived_observable_data
            .clustered_jets
            .get(handle)
            .and_then(|result| result.as_ref())
    }

    pub fn to_f64(&self) -> Event {
        GenericEvent {
            kinematic_configuration: (
                self.kinematic_configuration
                    .0
                    .iter()
                    .map(FourMomentum::to_f64)
                    .collect(),
                self.kinematic_configuration
                    .1
                    .iter()
                    .map(FourMomentum::to_f64)
                    .collect(),
            ),
            cut_info: self.cut_info.clone(),
            weight: into_complex_ff64(&self.weight),
            additional_weights: self.additional_weights.to_f64(),
            derived_observable_data: self.derived_observable_data.to_f64(),
        }
    }

    pub fn from_f64(event: &Event) -> Self {
        GenericEvent {
            kinematic_configuration: (
                event
                    .kinematic_configuration
                    .0
                    .iter()
                    .map(FourMomentum::from_ff64)
                    .collect(),
                event
                    .kinematic_configuration
                    .1
                    .iter()
                    .map(FourMomentum::from_ff64)
                    .collect(),
            ),
            cut_info: event.cut_info.clone(),
            weight: Complex::new(F::from_ff64(event.weight.re), F::from_ff64(event.weight.im)),
            additional_weights: GenericAdditionalWeightInfo::from_f64(&event.additional_weights),
            derived_observable_data: GenericDerivedEventData::from_f64(
                &event.derived_observable_data,
            ),
        }
    }
}

#[derive(Tabled)]
struct EventSummaryRow {
    field: String,
    value: String,
}

#[derive(Tabled)]
struct MomentumRow {
    #[tabled(rename = "state")]
    state: String,
    #[tabled(rename = "PDG")]
    pdg: String,
    #[tabled(rename = "E")]
    e: String,
    px: String,
    py: String,
    pz: String,
    #[tabled(rename = "sqrt(p^2)")]
    p2: String,
}

#[derive(Tabled)]
struct AdditionalWeightRow {
    key: String,
    value: String,
}

fn display_decimal_precision<T: FloatLike>(value: &F<T>) -> usize {
    let precision_bits = value.0.get_precision().max(1) as f64;
    (precision_bits * std::f64::consts::LOG10_2).ceil().max(1.0) as usize
}

pub(crate) fn format_real_generic<T: FloatLike>(value: &F<T>) -> String {
    format!("{:+.*e}", display_decimal_precision(value), value)
}

fn format_count(value: usize) -> String {
    if value < 1_000 {
        return value.to_string();
    }

    let value = value as f64;
    for (scale, suffix) in [
        (1_000_000_000_f64, "B"),
        (1_000_000_f64, "M"),
        (1_000_f64, "K"),
    ] {
        if value >= scale {
            let scaled = value / scale;
            let precision = if scaled >= 100.0 {
                0
            } else if scaled >= 10.0 {
                1
            } else {
                2
            };
            return format!("{scaled:.precision$}{suffix}");
        }
    }

    value.round().to_string()
}

pub(crate) fn format_optional_real_generic<T: FloatLike>(value: Option<&F<T>>) -> String {
    value
        .map(format_real_generic)
        .unwrap_or_else(|| "None".red().to_string())
}

pub(crate) fn format_complex_generic<T: FloatLike>(value: &Complex<F<T>>) -> String {
    let precision = display_decimal_precision(&value.re).max(display_decimal_precision(&value.im));
    format!("{:+.*e} {:+.*e}i", precision, value.re, precision, value.im)
}

fn format_pdg(pdg: Option<isize>) -> String {
    pdg.map(|value| value.to_string())
        .unwrap_or_else(|| "N/A".red().to_string())
}

fn format_pdg_with_state_color(pdg: Option<isize>, incoming: bool) -> String {
    let pdg = format_pdg(pdg);
    if incoming {
        pdg.bright_blue().to_string()
    } else {
        pdg.bright_green().to_string()
    }
}

fn format_lmb_channel_edge_ids(edge_ids: Option<&[usize]>) -> String {
    edge_ids
        .map(|edge_ids| {
            format!(
                "({})",
                edge_ids
                    .iter()
                    .map(|edge_id| edge_id.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .unwrap_or_else(|| "None".to_string())
}

fn momentum_mass_squared<T: FloatLike>(momentum: &FourMomentum<F<T>>) -> F<T> {
    momentum.temporal.value.clone() * momentum.temporal.value.clone()
        - momentum.spatial.px.clone() * momentum.spatial.px.clone()
        - momentum.spatial.py.clone() * momentum.spatial.py.clone()
        - momentum.spatial.pz.clone() * momentum.spatial.pz.clone()
}

fn momentum_mass<T: FloatLike>(momentum: &FourMomentum<F<T>>) -> F<T> {
    momentum_mass_squared(momentum).abs().sqrt()
}

fn format_momentum_rows<T: FloatLike>(
    state: String,
    momenta: &[FourMomentum<F<T>>],
    pdgs: &[isize],
    incoming: bool,
) -> Vec<MomentumRow> {
    momenta
        .iter()
        .enumerate()
        .map(|(index, momentum)| MomentumRow {
            state: state.clone(),
            pdg: format_pdg_with_state_color(pdgs.get(index).copied(), incoming),
            e: format_real_generic(&momentum.temporal.value),
            px: format_real_generic(&momentum.spatial.px),
            py: format_real_generic(&momentum.spatial.py),
            pz: format_real_generic(&momentum.spatial.pz),
            p2: format_real_generic(&momentum_mass(momentum))
                .bright_yellow()
                .to_string(),
        })
        .collect()
}

fn event_zero<T: FloatLike>(event: &GenericEvent<T>) -> F<T> {
    event
        .kinematic_configuration
        .0
        .first()
        .map(|momentum| momentum.temporal.value.zero())
        .or_else(|| {
            event
                .kinematic_configuration
                .1
                .first()
                .map(|momentum| momentum.temporal.value.zero())
        })
        .unwrap_or_else(|| F(T::new_zero()))
}

fn conservation_row<T: FloatLike>(event: &GenericEvent<T>) -> MomentumRow {
    let zero = event_zero(event);
    let incoming = event.kinematic_configuration.0.iter().fold(
        (zero.clone(), zero.clone(), zero.clone(), zero.clone()),
        |acc, momentum| {
            (
                acc.0 + momentum.temporal.value.clone(),
                acc.1 + momentum.spatial.px.clone(),
                acc.2 + momentum.spatial.py.clone(),
                acc.3 + momentum.spatial.pz.clone(),
            )
        },
    );
    let outgoing = event.kinematic_configuration.1.iter().fold(
        (zero.clone(), zero.clone(), zero.clone(), zero),
        |acc, momentum| {
            (
                acc.0 + momentum.temporal.value.clone(),
                acc.1 + momentum.spatial.px.clone(),
                acc.2 + momentum.spatial.py.clone(),
                acc.3 + momentum.spatial.pz.clone(),
            )
        },
    );
    let delta = FourMomentum {
        temporal: crate::momentum::Energy::new(incoming.0 - outgoing.0),
        spatial: crate::momentum::ThreeMomentum::new(
            incoming.1 - outgoing.1,
            incoming.2 - outgoing.2,
            incoming.3 - outgoing.3,
        ),
    };

    MomentumRow {
        state: "CHECK".bold().bright_yellow().to_string(),
        pdg: "N/A".bright_yellow().to_string(),
        e: format_real_generic(&delta.temporal.value)
            .bright_yellow()
            .to_string(),
        px: format_real_generic(&delta.spatial.px)
            .bright_yellow()
            .to_string(),
        py: format_real_generic(&delta.spatial.py)
            .bright_yellow()
            .to_string(),
        pz: format_real_generic(&delta.spatial.pz)
            .bright_yellow()
            .to_string(),
        p2: String::new().bright_yellow().to_string(),
    }
}

fn indent_block(block: &str, prefix: &str) -> String {
    let mut result = String::new();
    for (index, line) in block.lines().enumerate() {
        if index > 0 {
            result.push('\n');
        }
        result.push_str(prefix);
        result.push_str(line);
    }
    result
}

impl fmt::Display for AdditionalWeightKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AdditionalWeightKey::FullMultiplicativeFactor => {
                write!(f, "full multiplicative factor")
            }
            AdditionalWeightKey::Original => write!(f, "original"),
            AdditionalWeightKey::ThresholdCounterterm { subset_index } => {
                write!(f, "threshold_counterterm:{subset_index}")
            }
            AdditionalWeightKey::AmplitudeThresholdCounterterm {
                esurface_id,
                overlap_group,
            } => {
                write!(
                    f,
                    "threshold_counterterm:{esurface_id}:{overlap_group}"
                )
            }
        }
    }
}

impl fmt::Display for CutInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let rows = vec![
            EventSummaryRow {
                field: "graph".to_string(),
                value: self.graph_id.to_string(),
            },
            EventSummaryRow {
                field: "graph group".to_string(),
                value: self
                    .graph_group_id
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "None".to_string()),
            },
            EventSummaryRow {
                field: "orientation".to_string(),
                value: self
                    .orientation_id
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "None".to_string()),
            },
            EventSummaryRow {
                field: "cut".to_string(),
                value: self.cut_id.to_string(),
            },
            EventSummaryRow {
                field: "lmb channel id".to_string(),
                value: self
                    .lmb_channel_id
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "None".to_string()),
            },
            EventSummaryRow {
                field: "lmb channel".to_string(),
                value: format_lmb_channel_edge_ids(self.lmb_channel_edge_ids.as_deref()),
            },
        ];
        write!(f, "{}", Table::new(rows).with(Style::rounded()))
    }
}

impl<T: FloatLike> fmt::Display for GenericEvent<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let summary = vec![
            EventSummaryRow {
                field: "graph".to_string(),
                value: self.cut_info.graph_id.to_string(),
            },
            EventSummaryRow {
                field: "graph group".to_string(),
                value: self
                    .cut_info
                    .graph_group_id
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "None".to_string()),
            },
            EventSummaryRow {
                field: "orientation".to_string(),
                value: self
                    .cut_info
                    .orientation_id
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "None".to_string()),
            },
            EventSummaryRow {
                field: "cut".to_string(),
                value: self.cut_info.cut_id.to_string(),
            },
            EventSummaryRow {
                field: "lmb channel id".to_string(),
                value: self
                    .cut_info
                    .lmb_channel_id
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "None".to_string()),
            },
            EventSummaryRow {
                field: "lmb channel".to_string(),
                value: format_lmb_channel_edge_ids(self.cut_info.lmb_channel_edge_ids.as_deref()),
            },
            EventSummaryRow {
                field: "weight".to_string(),
                value: format_complex_generic(&self.weight),
            },
        ];

        writeln!(f, "{}", "Event".bold().bright_cyan())?;
        writeln!(f, "{}", Table::new(summary).with(Style::rounded()))?;

        let mut momentum_rows = format_momentum_rows(
            "IN".bold().bright_blue().to_string(),
            &self.kinematic_configuration.0,
            &self.cut_info.particle_pdgs.0,
            true,
        );
        let incoming_row_count = momentum_rows.len();
        momentum_rows.extend(format_momentum_rows(
            "OUT".bold().bright_green().to_string(),
            &self.kinematic_configuration.1,
            &self.cut_info.particle_pdgs.1,
            false,
        ));
        if !momentum_rows.is_empty() {
            momentum_rows.push(conservation_row(self));
            let check_separator_index = momentum_rows.len();
            let header_separator = (
                1,
                HorizontalLine::new('─')
                    .intersection('┼')
                    .left('├')
                    .right('┤'),
            );
            let check_separator = (
                check_separator_index,
                HorizontalLine::new('─')
                    .intersection('┼')
                    .left('├')
                    .right('┤'),
            );
            writeln!(f)?;
            writeln!(f, "{}", "Kinematics".bold().bright_blue())?;
            if incoming_row_count > 0 && incoming_row_count < check_separator_index - 1 {
                let in_out_separator = (
                    incoming_row_count + 1,
                    HorizontalLine::new('─')
                        .intersection('┼')
                        .left('├')
                        .right('┤'),
                );
                writeln!(
                    f,
                    "{}",
                    Table::new(momentum_rows).with(Style::rounded().horizontals([
                        header_separator,
                        in_out_separator,
                        check_separator,
                    ]))
                )?;
            } else {
                writeln!(
                    f,
                    "{}",
                    Table::new(momentum_rows)
                        .with(Style::rounded().horizontals([header_separator, check_separator]))
                )?;
            }
        }

        if !self.additional_weights.weights.is_empty() {
            let additional_weights = self
                .additional_weights
                .weights
                .iter()
                .map(|(key, value)| AdditionalWeightRow {
                    key: key.to_string(),
                    value: format_complex_generic(value),
                })
                .collect::<Vec<_>>();

            writeln!(f)?;
            writeln!(f, "{}", "Additional weights".bold().bright_magenta())?;
            write!(
                f,
                "{}",
                Table::new(additional_weights).with(Style::rounded())
            )?;
        }

        Ok(())
    }
}

impl<T: FloatLike> fmt::Display for GenericEventGroup<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "{}",
            format!("Event group ({} event(s))", format_count(self.len()))
                .bold()
                .bright_green()
        )?;
        for (index, event) in self.iter().enumerate() {
            writeln!(f, "  {}", format!("Event {index}:").bold().bright_yellow())?;
            let event_str = event.to_string();
            writeln!(f, "{}", indent_block(&event_str, "    "))?;
        }
        Ok(())
    }
}

impl<T: FloatLike> fmt::Display for GenericEventGroupList<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (index, group) in self.iter().enumerate() {
            if index > 0 {
                writeln!(f)?;
            }
            writeln!(
                f,
                "{}",
                format!("Event group {index}:").bold().bright_green()
            )?;
            writeln!(f, "{}", indent_block(&group.to_string(), "  "))?;
        }
        Ok(())
    }
}
