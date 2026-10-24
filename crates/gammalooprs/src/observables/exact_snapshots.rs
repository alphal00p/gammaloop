//! Strict snapshot entrypoints for the negotiated exact represented-f64 mode.
//! Bare legacy files, live filling, rescaling and spatial rebinning do not
//! silently erase or downgrade the exact payload.
//! V1 histograms have a fixed declared continuous partition or discrete ID set;
//! every batch includes every bin with implicit cut zeros already counted.
//! Dynamic observed-ID unions require a future explicit metadata/zero-extension
//! contract. This restriction does not govern independent scalar projection bins.

use super::{
    HistogramBinSnapshot, HistogramSnapshot, HistogramSnapshotKind, ObservableSnapshotBundle,
    checked_histogram_count,
    exact_moments::{
        ExactF64Moments, MAX_EXACT_OBJECTS, MAX_JSON_FRAME_BYTES, MomentOrder, ReportPolicy,
        certified_round_f64,
    },
};
use eyre::{Result, eyre};
use serde::{
    Deserialize, Deserializer, Serialize,
    de::{DeserializeSeed, MapAccess, SeqAccess, Visitor},
};
use std::{
    cell::Cell,
    collections::{BTreeMap, BTreeSet},
    fmt,
};
use symbolica::domains::rational::Rational;

impl HistogramBinSnapshot {
    pub fn validate_exact(&self, sample_count: u64, order: MomentOrder) -> Result<()> {
        let exact = self
            .exact_moments
            .as_ref()
            .ok_or_else(|| eyre!("Histogram bin is missing exact moments"))?;
        exact.validate_count(sample_count)?;
        if exact.order() != order {
            return Err(eyre!(
                "Histogram bin exact moment order differs from requested order"
            ));
        }
        if self.x_min.is_some_and(|value| !value.is_finite())
            || self.x_max.is_some_and(|value| !value.is_finite())
        {
            return Err(eyre!("Histogram bin has nonfinite coordinate metadata"));
        }
        let expected_sum = certified_round_f64(exact.sum(), ReportPolicy::ExactModeRawDiagnostic)?;
        let expected_square =
            certified_round_f64(exact.sum_squared(), ReportPolicy::ExactModeRawDiagnostic)?;
        if !self.sum_weights.is_finite()
            || !self.sum_weights_squared.is_finite()
            || self.sum_weights != expected_sum
            || self.sum_weights_squared != expected_square
        {
            return Err(eyre!(
                "Histogram raw diagnostics differ from certified exact sums or exceed finite raw range"
            ));
        }
        if sample_count == 0 && (self.entry_count != 0 || self.mitigated_fill_count != 0) {
            return Err(eyre!("Empty histogram bin has nonzero entry diagnostics"));
        }
        if let Some(centered) = self.centered_sum_weights_squared {
            let expected = if sample_count == 0 {
                0.0
            } else {
                certified_round_f64(
                    &(exact.variance_numerator() / Rational::from(sample_count)),
                    ReportPolicy::ExactModeRawDiagnostic,
                )?
            };
            if !centered.is_finite() || centered != expected {
                return Err(eyre!(
                    "Optional centered diagnostic disagrees with exact moments"
                ));
            }
        }
        Ok(())
    }

    pub fn average_exact(&self, sample_count: u64) -> Result<f64> {
        let exact = self
            .exact_moments
            .as_ref()
            .ok_or_else(|| eyre!("Histogram bin is missing exact moments"))?;
        self.validate_exact(sample_count, exact.order())?;
        Ok(exact.mean()?)
    }

    pub fn error_exact(&self, sample_count: u64) -> Result<f64> {
        let exact = self
            .exact_moments
            .as_ref()
            .ok_or_else(|| eyre!("Histogram bin is missing exact moments"))?;
        self.validate_exact(sample_count, exact.order())?;
        Ok(exact.histogram_sem()?)
    }

    pub(super) fn merge_exact_samples(
        &mut self,
        other: &Self,
        count: usize,
        other_count: usize,
    ) -> Result<()> {
        let order = self
            .exact_moments
            .as_ref()
            .or(other.exact_moments.as_ref())
            .ok_or_else(|| eyre!("Exact merge requires exact moments"))?
            .order();
        self.validate_exact(u64::try_from(count)?, order)?;
        other.validate_exact(u64::try_from(other_count)?, order)?;
        let mut exact = self.exact_moments.clone().expect("validated exact state");
        exact.merge(other.exact_moments.as_ref().expect("validated exact state"))?;
        let mut staged = self.clone();
        staged.sum_weights =
            certified_round_f64(exact.sum(), ReportPolicy::ExactModeRawDiagnostic)?;
        staged.sum_weights_squared =
            certified_round_f64(exact.sum_squared(), ReportPolicy::ExactModeRawDiagnostic)?;
        staged.entry_count = checked_histogram_count(self.entry_count, other.entry_count)?;
        staged.mitigated_fill_count =
            checked_histogram_count(self.mitigated_fill_count, other.mitigated_fill_count)?;
        // The exact payload is authoritative; no rounded M2 is manufactured.
        staged.centered_sum_weights_squared = None;
        staged.exact_moments = Some(exact);
        *self = staged;
        Ok(())
    }
}

impl HistogramSnapshot {
    pub(super) fn first_exact_order(&self) -> Option<MomentOrder> {
        self.bins
            .iter()
            .chain([&self.underflow_bin, &self.overflow_bin])
            .find_map(|bin| bin.exact_moments.as_ref().map(ExactF64Moments::order))
    }

    pub fn has_exact_moments(&self) -> bool {
        self.bins
            .iter()
            .chain([&self.underflow_bin, &self.overflow_bin])
            .any(|bin| bin.exact_moments.is_some())
    }

    pub(super) fn require_legacy(&self, operation: &str) -> Result<()> {
        if self.has_exact_moments() {
            return Err(eyre!("Unsupported exact snapshot operation: {operation}"));
        }
        Ok(())
    }

    pub fn validate_exact(&self, order: MomentOrder) -> Result<()> {
        let count = u64::try_from(self.sample_count)?;
        let object_count = self
            .bins
            .len()
            .checked_add(2)
            .ok_or_else(|| eyre!("Exact object count overflow"))?;
        if object_count > MAX_EXACT_OBJECTS {
            return Err(eyre!("Exact snapshot object limit exceeded"));
        }
        if self.x_min.is_some_and(|value| !value.is_finite())
            || self.x_max.is_some_and(|value| !value.is_finite())
        {
            return Err(eyre!("Histogram has nonfinite coordinate metadata"));
        }
        for bin in self
            .bins
            .iter()
            .chain([&self.underflow_bin, &self.overflow_bin])
        {
            bin.validate_exact(count, order)?;
        }
        if self.sample_count == 0
            && (self.statistics.in_range_entry_count != 0
                || self.statistics.nan_value_count != 0
                || self.statistics.mitigated_pair_count != 0)
        {
            return Err(eyre!("Empty histogram has nonzero entry diagnostics"));
        }
        match self.kind {
            HistogramSnapshotKind::Continuous => {
                let lower = self
                    .x_min
                    .ok_or_else(|| eyre!("Continuous exact histogram lacks lower bound"))?;
                let upper = self
                    .x_max
                    .ok_or_else(|| eyre!("Continuous exact histogram lacks upper bound"))?;
                if lower >= upper || self.bins.is_empty() {
                    return Err(eyre!("Invalid continuous exact histogram range"));
                }
                let mut previous = lower;
                for bin in &self.bins {
                    if bin.x_min != Some(previous)
                        || bin.x_max.is_none_or(|value| value <= previous)
                        || bin.bin_id.is_some()
                    {
                        return Err(eyre!(
                            "Exact histogram bins are not a contiguous ordered continuous partition"
                        ));
                    }
                    previous = bin.x_max.expect("checked upper bound");
                }
                if previous != upper {
                    return Err(eyre!("Exact histogram final edge differs from upper bound"));
                }
            }
            HistogramSnapshotKind::Discrete => {
                let mut ids = BTreeSet::new();
                for bin in &self.bins {
                    let id = bin
                        .bin_id
                        .ok_or_else(|| eyre!("Discrete exact bin lacks bin_id"))?;
                    if !ids.insert(id) {
                        return Err(eyre!("Duplicate discrete exact bin_id"));
                    }
                }
            }
        }
        Ok(())
    }

    pub fn merge_exact(&self, other: &Self, order: MomentOrder) -> Result<Self> {
        self.validate_exact(order)?;
        other.validate_exact(order)?;
        if self.discrete_ordering != other.discrete_ordering {
            return Err(eyre!("Exact histogram ordering metadata differs"));
        }
        for (left, right) in [
            (&self.underflow_bin, &other.underflow_bin),
            (&self.overflow_bin, &other.overflow_bin),
        ] {
            if !same_bin_metadata(left, right) {
                return Err(eyre!("Exact underflow/overflow metadata differs"));
            }
        }
        if self.kind == HistogramSnapshotKind::Continuous {
            if self
                .bins
                .iter()
                .zip(&other.bins)
                .any(|(left, right)| left.label != right.label)
            {
                return Err(eyre!("Exact histogram bin labels differ"));
            }
        } else {
            let bins: BTreeMap<_, _> = other.bins.iter().map(|bin| (bin.bin_id, bin)).collect();
            if self.bins.iter().any(|bin| {
                bins.get(&bin.bin_id)
                    .is_none_or(|other| !same_bin_metadata(bin, other))
            }) {
                return Err(eyre!("Exact discrete bin metadata differs"));
            }
        }
        let mut merged = self.clone();
        merged.merge_samples_in_place(other)?;
        merged.validate_exact(order)?;
        Ok(merged)
    }

    pub fn to_exact_json(&self, order: MomentOrder) -> Result<Vec<u8>> {
        self.validate_exact(order)?;
        bounded_json(&SnapshotEnvelope {
            kind: "gammaloop_exact_histogram_v1".into(),
            moment_order: order_number(order),
            payload: self.clone(),
        })
    }

    pub fn from_exact_json(bytes: &[u8], order: MomentOrder) -> Result<Self> {
        preflight_exact_objects(bytes)?;
        let envelope: SnapshotEnvelope = serde_json::from_slice(bytes)?;
        if envelope.kind != "gammaloop_exact_histogram_v1"
            || envelope.moment_order != order_number(order)
        {
            return Err(eyre!("Exact histogram envelope capability mismatch"));
        }
        envelope.payload.validate_exact(order)?;
        Ok(envelope.payload)
    }
}

impl ObservableSnapshotBundle {
    pub fn has_exact_moments(&self) -> bool {
        self.histograms
            .values()
            .any(HistogramSnapshot::has_exact_moments)
    }

    pub(super) fn require_legacy(&self, operation: &str) -> Result<()> {
        if self.has_exact_moments() {
            return Err(eyre!("Unsupported exact bundle operation: {operation}"));
        }
        Ok(())
    }

    pub fn validate_exact(&self, order: MomentOrder) -> Result<()> {
        let mut total_objects = 0usize;
        let mut common_count = None;
        for histogram in self.histograms.values() {
            total_objects = total_objects
                .checked_add(histogram.bins.len())
                .and_then(|count| count.checked_add(2))
                .ok_or_else(|| eyre!("Exact bundle object count overflow"))?;
            if total_objects > MAX_EXACT_OBJECTS {
                return Err(eyre!("Exact bundle object limit exceeded"));
            }
            histogram.validate_exact(order)?;
            if common_count.is_some_and(|count| count != histogram.sample_count) {
                return Err(eyre!("Exact bundle histogram sample counts differ"));
            }
            common_count = Some(histogram.sample_count);
        }
        Ok(())
    }

    /// An empty map is the layout-free identity. A nonempty legacy bundle never
    /// silently adopts an exact label, even if one rounded raw total is zero.
    pub fn merge_exact_in_place(&mut self, other: &Self, order: MomentOrder) -> Result<()> {
        self.validate_exact(order)?;
        other.validate_exact(order)?;
        if self.histograms.is_empty() {
            *self = other.clone();
            return Ok(());
        }
        if other.histograms.is_empty() {
            return Ok(());
        }
        if !self.histograms.keys().eq(other.histograms.keys()) {
            return Err(eyre!("Exact bundle histogram names differ"));
        }
        let histograms = self
            .histograms
            .iter()
            .map(|(name, histogram)| {
                Ok((
                    name.clone(),
                    histogram.merge_exact(&other.histograms[name], order)?,
                ))
            })
            .collect::<Result<_>>()?;
        let staged = Self { histograms };
        staged.validate_exact(order)?;
        *self = staged;
        Ok(())
    }

    pub fn to_exact_json(&self, order: MomentOrder) -> Result<Vec<u8>> {
        self.validate_exact(order)?;
        bounded_json(&BundleEnvelope {
            kind: "gammaloop_exact_observables_v1".into(),
            moment_order: order_number(order),
            payload: self.clone(),
        })
    }

    pub fn from_exact_json(bytes: &[u8], order: MomentOrder) -> Result<Self> {
        preflight_exact_objects(bytes)?;
        let envelope: BundleEnvelope = serde_json::from_slice(bytes)?;
        if envelope.kind != "gammaloop_exact_observables_v1"
            || envelope.moment_order != order_number(order)
        {
            return Err(eyre!("Exact observables envelope capability mismatch"));
        }
        envelope.payload.validate_exact(order)?;
        Ok(envelope.payload)
    }
}

fn order_number(order: MomentOrder) -> u8 {
    if order == MomentOrder::Second { 2 } else { 4 }
}
fn same_bin_metadata(left: &HistogramBinSnapshot, right: &HistogramBinSnapshot) -> bool {
    left.x_min == right.x_min
        && left.x_max == right.x_max
        && left.bin_id == right.bin_id
        && left.label == right.label
}
fn check_frame(bytes: &[u8]) -> Result<()> {
    if bytes.len() > MAX_JSON_FRAME_BYTES {
        return Err(eyre!("Exact JSON frame limit exceeded"));
    }
    Ok(())
}
fn preflight_exact_objects(bytes: &[u8]) -> Result<()> {
    check_frame(bytes)?;
    // Walk original JSON tokens before constructing any Rational or Value. In
    // particular, a later duplicate payload must not erase a large first payload
    // from the resource count. Reject all duplicate map keys at this boundary.
    let count = Cell::new(0usize);
    let mut decoder = serde_json::Deserializer::from_slice(bytes);
    ExactResourcePreflight(&count).deserialize(&mut decoder)?;
    decoder.end()?;
    Ok(())
}

#[derive(Clone, Copy)]
struct ExactResourcePreflight<'a>(&'a Cell<usize>);
impl<'de> DeserializeSeed<'de> for ExactResourcePreflight<'_> {
    type Value = ();
    fn deserialize<D: Deserializer<'de>>(self, decoder: D) -> std::result::Result<(), D::Error> {
        decoder.deserialize_any(self)
    }
}
impl<'de> Visitor<'de> for ExactResourcePreflight<'_> {
    type Value = ();
    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("bounded JSON without duplicate keys")
    }
    fn visit_unit<E: serde::de::Error>(self) -> std::result::Result<(), E> {
        Ok(())
    }
    fn visit_bool<E: serde::de::Error>(self, _: bool) -> std::result::Result<(), E> {
        Ok(())
    }
    fn visit_i64<E: serde::de::Error>(self, _: i64) -> std::result::Result<(), E> {
        Ok(())
    }
    fn visit_u64<E: serde::de::Error>(self, _: u64) -> std::result::Result<(), E> {
        Ok(())
    }
    fn visit_f64<E: serde::de::Error>(self, _: f64) -> std::result::Result<(), E> {
        Ok(())
    }
    fn visit_str<E: serde::de::Error>(self, _: &str) -> std::result::Result<(), E> {
        Ok(())
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut sequence: A) -> std::result::Result<(), A::Error> {
        while sequence.next_element_seed(self)?.is_some() {}
        Ok(())
    }
    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> std::result::Result<(), A::Error> {
        let mut keys = BTreeSet::new();
        while let Some(key) = map.next_key::<String>()? {
            if !keys.insert(key.clone()) {
                return Err(serde::de::Error::custom(
                    "duplicate JSON key in exact snapshot",
                ));
            }
            if key == "exact_moments" {
                let count = self.0.get() + 1;
                if count > MAX_EXACT_OBJECTS {
                    return Err(serde::de::Error::custom(
                        "exact object limit exceeded before rational decoding",
                    ));
                }
                self.0.set(count);
            }
            map.next_value_seed(self)?;
        }
        Ok(())
    }
}
fn bounded_json(value: &impl Serialize) -> Result<Vec<u8>> {
    let bytes = serde_json::to_vec(value)?;
    check_frame(&bytes)?;
    Ok(bytes)
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct SnapshotEnvelope {
    kind: String,
    moment_order: u8,
    payload: HistogramSnapshot,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct BundleEnvelope {
    kind: String,
    moment_order: u8,
    payload: ObservableSnapshotBundle,
}

/// Typed struct decoding already rejects duplicate fields. This map visitor
/// also rejects repeated histogram names, before BTreeMap can erase them.
pub(super) fn deserialize_unique_histograms<'de, D: Deserializer<'de>>(
    deserializer: D,
) -> std::result::Result<BTreeMap<String, HistogramSnapshot>, D::Error> {
    struct UniqueHistograms;
    impl<'de> Visitor<'de> for UniqueHistograms {
        type Value = BTreeMap<String, HistogramSnapshot>;
        fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
            formatter.write_str("a map of unique histogram names")
        }
        fn visit_map<A: MapAccess<'de>>(
            self,
            mut map: A,
        ) -> std::result::Result<Self::Value, A::Error> {
            let mut histograms = BTreeMap::new();
            while let Some(name) = map.next_key::<String>()? {
                if histograms.contains_key(&name) {
                    return Err(serde::de::Error::custom("duplicate histogram name"));
                }
                let histogram = map.next_value()?;
                histograms.insert(name, histogram);
                if histograms.len() > MAX_EXACT_OBJECTS {
                    return Err(serde::de::Error::custom(
                        "histogram map resource limit exceeded",
                    ));
                }
            }
            Ok(histograms)
        }
    }
    deserializer.deserialize_map(UniqueHistograms)
}

#[cfg(test)]
mod tests {
    use super::super::exact_moments::MAX_SAMPLE_COUNT;
    use super::super::{HistogramStatisticsSnapshot, tests::sample_histogram_snapshot};
    use super::*;

    fn set_moments(bin: &mut HistogramBinSnapshot, exact: ExactF64Moments) {
        bin.sum_weights =
            certified_round_f64(exact.sum(), ReportPolicy::ExactModeRawDiagnostic).unwrap();
        bin.sum_weights_squared =
            certified_round_f64(exact.sum_squared(), ReportPolicy::ExactModeRawDiagnostic).unwrap();
        bin.centered_sum_weights_squared = None;
        bin.exact_moments = Some(exact);
        bin.entry_count = 0;
        bin.mitigated_fill_count = 0;
    }

    fn snapshot(samples: &[f64], order: MomentOrder) -> HistogramSnapshot {
        let mut snapshot = sample_histogram_snapshot();
        snapshot.sample_count = samples.len();
        snapshot.statistics = HistogramStatisticsSnapshot {
            in_range_entry_count: samples.len(),
            nan_value_count: 0,
            mitigated_pair_count: 0,
        };
        let mut zeros = ExactF64Moments::empty(order);
        zeros.extend_zeros(samples.len() as u64).unwrap();
        for bin in snapshot
            .bins
            .iter_mut()
            .chain([&mut snapshot.underflow_bin, &mut snapshot.overflow_bin])
        {
            set_moments(bin, zeros.clone());
        }
        let mut exact = ExactF64Moments::empty(order);
        for &sample in samples {
            exact.push(sample).unwrap();
        }
        set_moments(&mut snapshot.bins[0], exact);
        snapshot.bins[0].entry_count = samples.len();
        snapshot.validate_exact(order).unwrap();
        snapshot
    }

    fn bundle(samples: &[f64], order: MomentOrder) -> ObservableSnapshotBundle {
        ObservableSnapshotBundle {
            histograms: ["nominal", "scale-variation"]
                .into_iter()
                .map(|name| (name.into(), snapshot(samples, order)))
                .collect(),
        }
    }

    #[test]
    fn strict_snapshot_merge_preserves_pair_and_three_adjacent_statistics() {
        let order = MomentOrder::Fourth;
        let left = snapshot(&[134_217_728.0], order);
        let right = snapshot(&[134_217_729.0], order);
        for merged in [
            left.merge_exact(&right, order).unwrap(),
            left.merge(&right).unwrap(),
        ] {
            assert_eq!(merged.bins[0].average_exact(2).unwrap(), 134_217_728.5);
            assert_eq!(merged.bins[0].error_exact(2).unwrap(), 0.5);
            assert_eq!(merged.underflow_bin.error_exact(2).unwrap(), 0.0);
            assert_eq!(
                merged
                    .overflow_bin
                    .exact_moments
                    .as_ref()
                    .unwrap()
                    .sample_count(),
                2
            );
        }
        let a = 1.0_f64;
        let b = a.next_up();
        let c = b.next_up();
        let partitioned = snapshot(&[a, c], order)
            .merge_exact(&snapshot(&[b], order), order)
            .unwrap();
        assert_eq!(partitioned, snapshot(&[a, b, c], order));
    }

    #[test]
    fn strict_bundle_merge_checks_every_scale_and_is_atomic() {
        let order = MomentOrder::Second;
        let left = bundle(&[1.0, 2.0], order);
        let right = bundle(&[3.0, 4.0], order);
        let mut merged = left.clone();
        merged.merge_exact_in_place(&right, order).unwrap();
        assert_eq!(merged, bundle(&[1.0, 2.0, 3.0, 4.0], order));
        let mut empty = ObservableSnapshotBundle::default();
        empty.merge_exact_in_place(&left, order).unwrap();
        assert_eq!(empty, left);
        for bad in [
            {
                let mut value = right.clone();
                value
                    .histograms
                    .get_mut("scale-variation")
                    .unwrap()
                    .overflow_bin
                    .exact_moments = None;
                value
            },
            {
                let mut value = right.clone();
                value.histograms.get_mut("scale-variation").unwrap().bins[0].sum_weights += 1.0;
                value
            },
            {
                let mut value = right.clone();
                value.histograms.remove("scale-variation");
                value
            },
            {
                let mut value = right.clone();
                value
                    .histograms
                    .get_mut("scale-variation")
                    .unwrap()
                    .underflow_bin
                    .label = Some("wrong layout".into());
                value
            },
        ] {
            let mut unchanged = left.clone();
            assert!(unchanged.merge_exact_in_place(&bad, order).is_err());
            assert_eq!(unchanged, left);
            assert!(unchanged.merge_in_place(&bad).is_err());
            assert_eq!(unchanged, left);
        }
    }

    #[test]
    fn tagged_roundtrip_rejects_old_shape_duplicates_and_partial_exact_data() {
        let order = MomentOrder::Fourth;
        let state = bundle(&[1.0, -1.0, 0.0], order);
        let bytes = state.to_exact_json(order).unwrap();
        assert_eq!(
            ObservableSnapshotBundle::from_exact_json(&bytes, order).unwrap(),
            state
        );
        assert!(serde_json::from_slice::<ObservableSnapshotBundle>(&bytes).is_err());
        assert!(ObservableSnapshotBundle::from_exact_json(&bytes, MomentOrder::Second).is_err());
        let text = String::from_utf8(bytes).unwrap();
        let duplicate = text.replacen("\"version\":1", "\"version\":1,\"version\":1", 1);
        assert!(ObservableSnapshotBundle::from_exact_json(duplicate.as_bytes(), order).is_err());
        let hist = &state.histograms["nominal"];
        let bytes = hist.to_exact_json(order).unwrap();
        assert_eq!(
            HistogramSnapshot::from_exact_json(&bytes, order).unwrap(),
            *hist
        );
        assert!(serde_json::from_slice::<HistogramSnapshot>(&bytes).is_err());
        let mut partial = state.clone();
        partial
            .histograms
            .get_mut("nominal")
            .unwrap()
            .underflow_bin
            .exact_moments = None;
        assert!(partial.to_exact_json(order).is_err());
        let raw_hist = serde_json::to_string(hist).unwrap();
        let duplicate_name = format!(
            "{{\"kind\":\"gammaloop_exact_observables_v1\",\"moment_order\":4,\"payload\":{{\"histograms\":{{\"x\":{raw_hist},\"x\":{raw_hist}}}}}}}"
        );
        assert!(
            ObservableSnapshotBundle::from_exact_json(duplicate_name.as_bytes(), order).is_err()
        );
    }

    #[test]
    fn exact_snapshot_refuses_unsupported_transforms_and_legacy_live_import() {
        let order = MomentOrder::Second;
        let state = snapshot(&[1.0, 2.0], order);
        assert!(state.clone().into_accumulator_state().is_err());
        assert!(
            bundle(&[1.0, 2.0], order)
                .into_accumulator_bundle()
                .is_err()
        );
        assert!(state.rescaled(1.0).is_err());
        assert!(state.rebin(1).is_err());
        let mut legacy = state.clone();
        for bin in legacy
            .bins
            .iter_mut()
            .chain([&mut legacy.underflow_bin, &mut legacy.overflow_bin])
        {
            bin.exact_moments = None;
        }
        assert!(state.merge_exact(&legacy, order).is_err());
        assert!(state.merge(&legacy).is_err());
        assert!(legacy.into_accumulator_state().is_ok());
    }

    #[test]
    fn exact_tiny_raw_square_can_be_zero_but_sem_remains_representable() {
        let tiny = f64::from_bits(1);
        let state = snapshot(&[tiny, 3.0 * tiny], MomentOrder::Fourth);
        assert_eq!(state.bins[0].sum_weights_squared, 0.0);
        assert_eq!(state.bins[0].error_exact(2).unwrap(), tiny);
        let restored = HistogramSnapshot::from_exact_json(
            &state.to_exact_json(MomentOrder::Fourth).unwrap(),
            MomentOrder::Fourth,
        )
        .unwrap();
        assert_eq!(restored.bins[0].error_exact(2).unwrap(), tiny);
        assert!(restored.bins[0].error_exact(1).is_err());
        let single = snapshot(&[tiny], MomentOrder::Second);
        assert!(single.bins[0].error_exact(1).is_err());
    }

    #[test]
    fn exact_validation_rejects_raw_overflow_bad_counts_and_metadata() {
        let order = MomentOrder::Second;
        let base = snapshot(&[1.0, 2.0], order);
        for bad in [
            {
                let mut value = base.clone();
                value.sample_count += 1;
                value
            },
            {
                let mut value = base.clone();
                value.bins[0].sum_weights = f64::INFINITY;
                value
            },
            {
                let mut value = base.clone();
                value.bins[0].x_max = Some(f64::NAN);
                value
            },
            {
                let mut value = base.clone();
                value.bins[0].x_min = Some(1.0);
                value
            },
            {
                let mut value = base.clone();
                value.bins[0].centered_sum_weights_squared = Some(0.0);
                value
            },
        ] {
            assert!(bad.validate_exact(order).is_err());
        }
        let mut overflow = base.clone();
        let mut exact = ExactF64Moments::empty(order);
        exact.push(1e200).unwrap();
        exact.push(-1e200).unwrap();
        overflow.bins[0].exact_moments = Some(exact);
        assert!(overflow.validate_exact(order).is_err());
        let mut different_counts = bundle(&[1.0, 2.0], order);
        different_counts
            .histograms
            .insert("scale-variation".into(), snapshot(&[1.0, 2.0, 3.0], order));
        assert!(different_counts.validate_exact(order).is_err());
    }

    #[test]
    fn exact_count_overflow_rejects_complete_bundle_without_mutation() {
        let order = MomentOrder::Second;
        let mut state = bundle(&[], order);
        for histogram in state.histograms.values_mut() {
            histogram.sample_count = usize::try_from(MAX_SAMPLE_COUNT).unwrap();
            for bin in histogram
                .bins
                .iter_mut()
                .chain([&mut histogram.underflow_bin, &mut histogram.overflow_bin])
            {
                let mut exact = ExactF64Moments::empty(order);
                exact.extend_zeros(MAX_SAMPLE_COUNT).unwrap();
                set_moments(bin, exact);
            }
        }
        state.validate_exact(order).unwrap();
        let before = state.clone();
        assert!(
            state
                .merge_exact_in_place(&bundle(&[0.0], order), order)
                .is_err()
        );
        assert_eq!(state, before);
    }

    #[test]
    fn exact_hwu_uses_certified_mean_and_sem_and_rejects_missing_uncertainty() {
        let state = snapshot(&[134_217_728.0, 134_217_729.0], MomentOrder::Second);
        let mut bytes = Vec::new();
        state.write_hwu_block(&mut bytes).unwrap();
        let text = String::from_utf8(bytes).unwrap();
        assert!(text.contains("5.00000000e-1"), "{text}");
        let mut bytes = Vec::new();
        assert!(
            snapshot(&[1.0], MomentOrder::Second)
                .write_hwu_block(&mut bytes)
                .is_err()
        );
        assert!(bytes.is_empty());
    }

    #[test]
    fn resource_preflight_counts_original_duplicate_payload_before_rationals() {
        let entries = vec!["{\"exact_moments\":{}}"; MAX_EXACT_OBJECTS + 1].join(",");
        let bytes = format!("{{\"payload\":{{\"bins\":[{entries}]}},\"payload\":{{}}}}");
        let error =
            HistogramSnapshot::from_exact_json(bytes.as_bytes(), MomentOrder::Second).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("object limit exceeded before rational decoding"),
            "{error}"
        );
        let duplicate = br#"{"unrecognized":{"same":1,"same":2}}"#;
        assert!(
            preflight_exact_objects(duplicate)
                .unwrap_err()
                .to_string()
                .contains("duplicate JSON key")
        );
    }

    #[test]
    fn mixed_exact_live_restore_rejects_before_mutating_earlier_legacy_histogram() {
        use super::super::{
            EntrySelection, EventProcessingRuntime, GroupBinContribution, HistogramObservable,
            IntegralDefinition, ObservableDefinition, Observables,
        };
        let mut runtime = EventProcessingRuntime::default();
        for name in ["a", "z"] {
            let observable = HistogramObservable {
                definition: ObservableDefinition::Integral(IntegralDefinition),
                selections: Vec::new(),
                entry_selection: EntrySelection::All,
                entry_index: 0,
                misbinning_max_normalized_distance: None,
                state: sample_histogram_snapshot()
                    .into_accumulator_state()
                    .unwrap(),
                pending_contributions: Vec::new(),
                candidate_pairs: Vec::new(),
                used_contributions: Vec::new(),
                grouped_contributions: vec![GroupBinContribution::default()],
                grouped_underflow: GroupBinContribution::default(),
                grouped_overflow: GroupBinContribution::default(),
                touched_bins: Vec::new(),
            };
            runtime
                .observables
                .insert(name.into(), Observables::Histogram(observable));
        }
        let before = runtime.snapshot_bundle();
        let mut mixed = before.clone();
        mixed.histograms.get_mut("a").unwrap().bins[0].sum_weights += 1.0;
        mixed
            .histograms
            .insert("z".into(), snapshot(&[1.0, 2.0], MomentOrder::Second));
        assert!(runtime.restore_snapshot_bundle(&mixed).is_err());
        assert_eq!(runtime.snapshot_bundle(), before);
    }
}
