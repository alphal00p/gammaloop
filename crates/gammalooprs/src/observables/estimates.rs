//! Final estimates from independent populations, keeping their raw sample counts separate.

use super::{
    HistogramBinSnapshot, HistogramSnapshot, HistogramSnapshotKind, ObservableSnapshotBundle,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use thiserror::Error;

/// A sum of independent estimates; no fictitious pooled sample count is assigned.
/// See [`ObservableSnapshotBundle::sum_independent_estimates`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct IndependentEstimate {
    pub mean: Option<f64>,
    pub error: Option<f64>,
}

impl IndependentEstimate {
    fn from_bins<'a>(bins: impl IntoIterator<Item = (&'a HistogramBinSnapshot, usize)>) -> Self {
        let mut mean = Some(0.0);
        let mut error = Some(0.0_f64);
        for (bin, count) in bins {
            mean = mean.and_then(|sum| (count > 0).then(|| sum + bin.average(count)));
            error = error.and_then(|sum| (count > 1).then(|| sum.hypot(bin.error(count))));
        }
        Self { mean, error }
    }
}

/// Bin estimates in the first input snapshot's order, with counts for each population.
/// Coordinates and labels are supplied by that compatible raw snapshot.
/// See [`ObservableSnapshotBundle::sum_independent_estimates`].
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct HistogramEstimates {
    pub sample_counts: Vec<usize>,
    pub bins: Vec<IndependentEstimate>,
    pub underflow: IndependentEstimate,
    pub overflow: IndependentEstimate,
}

impl HistogramEstimates {
    fn from_histograms(histograms: &[&HistogramSnapshot]) -> Self {
        let first = histograms[0];
        let discrete_bins = histograms
            .iter()
            .map(|histogram| {
                histogram
                    .bins
                    .iter()
                    .filter_map(|bin| bin.bin_id.map(|id| (id, bin)))
                    .collect::<BTreeMap<_, _>>()
            })
            .collect::<Vec<_>>();
        let bins = first
            .bins
            .iter()
            .enumerate()
            .map(|(index, bin)| {
                IndependentEstimate::from_bins(histograms.iter().enumerate().map(
                    |(population, histogram)| {
                        let other = match first.kind {
                            HistogramSnapshotKind::Discrete => {
                                let id = bin
                                    .bin_id
                                    .expect("validated discrete histogram has bin identifiers");
                                discrete_bins[population][&id]
                            }
                            HistogramSnapshotKind::Continuous => &histogram.bins[index],
                        };
                        (other, histogram.sample_count)
                    },
                ))
            })
            .collect();
        Self {
            sample_counts: histograms
                .iter()
                .map(|histogram| histogram.sample_count)
                .collect(),
            bins,
            underflow: IndependentEstimate::from_bins(
                histograms
                    .iter()
                    .map(|histogram| (&histogram.underflow_bin, histogram.sample_count)),
            ),
            overflow: IndependentEstimate::from_bins(
                histograms
                    .iter()
                    .map(|histogram| (&histogram.overflow_bin, histogram.sample_count)),
            ),
        }
    }
}

/// Presentation estimates, distinct from mergeable/reconstructible raw snapshots.
/// See [`ObservableSnapshotBundle::sum_independent_estimates`].
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ObservableEstimates {
    pub histograms: BTreeMap<String, HistogramEstimates>,
}

/// Independent estimates must describe the same observables and bins.
/// See [`ObservableSnapshotBundle::sum_independent_estimates`].
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ObservableEstimateError {
    #[error("independent estimates have different observable names")]
    ObservableNames,
    #[error("independent estimates have incompatible histogram '{0}'")]
    HistogramLayout(String),
}

impl ObservableSnapshotBundle {
    /// Sum independent population means and variances, without pooling their counts.
    ///
    /// For stratification, each input is one stratum's conditional-weight statistics.
    /// Merge worker batches within each stratum before calling this operation.
    /// A population with no attempts makes the mean unavailable; fewer than two
    /// attempts makes the error unavailable. Covariances between bins are not implied.
    ///
    /// # Errors
    /// Returns an error for different observable names or incompatible bin layouts.
    ///
    /// # Examples
    /// ```
    /// use gammalooprs::observables::ObservableSnapshotBundle;
    /// let strata = [ObservableSnapshotBundle::default(), ObservableSnapshotBundle::default()];
    /// let estimates = ObservableSnapshotBundle::sum_independent_estimates(&strata)?;
    /// assert!(estimates.histograms.is_empty());
    /// # Ok::<(), gammalooprs::observables::ObservableEstimateError>(())
    /// ```
    pub fn sum_independent_estimates<'a>(
        populations: impl IntoIterator<Item = &'a Self>,
    ) -> Result<ObservableEstimates, ObservableEstimateError> {
        let populations = populations.into_iter().collect::<Vec<_>>();
        let Some(first) = populations.first() else {
            return Ok(ObservableEstimates::default());
        };
        if populations
            .iter()
            .any(|population| !population.histograms.keys().eq(first.histograms.keys()))
        {
            return Err(ObservableEstimateError::ObservableNames);
        }
        let mut histograms = BTreeMap::new();
        for (name, reference) in &first.histograms {
            let inputs = populations
                .iter()
                .map(|population| &population.histograms[name])
                .collect::<Vec<_>>();
            if inputs
                .iter()
                .any(|histogram| !reference.has_compatible_layout(histogram))
            {
                return Err(ObservableEstimateError::HistogramLayout(name.clone()));
            }
            histograms.insert(name.clone(), HistogramEstimates::from_histograms(&inputs));
        }
        Ok(ObservableEstimates { histograms })
    }
}
