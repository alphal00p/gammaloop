//! Errors at the boundary between one evaluated sample and observable statistics.

use super as observables;
use super::{
    GroupBinContribution, HistogramAccumulatorState, HistogramBinPosition, HistogramSnapshotKind,
};
use eyre::Report;
use std::collections::BTreeMap;
use std::error::Error as StdError;
use thiserror::Error;

/// A sample cannot be accumulated with an invalid outer integration weight.
///
/// See [`super::EventProcessingRuntime::accumulate_sample`] for a usage example.
#[derive(Debug, Error)]
pub enum ObservableAccumulationError {
    /// NaN and infinite coefficients cannot define observable moments.
    #[error("observable sample coefficient must be finite")]
    NonFiniteCoefficient,
    /// Configured observables have not been compiled by integrand warm-up.
    #[error("warm up the integrand before creating an observable accumulator")]
    NotInitialized,
    /// A fill operation does not match the histogram's coordinate kind.
    #[error("histogram '{title}' requires {expected:?} coordinates")]
    HistogramKind {
        title: String,
        expected: HistogramSnapshotKind,
    },
    /// A raw snapshot cannot be merged into this accumulator's definitions.
    #[error("incompatible observable snapshot: {source}")]
    IncompatibleSnapshot {
        #[source]
        source: Box<dyn StdError + Send + Sync>,
    },
}

impl From<Report> for ObservableAccumulationError {
    fn from(source: Report) -> Self {
        Self::IncompatibleSnapshot {
            source: source.into(),
        }
    }
}

impl HistogramAccumulatorState {
    /// Add all entries of one independent sample, summing correlated bin fills.
    /// Empty entries still advance the histogram's sample count.
    /// See [`Self::fill_discrete_sample`] for the discrete counterpart.
    ///
    /// # Errors
    /// Returns an error for a discrete histogram.
    ///
    /// # Examples
    /// ```
    /// use gammalooprs::observables::{HistogramAccumulatorState, ObservablePhase, ObservableValueTransform};
    /// let mut histogram = HistogramAccumulatorState::continuous(
    ///     "x".into(), "AL".into(), ObservablePhase::Real,
    ///     ObservableValueTransform::Identity, 0.0, 1.0, false, false, 1,
    /// );
    /// histogram.fill_continuous_sample(&[(0.2, 10.0), (0.3, -9.0)])?;
    /// assert_eq!(histogram.snapshot().bins[0].sum_weights_squared, 1.0);
    /// # Ok::<(), gammalooprs::observables::ObservableAccumulationError>(())
    /// ```
    pub fn fill_continuous_sample(
        &mut self,
        entries: &[(f64, f64)],
    ) -> Result<(), ObservableAccumulationError> {
        if self.kind != HistogramSnapshotKind::Continuous {
            return Err(ObservableAccumulationError::HistogramKind {
                title: self.title.clone(),
                expected: self.kind,
            });
        }
        let (min, max, bins) = (self.x_min, self.x_max, self.bins.len());
        self.fill_sample_bins(entries.iter().map(|&(value, weight)| {
            let position = if value.is_finite() {
                observables::histogram_bin_position(value, min, max, bins)
            } else {
                None
            };
            (position, weight)
        }));
        Ok(())
    }

    /// Add one sample's correlated entries to a discrete histogram.
    /// Empty entries still count as one attempted sample.
    /// See [`Self::fill_continuous_sample`] for an example of sample-level grouping.
    ///
    /// # Errors
    /// Returns an error for a continuous histogram.
    pub fn fill_discrete_sample(
        &mut self,
        entries: &[(isize, f64)],
    ) -> Result<(), ObservableAccumulationError> {
        if self.kind != HistogramSnapshotKind::Discrete {
            return Err(ObservableAccumulationError::HistogramKind {
                title: self.title.clone(),
                expected: self.kind,
            });
        }
        let min = self
            .discrete_min_bin_id
            .expect("discrete histogram missing min bin id");
        let bins = self.bins.len();
        self.fill_sample_bins(entries.iter().map(|&(bin, weight)| {
            (
                Some(observables::histogram_bin_position_discrete(bin, min, bins)),
                weight,
            )
        }));
        Ok(())
    }

    fn fill_sample_bins(
        &mut self,
        entries: impl IntoIterator<Item = (Option<HistogramBinPosition>, f64)>,
    ) {
        self.new_sample_count += 1;
        let mut grouped = BTreeMap::<HistogramBinPosition, GroupBinContribution>::new();
        for (position, weight) in entries {
            let Some(position) = position.filter(|_| weight.is_finite()) else {
                self.statistics.register_nan();
                continue;
            };
            let contribution = grouped.entry(position).or_default();
            contribution.projected_weight += weight;
            contribution.entry_count += 1;
        }
        for (position, contribution) in grouped {
            let bin = match position {
                HistogramBinPosition::Underflow => &mut self.underflow_bin,
                HistogramBinPosition::Overflow => &mut self.overflow_bin,
                HistogramBinPosition::InRange(index) => {
                    self.statistics.new_in_range_entry_count += contribution.entry_count;
                    &mut self.bins[index]
                }
            };
            bin.add_sample(contribution.projected_weight, contribution.entry_count, 0);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::observables::{
        DiscreteBinOrdering, HistogramAccumulatorState, ObservablePhase, ObservableSnapshotBundle,
        ObservableValueTransform,
    };
    use std::collections::BTreeMap;

    fn histogram() -> HistogramAccumulatorState {
        HistogramAccumulatorState::new_continuous(
            "test".into(),
            "AL".into(),
            ObservablePhase::Real,
            ObservableValueTransform::Identity,
            false,
            0.0,
            1.0,
            false,
            false,
            1,
        )
    }

    #[test]
    fn correlated_fills_are_summed_before_squaring_including_flow_bins() {
        let mut state = histogram();
        state
            .fill_continuous_sample(&[
                (-1.0, 10.0),
                (-1.0, -9.0),
                (0.25, 10.0),
                (0.25, -9.0),
                (2.0, 5.0),
                (2.0, -4.0),
            ])
            .unwrap();
        state.fill_continuous_sample(&[]).unwrap();
        let snapshot = state.snapshot();
        assert_eq!(snapshot.sample_count, 2);
        for bin in [
            &snapshot.bins[0],
            &snapshot.underflow_bin,
            &snapshot.overflow_bin,
        ] {
            assert_eq!(bin.entry_count, 2);
            assert_eq!(bin.sum_weights, 1.0);
            assert_eq!(bin.sum_weights_squared, 1.0);
            assert_eq!(bin.average(2), 0.5);
            assert_eq!(bin.error(2), 0.5);
        }
        let mut discrete = HistogramAccumulatorState::new_discrete(
            "discrete".into(),
            "AL".into(),
            ObservablePhase::Real,
            ObservableValueTransform::Identity,
            0,
            0,
            DiscreteBinOrdering::AscendingBinId,
            false,
            vec![],
        )
        .unwrap();
        discrete
            .fill_discrete_sample(&[(0, 10.0), (0, -9.0)])
            .unwrap();
        assert_eq!(discrete.snapshot().bins[0].sum_weights_squared, 1.0);
    }

    #[test]
    fn independent_strata_sum_estimates_instead_of_pooling_sample_counts() {
        let mut positive = histogram();
        positive.fill_continuous_sample(&[(0.25, 2.0)]).unwrap();
        positive.fill_continuous_sample(&[(0.25, 4.0)]).unwrap();
        let mut negative = histogram();
        for _ in 0..3 {
            negative.fill_continuous_sample(&[(0.25, -3.0)]).unwrap();
        }
        let strata = [positive, negative].map(|state| ObservableSnapshotBundle {
            histograms: BTreeMap::from([("test".into(), state.snapshot())]),
        });
        let result = ObservableSnapshotBundle::sum_independent_estimates(&strata).unwrap();
        assert_eq!(result.histograms["test"].sample_counts, [2, 3]);
        assert_eq!(result.histograms["test"].bins[0].mean, Some(0.0));
        assert_eq!(result.histograms["test"].bins[0].error, Some(1.0));
        let mut incomplete = strata[1].clone();
        incomplete.histograms.get_mut("test").unwrap().sample_count = 0;
        let unknown =
            ObservableSnapshotBundle::sum_independent_estimates([&strata[0], &incomplete]).unwrap();
        assert_eq!(unknown.histograms["test"].bins[0].mean, None);
        assert_eq!(unknown.histograms["test"].bins[0].error, None);
    }

    #[test]
    fn finalized_histograms_merge_all_samples_without_losing_checkpointed_counts() {
        let mut left = histogram();
        left.fill_continuous_sample(&[(0.25, 2.0)]).unwrap();
        left.update_results();
        let mut right = histogram();
        right.fill_continuous_sample(&[(0.25, 4.0)]).unwrap();
        right.update_results();
        let mut restored = right.snapshot().into_accumulator_state();
        left.merge_in_place(&mut restored).unwrap();
        let snapshot = left.snapshot();
        assert_eq!(snapshot.sample_count, 2);
        assert_eq!(snapshot.bins[0].sum_weights, 6.0);
        assert_eq!(snapshot.bins[0].sum_weights_squared, 20.0);
        assert_eq!(restored.snapshot().sample_count, 0);
    }
}
