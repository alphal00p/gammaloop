//! Positive channel scores and numerically stable multichannel partitions.
//!
//! A map owns its integration Jacobian.  This module only combines the
//! positive scores supplied by maps into the partition of the raw sampling
//! point.  Exact map densities are the default; singularity proxies are an
//! explicit alternative and are never silently substituted for a missing map
//! density.

use std::fmt;
use std::sync::Arc;

use color_eyre::Result;
use eyre::eyre;

/// Which positive score is used to form a multichannel partition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SamplingPartitionMode {
    /// Use the exact push-forward density of each sampling map.
    MapDensity,
    /// Use an explicitly supplied positive singularity proxy.
    SingularityProxy,
}

/// A score evaluator returning a log-density.
///
/// `None` denotes that the channel is not supported at this raw point (and is
/// therefore assigned zero partition weight).  A supported score must be
/// finite; the partition checks that at least one channel is supported.
pub trait SamplingScoreEvaluator: Send + Sync {
    fn log_score(&self, raw_coordinates: &[f64]) -> Result<Option<f64>>;
}

struct ClosureScoreEvaluator<F>(F);

impl<F> SamplingScoreEvaluator for ClosureScoreEvaluator<F>
where
    F: Fn(&[f64]) -> Result<Option<f64>> + Send + Sync,
{
    fn log_score(&self, raw_coordinates: &[f64]) -> Result<Option<f64>> {
        (self.0)(raw_coordinates)
    }
}

/// A shareable score evaluator constructed from a Rust callback.
///
/// Map compilers can use this type to wrap their warmup-compiled Symbolica
/// evaluator.  The callback returns a logarithmic score so large powers are
/// handled without overflowing before the log-sum-exp reduction.
pub struct SamplingScoreFunction {
    evaluator: Arc<dyn SamplingScoreEvaluator>,
}

impl Clone for SamplingScoreFunction {
    fn clone(&self) -> Self {
        Self {
            evaluator: Arc::clone(&self.evaluator),
        }
    }
}

impl fmt::Debug for SamplingScoreFunction {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("SamplingScoreFunction(..)")
    }
}

impl SamplingScoreFunction {
    /// Construct a score from a callback already returning `log(score)`.
    pub fn from_log_function<F>(function: F) -> Self
    where
        F: Fn(&[f64]) -> Result<Option<f64>> + Send + Sync + 'static,
    {
        Self {
            evaluator: Arc::new(ClosureScoreEvaluator(function)),
        }
    }

    /// Construct a score from a callback returning a positive score.
    ///
    /// A `None` value still denotes an unsupported branch.  Zero and negative
    /// values are rejected because allowing them would invalidate the
    /// positivity invariant of the partition denominator.
    pub fn from_positive_function<F>(function: F) -> Self
    where
        F: Fn(&[f64]) -> Result<Option<f64>> + Send + Sync + 'static,
    {
        Self::from_log_function(move |coordinates| {
            function(coordinates)?
                .map(|score| {
                    if !score.is_finite() || score <= 0.0 {
                        return Err(eyre!(
                        "sampling channel score must be finite and strictly positive, got {score}"
                    ));
                    }
                    Ok(score.ln())
                })
                .transpose()
        })
    }

    fn evaluate(&self, raw_coordinates: &[f64]) -> Result<Option<f64>> {
        self.evaluator.log_score(raw_coordinates)
    }
}

/// One named map and its scores in the two supported partition modes.
#[derive(Clone, Debug)]
pub struct SamplingChannelScore {
    pub name: String,
    pub map_density: SamplingScoreFunction,
    pub singularity_proxy: Option<SamplingScoreFunction>,
}

impl SamplingChannelScore {
    /// Build a channel with an exact map density and no proxy.
    pub fn map_density(name: impl Into<String>, score: SamplingScoreFunction) -> Self {
        Self {
            name: name.into(),
            map_density: score,
            singularity_proxy: None,
        }
    }

    /// Build a channel with both exact density and an explicit proxy.
    pub fn with_proxy(
        name: impl Into<String>,
        map_density: SamplingScoreFunction,
        singularity_proxy: SamplingScoreFunction,
    ) -> Self {
        Self {
            name: name.into(),
            map_density,
            singularity_proxy: Some(singularity_proxy),
        }
    }

    fn score(&self, mode: SamplingPartitionMode, coordinates: &[f64]) -> Result<Option<f64>> {
        match mode {
            SamplingPartitionMode::MapDensity => self.map_density.evaluate(coordinates),
            SamplingPartitionMode::SingularityProxy => self
                .singularity_proxy
                .as_ref()
                .ok_or_else(|| {
                    eyre!(
                        "channel '{}' has no singularity_proxy; refusing to substitute its map density",
                        self.name
                    )
                })?
                .evaluate(coordinates),
        }
    }
}

/// The normalized channel weights and stable denominator at one raw point.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingPartition {
    pub mode: SamplingPartitionMode,
    pub log_scores: Vec<Option<f64>>,
    pub weights: Vec<f64>,
    pub log_denominator: f64,
}

impl SamplingPartition {
    /// Construct a partition over non-empty, uniquely named channels.
    pub fn new(
        mode: SamplingPartitionMode,
        channels: &[SamplingChannelScore],
        raw_coordinates: &[f64],
    ) -> Result<Self> {
        if channels.is_empty() {
            return Err(eyre!("sampling partition requires at least one channel"));
        }
        if raw_coordinates
            .iter()
            .any(|coordinate| !coordinate.is_finite())
        {
            return Err(eyre!("raw sampling coordinates must be finite"));
        }
        let mut names = channels.iter().map(|channel| channel.name.as_str());
        if let Some(empty) = names.find(|name| name.trim().is_empty()) {
            return Err(eyre!(
                "sampling partition contains an empty channel name: {empty:?}"
            ));
        }
        let mut sorted_names = channels
            .iter()
            .map(|channel| channel.name.as_str())
            .collect::<Vec<_>>();
        sorted_names.sort_unstable();
        if let Some(duplicate) = sorted_names.windows(2).find(|pair| pair[0] == pair[1]) {
            return Err(eyre!(
                "sampling partition contains duplicate channel name '{}'",
                duplicate[0]
            ));
        }

        let log_scores = channels
            .iter()
            .map(|channel| {
                channel.score(mode, raw_coordinates).and_then(|score| {
                    score
                        .map(|log_score| {
                            if !log_score.is_finite() {
                                return Err(eyre!(
                                    "channel '{}' returned a non-finite logarithmic score {log_score}",
                                    channel.name
                                ));
                            }
                            Ok(log_score)
                        })
                        .transpose()
                })
            })
            .collect::<Result<Vec<_>>>()?;

        let finite_scores = log_scores.iter().flatten().copied().collect::<Vec<_>>();
        let Some(&maximum) = finite_scores.iter().max_by(|a, b| a.total_cmp(b)) else {
            return Err(eyre!(
                "sampling partition denominator vanishes: all {} channels have zero support",
                channels.len()
            ));
        };
        let denominator_scaled = finite_scores
            .iter()
            .map(|score| (score - maximum).exp())
            .sum::<f64>();
        if !denominator_scaled.is_finite() || denominator_scaled <= 0.0 {
            return Err(eyre!(
                "sampling partition denominator is not finite and positive after log-sum-exp"
            ));
        }
        let log_denominator = maximum + denominator_scaled.ln();
        if !log_denominator.is_finite() {
            return Err(eyre!(
                "sampling partition logarithmic denominator is not finite"
            ));
        }
        let mut weights = log_scores
            .iter()
            .map(|score| score.map_or(0.0, |score| (score - log_denominator).exp()))
            .collect::<Vec<_>>();
        // Normalize once more to keep the invariant exact enough for callers
        // that inspect the weights; this does not change the stable reduction.
        let weight_sum = weights.iter().sum::<f64>();
        if !weight_sum.is_finite() || weight_sum <= 0.0 {
            return Err(eyre!(
                "sampling partition weights do not have a positive finite sum"
            ));
        }
        for weight in &mut weights {
            *weight /= weight_sum;
        }
        Ok(Self {
            mode,
            log_scores,
            weights,
            log_denominator,
        })
    }

    pub fn weight(&self, channel: usize) -> Option<f64> {
        self.weights.get(channel).copied()
    }

    pub fn weight_sum(&self) -> f64 {
        self.weights.iter().sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn positive(score: f64) -> SamplingScoreFunction {
        SamplingScoreFunction::from_positive_function(move |_| Ok(Some(score)))
    }

    #[test]
    fn map_density_weights_are_positive_and_sum_to_one() {
        let channels = vec![
            SamplingChannelScore::map_density("left", positive(2.0)),
            SamplingChannelScore::map_density("right", positive(3.0)),
        ];
        let partition =
            SamplingPartition::new(SamplingPartitionMode::MapDensity, &channels, &[0.25, -3.0])
                .unwrap();
        assert!((partition.weight_sum() - 1.0).abs() < 1.0e-14);
        assert!(partition.weights.iter().all(|weight| *weight > 0.0));
        assert!((partition.weights[0] - 0.4).abs() < 1.0e-14);
        assert!((partition.weights[1] - 0.6).abs() < 1.0e-14);
    }

    #[test]
    fn log_sum_exp_handles_extreme_scores_without_overflow() {
        let channels = vec![
            SamplingChannelScore::map_density(
                "small",
                SamplingScoreFunction::from_log_function(|_| Ok(Some(-1000.0))),
            ),
            SamplingChannelScore::map_density(
                "large",
                SamplingScoreFunction::from_log_function(|_| Ok(Some(1000.0))),
            ),
        ];
        let partition =
            SamplingPartition::new(SamplingPartitionMode::MapDensity, &channels, &[]).unwrap();
        assert_eq!(partition.weights, vec![0.0, 1.0]);
        assert!(partition.log_denominator.is_finite());
    }

    #[test]
    fn unsupported_channels_are_allowed_if_another_channel_has_support() {
        let channels = vec![
            SamplingChannelScore::map_density(
                "surface",
                SamplingScoreFunction::from_log_function(|_| Ok(None)),
            ),
            SamplingChannelScore::map_density("fallback", positive(1.0)),
        ];
        let partition =
            SamplingPartition::new(SamplingPartitionMode::MapDensity, &channels, &[]).unwrap();
        assert_eq!(partition.weights, vec![0.0, 1.0]);
    }

    #[test]
    fn all_zero_support_is_rejected() {
        let channels = vec![SamplingChannelScore::map_density(
            "surface",
            SamplingScoreFunction::from_log_function(|_| Ok(None)),
        )];
        let error =
            SamplingPartition::new(SamplingPartitionMode::MapDensity, &channels, &[]).unwrap_err();
        assert!(error.to_string().contains("denominator vanishes"));
    }

    #[test]
    fn proxy_mode_is_explicit_and_does_not_fallback() {
        let channels = vec![SamplingChannelScore::map_density("surface", positive(1.0))];
        let error = SamplingPartition::new(SamplingPartitionMode::SingularityProxy, &channels, &[])
            .unwrap_err();
        assert!(error.to_string().contains("no singularity_proxy"));
    }

    #[test]
    fn non_positive_scores_are_rejected() {
        let channels = vec![SamplingChannelScore::map_density(
            "bad",
            SamplingScoreFunction::from_positive_function(|_| Ok(Some(0.0))),
        )];
        let error =
            SamplingPartition::new(SamplingPartitionMode::MapDensity, &channels, &[]).unwrap_err();
        assert!(error.to_string().contains("strictly positive"));
    }
}
