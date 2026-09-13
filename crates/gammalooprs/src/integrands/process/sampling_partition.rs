//! Positive channel scores and numerically stable multichannel partitions.
//!
//! A map owns its integration Jacobian.  This module only combines the
//! positive scores supplied by maps into the partition of the raw sampling
//! point.  Exact map densities are the default; singularity proxies are an
//! explicit alternative and are never silently substituted for a missing map
//! density.

use std::fmt;
use std::sync::{Arc, Mutex};

use color_eyre::Result;
use eyre::eyre;
use symbolica::atom::Atom;

use super::SamplingExpressionEvaluator;

/// Which positive score is used to form a multichannel partition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SamplingPartitionMode {
    /// Use the exact push-forward density of each sampling map.
    MapDensity,
    /// Use an explicitly supplied positive singularity proxy.
    SingularityProxy,
}

type SamplingScoreCallback = dyn Fn(&[f64]) -> Result<Option<f64>> + Send + Sync;

/// A score evaluator returning a log-density.
///
/// `None` denotes that the channel is not supported at this raw point (and is
/// therefore assigned zero partition weight).  A supported score must be
/// finite; the partition checks that at least one channel is supported.
enum SamplingScoreEvaluator {
    Callback(Arc<SamplingScoreCallback>),
    Symbolica(Box<Mutex<SamplingExpressionEvaluator>>),
}

/// A score evaluator constructed from a Rust callback or a Symbolica expression.
///
/// Callback scores remain shareable. Cloning a compiled Symbolica score copies
/// its mutable evaluation buffers into an independent worker-local mutex without
/// compiling the expression again. Logarithmic callback scores handle large
/// powers without overflowing before the log-sum-exp reduction.
pub struct SamplingScoreFunction {
    evaluator: SamplingScoreEvaluator,
}

impl Clone for SamplingScoreFunction {
    fn clone(&self) -> Self {
        Self {
            evaluator: match &self.evaluator {
                SamplingScoreEvaluator::Callback(callback) => {
                    SamplingScoreEvaluator::Callback(Arc::clone(callback))
                }
                SamplingScoreEvaluator::Symbolica(evaluator) => {
                    SamplingScoreEvaluator::Symbolica(Box::new(Mutex::new(
                        evaluator
                            .lock()
                            .expect("sampling Symbolica score evaluator mutex was poisoned")
                            .clone(),
                    )))
                }
            },
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
            evaluator: SamplingScoreEvaluator::Callback(Arc::new(function)),
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

    /// Construct a positive score from one Symbolica expression.
    ///
    /// The expression is compiled eagerly when this function is called and
    /// evaluated at every partition point.  A score is only accepted when it
    /// is finite, real, and strictly positive.  This keeps proxy channels on
    /// the same positivity boundary as callback-based scores while allowing
    /// process warmup code to build the evaluator once from the user supplied
    /// Symbolica expression.
    ///
    /// The evaluator is protected by a mutex because Symbolica's eager
    /// evaluator is stateful.  It is still compiled exactly once and the
    /// resulting score is safe to share between channel partitions. Worker
    /// clones own independent mutexes and evaluator buffers.
    ///
    /// Positivity alone does not establish that a proxy captures every
    /// singular feature of a map.  The graph compiler must supply and audit
    /// those expressions explicitly; declaring a positive expression does not
    /// certify its power near every intended singular feature.
    pub fn from_symbolica_positive_expression(
        expression: Atom,
        parameters: impl IntoIterator<Item = Atom>,
    ) -> Result<Self> {
        let evaluator = SamplingExpressionEvaluator::new([expression], parameters, false)?;
        Ok(Self {
            evaluator: SamplingScoreEvaluator::Symbolica(Box::new(Mutex::new(evaluator))),
        })
    }

    pub(crate) fn evaluate(&self, raw_coordinates: &[f64]) -> Result<Option<f64>> {
        let evaluator = match &self.evaluator {
            SamplingScoreEvaluator::Callback(callback) => return callback(raw_coordinates),
            SamplingScoreEvaluator::Symbolica(evaluator) => evaluator,
        };
        let mut evaluator = evaluator
            .lock()
            .map_err(|_| eyre!("sampling Symbolica score evaluator mutex was poisoned"))?;
        let values = evaluator.evaluate(raw_coordinates)?;
        let value = values
            .first()
            .ok_or_else(|| eyre!("sampling Symbolica score evaluator returned no output"))?;
        let real = value.re.0;
        let imaginary = value.im.0;
        let tolerance = 1.0e-13 * real.abs().max(1.0);
        if !real.is_finite() || !imaginary.is_finite() {
            return Err(eyre!(
                "sampling Symbolica score evaluated to a non-finite value: {value:?}"
            ));
        }
        if imaginary.abs() > tolerance {
            return Err(eyre!(
                "sampling Symbolica score has a non-negligible imaginary part {imaginary}"
            ));
        }
        if real <= 0.0 {
            return Err(eyre!(
                "sampling channel score must be finite and strictly positive, got {real}"
            ));
        }
        Ok(Some(real.ln()))
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
        Self::from_log_scores(
            mode,
            raw_coordinates,
            channels.iter().map(|channel| {
                (channel.name.as_str(), || {
                    channel.score(mode, raw_coordinates)
                })
            }),
        )
    }

    /// Reduce borrowed channel scores after validating names and raw coordinates.
    /// Deferring each evaluation avoids cloning stateful map or proxy programs
    /// merely to form the partition at one point.
    pub(crate) fn from_log_scores<'a, S>(
        mode: SamplingPartitionMode,
        raw_coordinates: &[f64],
        channels: impl IntoIterator<Item = (&'a str, S)>,
    ) -> Result<Self>
    where
        S: FnOnce() -> Result<Option<f64>>,
    {
        if raw_coordinates
            .iter()
            .any(|coordinate| !coordinate.is_finite())
        {
            return Err(eyre!("raw sampling coordinates must be finite"));
        }
        let channels = channels.into_iter().collect::<Vec<_>>();
        if channels.is_empty() {
            return Err(eyre!("sampling partition requires at least one channel"));
        }
        let mut names = channels.iter().map(|(name, _)| *name);
        if let Some(empty) = names.find(|name| name.trim().is_empty()) {
            return Err(eyre!(
                "sampling partition contains an empty channel name: {empty:?}"
            ));
        }
        let mut sorted_names = channels.iter().map(|(name, _)| *name).collect::<Vec<_>>();
        sorted_names.sort_unstable();
        if let Some(duplicate) = sorted_names.windows(2).find(|pair| pair[0] == pair[1]) {
            return Err(eyre!(
                "sampling partition contains duplicate channel name '{}'",
                duplicate[0]
            ));
        }

        let log_scores = channels
            .into_iter()
            .map(|(name, score)| {
                score().and_then(|score| {
                    score
                        .map(|log_score| {
                            if !log_score.is_finite() {
                                return Err(eyre!(
                                    "channel '{}' returned a non-finite logarithmic score {log_score}",
                                    name
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
                log_scores.len()
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
    use crate::initialisation::test_initialise;
    use symbolica::symbol;

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

    #[test]
    fn symbolica_positive_proxy_is_compiled_once_and_partitioned() {
        test_initialise().unwrap();
        let radius = Atom::var(symbol!("sampling_proxy_radius"));
        let proxy = SamplingScoreFunction::from_symbolica_positive_expression(
            radius.clone() * radius.clone() + 1,
            [radius],
        )
        .unwrap();
        let channels = vec![
            SamplingChannelScore::with_proxy("surface", positive(1.0), proxy),
            SamplingChannelScore::with_proxy("fallback", positive(1.0), positive(1.0)),
        ];
        let partition =
            SamplingPartition::new(SamplingPartitionMode::SingularityProxy, &channels, &[2.0])
                .unwrap();
        assert!((partition.weights[0] - 5.0 / 6.0).abs() < 1.0e-14);
        assert!((partition.weights[1] - 1.0 / 6.0).abs() < 1.0e-14);
    }

    #[test]
    fn symbolica_positive_proxy_rejects_non_positive_values() {
        test_initialise().unwrap();
        let radius = Atom::var(symbol!("sampling_proxy_non_positive"));
        let proxy = SamplingScoreFunction::from_symbolica_positive_expression(
            radius.clone() - radius.clone(),
            [radius],
        )
        .unwrap();
        let channels = vec![SamplingChannelScore::with_proxy(
            "surface",
            positive(1.0),
            proxy,
        )];
        let error =
            SamplingPartition::new(SamplingPartitionMode::SingularityProxy, &channels, &[2.0])
                .unwrap_err();
        assert!(error.to_string().contains("strictly positive"));
    }

    #[test]
    fn cloned_symbolica_scores_have_independent_worker_buffers() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("sampling_worker_x"));
        let y = Atom::var(symbol!("sampling_worker_y"));
        let proxy = SamplingScoreFunction::from_symbolica_positive_expression(
            x.clone() * x.clone() + y.clone() * y.clone() + 1,
            [x, y],
        )
        .unwrap();
        let workers = (0..4).map(|_| proxy.clone()).collect::<Vec<_>>();
        let SamplingScoreEvaluator::Symbolica(original) = &proxy.evaluator else {
            panic!("expected a compiled Symbolica score")
        };
        // Holding the original lock must neither block worker evaluation nor
        // alias its scratch buffers. The try-lock fails promptly for a shared
        // mutex regression, before the threaded test could deadlock.
        let _original_guard = original.lock().unwrap();
        for worker in &workers {
            let SamplingScoreEvaluator::Symbolica(evaluator) = &worker.evaluator else {
                panic!("expected a compiled Symbolica score")
            };
            assert!(evaluator.try_lock().is_ok());
        }
        std::thread::scope(|scope| {
            let handles = workers
                .into_iter()
                .enumerate()
                .map(|(worker_index, worker)| {
                    scope.spawn(move || {
                        for sample in 0..128 {
                            let x = worker_index as f64 + sample as f64 / 17.0;
                            let y = 3.0 - sample as f64 / 23.0;
                            let expected = (1.0 + x * x + y * y).ln();
                            let actual = worker.evaluate(&[x, y]).unwrap().unwrap();
                            assert!((actual - expected).abs() < 1.0e-13);
                        }
                    })
                })
                .collect::<Vec<_>>();
            for handle in handles {
                handle.join().unwrap();
            }
        });
        let callback = positive(2.0);
        let cloned = callback.clone();
        let (SamplingScoreEvaluator::Callback(original), SamplingScoreEvaluator::Callback(cloned)) =
            (&callback.evaluator, &cloned.evaluator)
        else {
            panic!("expected callback scores")
        };
        assert!(Arc::ptr_eq(original, cloned));
    }

    #[test]
    fn borrowed_partition_validates_structure_before_evaluating_scores() {
        let channels = [
            SamplingChannelScore::map_density("left", positive(2.0)),
            SamplingChannelScore::map_density("right", positive(3.0)),
        ];
        let mode = SamplingPartitionMode::MapDensity;
        let borrowed = SamplingPartition::from_log_scores(
            mode,
            &[0.25],
            channels
                .iter()
                .map(|channel| (channel.name.as_str(), || channel.score(mode, &[0.25]))),
        )
        .unwrap();
        assert_eq!(
            borrowed,
            SamplingPartition::new(mode, &channels, &[0.25]).unwrap()
        );

        let evaluated = std::cell::Cell::new(0);
        for (names, coordinates, diagnostic) in [
            (["", "valid"], [0.25], "empty channel name"),
            (["same", "same"], [0.25], "duplicate channel name"),
            (["left", "right"], [f64::NAN], "coordinates must be finite"),
        ] {
            let error = SamplingPartition::from_log_scores(
                mode,
                &coordinates,
                names.into_iter().map(|name| {
                    (name, || {
                        evaluated.set(evaluated.get() + 1);
                        Ok(Some(0.0))
                    })
                }),
            )
            .unwrap_err();
            assert!(error.to_string().contains(diagnostic), "{error}");
            assert_eq!(evaluated.get(), 0);
        }
    }
}
