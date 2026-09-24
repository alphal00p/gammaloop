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

use super::maps::SamplingEvaluationError;
use crate::utils::{F, FloatLike};

use super::SamplingExpressionEvaluator;

/// Which positive score is used to form a multichannel partition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SamplingPartitionMode {
    /// Use the exact push-forward density of each sampling map.
    MapDensity,
    /// Product of basis-edge on-shell energies; this is a partition score,
    /// not a normalized probability density. Requires a complete LMB channel.
    Ose,
    /// Use an explicitly supplied positive singularity proxy.
    SingularityProxy,
}

type SamplingScoreCallback<T> = dyn Fn(&[T]) -> Result<Option<T>> + Send + Sync;

/// A score evaluator returning a log-density.
///
/// `None` denotes that the channel is not supported at this raw point (and is
/// therefore assigned zero partition weight).  A supported score must be
/// finite; the partition checks that at least one channel is supported.
enum SamplingScoreEvaluator<T: FloatLike> {
    Callback(Arc<SamplingScoreCallback<T>>),
    Symbolica(Box<Mutex<SamplingExpressionEvaluator>>),
}

/// A score evaluator constructed from a Rust callback or a Symbolica expression.
///
/// Callback scores remain shareable. Cloning a compiled Symbolica score copies
/// its mutable evaluation buffers into an independent worker-local mutex without
/// compiling the expression again. Logarithmic callback scores handle large
/// powers without overflowing before the log-sum-exp reduction.
pub struct SamplingScoreFunction<T: FloatLike = f64> {
    evaluator: SamplingScoreEvaluator<T>,
}

impl<T: FloatLike> Clone for SamplingScoreFunction<T> {
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

impl<T: FloatLike> fmt::Debug for SamplingScoreFunction<T> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("SamplingScoreFunction(..)")
    }
}

impl<T: FloatLike> SamplingScoreFunction<T> {
    /// Construct a score from a callback already returning `log(score)`.
    pub fn from_log_function<C>(function: C) -> Self
    where
        C: Fn(&[T]) -> Result<Option<T>> + Send + Sync + 'static,
    {
        Self {
            evaluator: SamplingScoreEvaluator::Callback(Arc::new(function)),
        }
    }

    /// Construct a score from a callback returning a positive score.
    ///
    /// A `None` value still denotes an unsupported branch. Zero may be numerical
    /// underflow and requests native reevaluation; negative values are rejected.
    /// Neither can silently invalidate the partition's positivity invariant.
    pub fn from_positive_function<C>(function: C) -> Self
    where
        C: Fn(&[T]) -> Result<Option<T>> + Send + Sync + 'static,
    {
        Self::from_log_function(move |coordinates| {
            function(coordinates)?
                .map(|score| {
                    let score = F(score);
                    if !score.0.is_finite() {
                        return Err(SamplingEvaluationError::Unrepresentable {
                            operation: "sampling callback score",
                            detail: format!("computed score is not finite: {score}"),
                        }.into());
                    }
                    if score == score.zero() {
                        return Err(SamplingEvaluationError::Unrepresentable {
                            operation: "sampling callback score",
                            detail: "declared-positive callback evaluated to zero at the active precision".into(),
                        }.into());
                    }
                    if score < score.zero() {
                        return Err(eyre!(
                        "sampling channel score must be finite and strictly positive, got {score}"
                    ));
                    }
                    Ok(score.ln().0)
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
        Self::from_compiled_symbolica_positive_expression(
            SamplingExpressionEvaluator::new_positive_proxy(expression, parameters)?,
        )
    }

    /// Bind an owned warmup program to this score's native precision.
    ///
    /// The supplied program is prepared by
    /// [`SamplingExpressionEvaluator::new_positive_proxy`]. Its native programs
    /// and mutable buffers are retained without recompilation; cloning the
    /// neutral program before binding gives every worker an independent mutex.
    pub(crate) fn from_compiled_symbolica_positive_expression(
        evaluator: SamplingExpressionEvaluator,
    ) -> Result<Self> {
        if evaluator.output_count() != 1 {
            return Err(eyre!(
                "sampling Symbolica score requires exactly one output, got {}",
                evaluator.output_count()
            ));
        }
        Ok(Self {
            evaluator: SamplingScoreEvaluator::Symbolica(Box::new(Mutex::new(evaluator))),
        })
    }

    pub(crate) fn evaluate(&self, raw_coordinates: &[T]) -> Result<Option<T>> {
        let evaluator = match &self.evaluator {
            SamplingScoreEvaluator::Callback(callback) => return callback(raw_coordinates),
            SamplingScoreEvaluator::Symbolica(evaluator) => evaluator,
        };
        let mut evaluator = evaluator
            .lock()
            .map_err(|_| eyre!("sampling Symbolica score evaluator mutex was poisoned"))?;
        let value = evaluator
            .evaluate(raw_coordinates)?
            .into_iter()
            .next()
            .ok_or_else(|| eyre!("sampling Symbolica score evaluator returned no output"))?;
        let real = F(SamplingExpressionEvaluator::real_component(
            value, "proxy", 0, None,
        )?);
        if real < real.zero() {
            return Err(eyre!(
                "sampling channel score must be finite and strictly positive, got {real}"
            ));
        }
        if real == real.zero() {
            // Declaring positivity is not a proof. A nonconstant or unevaluated
            // expression may underflow to zero; retry its native evaluation,
            // then fail explicitly if unresolved rather than dropping support.
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "sampling Symbolica score",
                detail: "declared-positive expression evaluated to zero at the active precision"
                    .into(),
            }
            .into());
        }
        Ok(Some(real.ln().0))
    }
}

/// One named map and its scores in the two supported partition modes.
#[derive(Clone, Debug)]
pub struct SamplingChannelScore<T: FloatLike = f64> {
    pub name: String,
    pub map_density: SamplingScoreFunction<T>,
    pub singularity_proxy: Option<SamplingScoreFunction<T>>,
}

impl<T: FloatLike> SamplingChannelScore<T> {
    /// Build a channel with an exact map density and no proxy.
    pub fn map_density(name: impl Into<String>, score: SamplingScoreFunction<T>) -> Self {
        Self {
            name: name.into(),
            map_density: score,
            singularity_proxy: None,
        }
    }

    /// Build a channel with both exact density and an explicit proxy.
    pub fn with_proxy(
        name: impl Into<String>,
        map_density: SamplingScoreFunction<T>,
        singularity_proxy: SamplingScoreFunction<T>,
    ) -> Self {
        Self {
            name: name.into(),
            map_density,
            singularity_proxy: Some(singularity_proxy),
        }
    }

    fn score(&self, mode: SamplingPartitionMode, coordinates: &[T]) -> Result<Option<T>> {
        match mode {
            SamplingPartitionMode::MapDensity => self.map_density.evaluate(coordinates),
            SamplingPartitionMode::Ose | SamplingPartitionMode::SingularityProxy => self
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
pub struct SamplingPartition<T: FloatLike = f64> {
    pub mode: SamplingPartitionMode,
    pub log_scores: Vec<Option<T>>,
    pub weights: Vec<T>,
    pub log_denominator: T,
}

impl<T: FloatLike> SamplingPartition<T> {
    /// Construct a partition over non-empty, uniquely named channels.
    pub fn new(
        mode: SamplingPartitionMode,
        channels: &[SamplingChannelScore<T>],
        raw_coordinates: &[T],
    ) -> Result<Self> {
        Self::from_log_scores(
            mode,
            raw_coordinates,
            channels.iter().map(|channel| channel.name.as_str()),
            |index| channels[index].score(mode, raw_coordinates),
        )
    }

    /// Reduce borrowed channel scores after validating names and raw coordinates.
    /// Deferring each evaluation avoids cloning stateful map or proxy programs
    /// merely to form the partition at one point. One sequential callback also
    /// permits ordinary mutable reborrows of the original-draw policy records.
    pub(crate) fn from_log_scores<'a>(
        mode: SamplingPartitionMode,
        raw_coordinates: &[T],
        channels: impl IntoIterator<Item = &'a str>,
        mut score: impl FnMut(usize) -> Result<Option<T>>,
    ) -> Result<Self> {
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
        let mut names = channels.iter().copied();
        if let Some(empty) = names.find(|name| name.trim().is_empty()) {
            return Err(eyre!(
                "sampling partition contains an empty channel name: {empty:?}"
            ));
        }
        let mut sorted_names = channels.clone();
        sorted_names.sort_unstable();
        if let Some(duplicate) = sorted_names.windows(2).find(|pair| pair[0] == pair[1]) {
            return Err(eyre!(
                "sampling partition contains duplicate channel name '{}'",
                duplicate[0]
            ));
        }

        let log_scores = channels
            .into_iter()
            .enumerate()
            .map(|(index, name)| {
                score(index).and_then(|score| {
                    score
                        .map(|log_score| {
                            if !log_score.is_finite() {
                                return Err(SamplingEvaluationError::Unrepresentable {
                                    operation: "sampling logarithmic score",
                                    detail: format!(
                                        "channel '{name}' returned non-finite log-score {log_score}"
                                    ),
                                }
                                .into());
                            }
                            Ok(log_score)
                        })
                        .transpose()
                })
            })
            .collect::<Result<Vec<_>>>()?;

        let finite_scores = log_scores
            .iter()
            .flatten()
            .cloned()
            .map(F)
            .collect::<Vec<_>>();
        let Some(maximum) = finite_scores
            .iter()
            .max_by(|a, b| a.partial_cmp(b).expect("finite scores"))
        else {
            return Err(eyre!(
                "sampling partition denominator vanishes: all {} channels have zero support",
                log_scores.len()
            ));
        };
        let denominator_scaled = finite_scores.iter().fold(maximum.zero(), |sum, score| {
            sum + F((score - maximum).0.exp())
        });
        if !denominator_scaled.0.is_finite() || denominator_scaled <= denominator_scaled.zero() {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "sampling partition denominator",
                detail: "denominator is not finite and positive after log-sum-exp".into(),
            }
            .into());
        }
        let log_denominator = maximum + denominator_scaled.ln();
        if !log_denominator.0.is_finite() {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "sampling partition denominator",
                detail: "logarithmic denominator is not finite".into(),
            }
            .into());
        }
        let mut weights = log_scores
            .iter()
            .map(|score| {
                let Some(score) = score else {
                    return Ok(maximum.zero());
                };
                let weight = F((F(score.clone()) - maximum).0.exp()) / &denominator_scaled;
                if !weight.0.is_finite() || weight <= weight.zero() {
                    return Err(SamplingEvaluationError::Unrepresentable {
                        operation: "sampling partition weight",
                        detail:
                            "positive supported weight is not representable at the active precision"
                                .into(),
                    }
                    .into());
                }
                Ok(weight)
            })
            .collect::<Result<Vec<_>>>()?;
        // Normalize once more to keep the invariant exact enough for callers
        // that inspect the weights; this does not change the stable reduction.
        let weight_sum = weights
            .iter()
            .fold(maximum.zero(), |sum, weight| sum + weight);
        if !weight_sum.0.is_finite() || weight_sum <= weight_sum.zero() {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "sampling partition weight sum",
                detail: "weights do not have a positive finite sum".into(),
            }
            .into());
        }
        for weight in &mut weights {
            *weight /= &weight_sum;
        }
        Ok(Self {
            mode,
            log_scores,
            weights: weights.into_iter().map(|weight| weight.0).collect(),
            log_denominator: log_denominator.0,
        })
    }

    pub fn weight(&self, channel: usize) -> Option<T> {
        self.weights.get(channel).cloned()
    }

    pub fn weight_sum(&self) -> T {
        self.weights
            .iter()
            .fold(F(self.log_denominator.clone()).zero(), |sum, weight| {
                sum + F(weight.clone())
            })
            .0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::initialisation::test_initialise;
    use symbolica::{atom::AtomCore, symbol};

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
    fn log_sum_exp_requires_native_precision_for_supported_tiny_weights() {
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
        let error =
            SamplingPartition::new(SamplingPartitionMode::MapDensity, &channels, &[]).unwrap_err();
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());

        fn check<T: FloatLike>(has_extended_exponent_range: bool) {
            let one = F::<T>::default().one();
            let large = one.from_usize(1000);
            let small = -&large;
            let channels = [
                SamplingChannelScore::map_density(
                    "small",
                    SamplingScoreFunction::from_log_function(move |_| Ok(Some(small.0.clone()))),
                ),
                SamplingChannelScore::map_density(
                    "large",
                    SamplingScoreFunction::from_log_function(move |_| Ok(Some(large.0.clone()))),
                ),
            ];
            let partition =
                SamplingPartition::new(SamplingPartitionMode::MapDensity, &channels, &[]);
            if !has_extended_exponent_range {
                let error = partition.unwrap_err();
                assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());
                return;
            }
            let partition = partition.unwrap();
            let tiny_weight = F(partition.weights[0].clone());
            assert!(tiny_weight > one.zero());
            assert!(
                (tiny_weight.ln() + one.from_usize(2000)).abs()
                    < one.epsilon() * one.from_usize(4096)
            );
            assert_eq!(F(partition.weight_sum()), one);
        }
        // QuadFloat improves precision but retains the binary64 exponent range.
        check::<crate::utils::QuadFloat>(false);
        check::<crate::utils::ArbPrec>(true);
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
        let channels: Vec<SamplingChannelScore> = vec![SamplingChannelScore::map_density(
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
            SamplingScoreFunction::from_positive_function(|_| Ok(Some(-1.0))),
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
        let error = SamplingScoreFunction::<f64>::from_symbolica_positive_expression(
            radius.clone() - radius.clone(),
            [radius],
        )
        .unwrap_err();
        assert!(error.to_string().contains("strictly positive"));
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_none());

        let parameter = Atom::var(symbol!("sampling_proxy_multiple_outputs"));
        let evaluator = SamplingExpressionEvaluator::new(
            [parameter.clone(), parameter.clone()],
            [parameter],
            &[],
        )
        .unwrap();
        let error =
            SamplingScoreFunction::<f64>::from_compiled_symbolica_positive_expression(evaluator)
                .unwrap_err();
        assert!(error.to_string().contains("exactly one output"));
    }

    #[test]
    fn native_symbolica_proxy_recovers_underflow_without_dropping_support() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("sampling_native_proxy_x"));
        let expression = (-x.clone()).exp();
        let program =
            SamplingExpressionEvaluator::new_positive_proxy(expression, [x.clone()]).unwrap();
        let proxy = SamplingScoreFunction::<f64>::from_compiled_symbolica_positive_expression(
            program.clone(),
        )
        .unwrap();
        let error = proxy.evaluate(&[1000.0]).unwrap_err();
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());
        let callback = SamplingScoreFunction::from_positive_function(|coordinates: &[f64]| {
            Ok(Some((-coordinates[0]).exp()))
        });
        let error = callback.evaluate(&[1000.0]).unwrap_err();
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());

        fn check<T: FloatLike>(
            program: &SamplingExpressionEvaluator,
            parameter: &Atom,
            has_extended_exponent_range: bool,
        ) {
            let one = F::<T>::default().one();
            let x = one.from_usize(1000);
            let proxy = SamplingScoreFunction::<T>::from_compiled_symbolica_positive_expression(
                program.clone(),
            )
            .unwrap();
            let SamplingScoreEvaluator::Symbolica(evaluator) = &proxy.evaluator else {
                panic!("expected a compiled Symbolica score")
            };
            assert!(evaluator.try_lock().is_ok());
            let callback = SamplingScoreFunction::from_positive_function(|coordinates: &[T]| {
                Ok(Some((-F(coordinates[0].clone())).0.exp()))
            });
            for result in [
                proxy.evaluate(std::slice::from_ref(&x.0)),
                callback.evaluate(std::slice::from_ref(&x.0)),
            ] {
                if has_extended_exponent_range {
                    let score = F(result.unwrap().expect("positive native support"));
                    assert!((score + &x).abs() < one.epsilon() * one.from_usize(4096));
                } else {
                    let error = result.unwrap_err();
                    assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());
                }
            }
            // Runtime negative values are invalid metadata rather than a reason
            // to keep retrying an arbitrary expression at higher precision.
            let invalid = SamplingScoreFunction::<T>::from_symbolica_positive_expression(
                parameter.clone(),
                [parameter.clone()],
            )
            .unwrap();
            let error = invalid.evaluate(&[(-one).0]).unwrap_err();
            assert!(error.to_string().contains("strictly positive"));
            assert!(error.downcast_ref::<SamplingEvaluationError>().is_none());
        }
        // Bind the same compiled payload at native precisions without reparsing
        // or optimization. The original worker's lock must remain independent.
        let SamplingScoreEvaluator::Symbolica(evaluator) = &proxy.evaluator else {
            panic!("expected a compiled Symbolica score")
        };
        let _original_guard = evaluator.lock().unwrap();
        check::<crate::utils::QuadFloat>(&program, &x, false);
        check::<crate::utils::ArbPrec>(&program, &x, true);
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
            channels.iter().map(|channel| channel.name.as_str()),
            |index| channels[index].score(mode, &[0.25]),
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
            let error = SamplingPartition::from_log_scores(mode, &coordinates, names, |_| {
                evaluated.set(evaluated.get() + 1);
                Ok(Some(0.0))
            })
            .unwrap_err();
            assert!(error.to_string().contains(diagnostic), "{error}");
            assert_eq!(evaluated.get(), 0);
        }
    }
}
