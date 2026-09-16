//! Native integration-grid construction, including a zero-dimensional point.

use symbolica::numerical_integration::ContinuousGrid;

use crate::utils::F;

/// Retain GammaLoop's point-domain sample for zero-loop integrands.
///
/// Numerica 3's constructor rejects zero dimensions, although its public grid
/// representation and sampling, statistics, cloning and merging operations
/// support the empty product: `Sample::Continuous(1, [])`. Validate all other
/// settings through the native constructor, then remove its temporary dimension
/// before any sampling or persistence. This does not add a dummy integration
/// coordinate or change the RNG stream. Replace this adapter when upstream
/// exposes an explicit point-domain constructor.
pub(super) fn continuous_grid(
    n_dims: usize,
    n_bins: usize,
    min_samples_for_update: usize,
    bin_number_evolution: Option<Vec<usize>>,
    train_on_avg: bool,
) -> Result<ContinuousGrid<F<f64>>, String> {
    let mut grid = ContinuousGrid::new(
        n_dims.max(1),
        n_bins,
        min_samples_for_update,
        bin_number_evolution,
        train_on_avg,
    )?;
    if n_dims == 0 {
        grid.continuous_dimensions.clear();
    }
    Ok(grid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::numerical_integration::{Grid, MonteCarloRng, Sample};

    fn assert_point_sample(sample: &Sample<F<f64>>) {
        let Sample::Continuous(weight, coordinates) = sample else {
            panic!("point domain changed sample shape: {sample:?}");
        };
        assert_eq!(*weight, F(1.0));
        assert!(coordinates.is_empty());
    }

    #[test]
    fn point_domain_has_unit_weight_empty_coordinates_and_no_rng_draws() {
        let mut grid = continuous_grid(0, 8, 10, None, false).unwrap();
        let mut rng = MonteCarloRng::new(123, 0);
        let untouched_rng = rng.clone();
        let mut sample = Sample::Continuous(F(7.0), vec![F(0.5)]);
        grid.sample(&mut rng, &mut sample);
        assert_point_sample(&sample);
        assert_eq!(rng, untouched_rng);
        assert_eq!(grid.probe(&[]).unwrap(), F(1.0));
        assert!(grid.probe(&[None]).is_err());
    }

    #[test]
    fn point_domain_uses_native_statistics_clone_and_merge() {
        let mut grid = continuous_grid(0, 8, 10, None, false).unwrap();
        let sample = Sample::Continuous(F(1.0), vec![]);
        grid.add_training_sample(&sample, F(2.0)).unwrap();
        let mut other = grid.clone_without_samples();
        other.add_training_sample(&sample, F(4.0)).unwrap();
        grid.merge(&other).unwrap();
        grid.update(F(0.5));
        assert_eq!(grid.accumulator.processed_samples, 2);
        assert_eq!(grid.accumulator.avg, F(3.0));
        assert!(grid.continuous_dimensions.is_empty());
        assert!(
            grid.merge(&continuous_grid(1, 8, 10, None, false).unwrap())
                .is_err()
        );
    }

    #[test]
    fn point_domain_roundtrip_preserves_native_grid_shape() {
        let mut grid = Grid::Continuous(continuous_grid(0, 8, 10, None, false).unwrap());
        let sample = Sample::Continuous(F(1.0), vec![]);
        grid.add_training_sample(&sample, F(2.0)).unwrap();
        let bytes = serde_json::to_vec(&grid).unwrap();
        let mut restored: Grid<F<f64>> = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(serde_json::to_vec(&restored).unwrap(), bytes);
        let mut rng = MonteCarloRng::new(456, 0);
        let untouched_rng = rng.clone();
        let mut restored_sample = Sample::new();
        restored.sample(&mut rng, &mut restored_sample);
        assert_point_sample(&restored_sample);
        assert_eq!(rng, untouched_rng);
        assert_eq!(
            restored
                .probe(&symbolica::numerical_integration::Probe::Continuous(vec![]))
                .unwrap(),
            F(1.0)
        );
    }

    #[test]
    fn positive_dimension_delegates_to_native_grid_unchanged() {
        let mut wrapped = continuous_grid(3, 8, 10, Some(vec![8, 16]), true).unwrap();
        let mut native: ContinuousGrid<F<f64>> =
            ContinuousGrid::new(3, 8, 10, Some(vec![8, 16]), true).unwrap();
        assert_eq!(
            serde_json::to_vec(&wrapped).unwrap(),
            serde_json::to_vec(&native).unwrap()
        );
        let mut wrapped_rng = MonteCarloRng::new(789, 0);
        let mut native_rng = wrapped_rng.clone();
        let mut wrapped_sample = Sample::new();
        let mut native_sample = Sample::new();
        wrapped.sample(&mut wrapped_rng, &mut wrapped_sample);
        native.sample(&mut native_rng, &mut native_sample);
        assert_eq!(
            serde_json::to_vec(&wrapped_sample).unwrap(),
            serde_json::to_vec(&native_sample).unwrap()
        );
        assert_eq!(wrapped_rng, native_rng);
    }

    #[test]
    fn point_domain_does_not_bypass_native_settings_validation() {
        for dimensions in [0, 1] {
            assert!(continuous_grid(dimensions, 0, 10, None, false).is_err());
            assert!(continuous_grid(dimensions, 8, 10, Some(vec![]), false).is_err());
            assert!(continuous_grid(dimensions, 8, 10, Some(vec![8, 0]), false).is_err());
        }
    }

    #[test]
    fn point_sample_preserves_zero_loop_parameterization() {
        let mut grid = continuous_grid(0, 8, 10, None, false).unwrap();
        let mut sample = Sample::new();
        grid.sample(&mut MonteCarloRng::new(321, 0), &mut sample);
        let Sample::Continuous(weight, coordinates) = sample else {
            panic!("zero-loop grid changed sample shape");
        };
        let (momenta, jacobian) = crate::utils::global_parameterize(
            &coordinates,
            F(1.0),
            &crate::settings::runtime::ParameterizationSettings::default(),
        );
        assert!(momenta.is_empty());
        assert_eq!(weight * jacobian, F(1.0));
    }
}
