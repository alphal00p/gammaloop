//! Range and chart regressions, with each precision reconstructed independently.

use std::array;

use numerica::domains::float::{Float, FloatLike};

use super::*;

fn close<T: KinematicScalar>(actual: &T, expected: &T, tolerance: &T) {
    assert!((actual.clone() - expected.clone()).norm() < *tolerance);
}

fn complex_close<T: KinematicScalar>(actual: &Complex<T>, expected: &Complex<T>, tolerance: &T) {
    close(&actual.re, &expected.re, tolerance);
    close(&actual.im, &expected.im, tolerance);
}

fn near_south<T: KinematicScalar>(prototype: &T) -> FourMomentum<T> {
    // Exact rational source: delta=2^-30 and p=(1+delta²,2delta,0,-1+delta²).
    // Build afresh in T; never promote a previously rounded four-vector.
    let delta = prototype.from_usize(2).pow(30).inv();
    let square = delta.clone() * delta.clone();
    FourMomentum::from_args(
        prototype.one() + square.clone(),
        prototype.from_usize(2) * delta,
        prototype.zero(),
        -prototype.one() + square,
    )
}

fn near_south_controls<T: KinematicScalar>(prototype: &T) {
    let momentum = near_south(prototype);
    let zero = prototype.zero();
    let delta = prototype.from_usize(2).pow(30).inv();
    let normalization = (prototype.one() + delta.clone() * &delta).sqrt().inv();
    let tolerance = prototype
        .from_usize(2)
        .pow(u64::from(prototype.get_precision() - 16))
        .inv();
    let positive = momentum.helicity_spinor(Sign::Positive).unwrap();
    let negative = momentum.helicity_spinor(Sign::Negative).unwrap();
    for (actual, expected) in [
        (&positive[0], delta.clone() * &normalization),
        (&positive[1], normalization.clone()),
        (&negative[0], -normalization.clone()),
        (&negative[1], delta * normalization),
    ] {
        assert_ne!(actual.re, zero);
        complex_close(
            actual,
            &Complex::new(expected.clone(), zero.clone()),
            &(tolerance.clone() * expected.norm()),
        );
    }
    let p = momentum
        .spatial
        .px
        .hypot(&momentum.spatial.py)
        .hypot(&momentum.spatial.pz);
    let sigma = [
        [
            Complex::new(momentum.spatial.pz.clone(), zero.clone()),
            Complex::new(momentum.spatial.px.clone(), -momentum.spatial.py.clone()),
        ],
        [
            Complex::new(momentum.spatial.px.clone(), momentum.spatial.py.clone()),
            Complex::new(-momentum.spatial.pz.clone(), zero.clone()),
        ],
    ];
    for (xi, eigenvalue) in [(&positive, p.clone()), (&negative, -p.clone())] {
        for i in 0..2 {
            let applied = &sigma[i][0] * &xi[0] + &sigma[i][1] * &xi[1];
            let target = &xi[i] * &eigenvalue;
            complex_close(&applied, &target, &(tolerance.clone() * xi[i].re.norm()));
        }
    }
    // Use a timelike energy for the independent slash projector: the generic
    // API intentionally keeps residual E-|p| roots of rounded null sources.
    let timelike = FourMomentum::from_args(
        prototype.from_usize(2),
        momentum.spatial.px.clone(),
        momentum.spatial.py.clone(),
        momentum.spatial.pz.clone(),
    );
    let mass =
        ((timelike.temporal.value.clone() - &p) * (timelike.temporal.value.clone() + &p)).sqrt();
    // Independent chiral slash matrix, without calling the chart implementation.
    let mut slash: [[Complex<T>; 4]; 4] =
        array::from_fn(|_| array::from_fn(|_| Complex::new(zero.clone(), zero.clone())));
    for i in 0..2 {
        for j in 0..2 {
            let energy = Complex::new(
                if i == j {
                    timelike.temporal.value.clone()
                } else {
                    zero.clone()
                },
                zero.clone(),
            );
            slash[i][j + 2] = &energy - &sigma[i][j];
            slash[i + 2][j] = energy + &sigma[i][j];
        }
    }
    for (kind, mass_sign) in [
        (WavefunctionKind::U, prototype.one()),
        (WavefunctionKind::V, -prototype.one()),
    ] {
        let mut projector: [[Complex<T>; 4]; 4] =
            array::from_fn(|_| array::from_fn(|_| Complex::new(zero.clone(), zero.clone())));
        for helicity in [Helicity::PLUS, Helicity::MINUS] {
            let state = timelike.wavefunction(kind, helicity).unwrap();
            let adjoint = timelike.wavefunction(kind.bar(), helicity).unwrap();
            assert_eq!(adjoint, state.bar());
            for (row, component) in projector.iter_mut().zip(state.components()) {
                for (entry, adjoint_component) in row.iter_mut().zip(adjoint.components()) {
                    *entry += component * adjoint_component;
                }
            }
        }
        for i in 0..4 {
            for j in 0..4 {
                let target = &slash[i][j]
                    + Complex::new(
                        if i == j {
                            mass.clone() * &mass_sign
                        } else {
                            zero.clone()
                        },
                        zero.clone(),
                    );
                complex_close(&projector[i][j], &target, &tolerance);
            }
        }
    }
    assert_ne!(momentum.transverse_basis().unwrap()[0][3], zero);
    for helicity in [Helicity::PLUS, Helicity::MINUS] {
        assert!(
            momentum
                .wavefunction(WavefunctionKind::Epsilon, helicity)
                .is_ok()
        );
    }
}

#[test]
fn rounded_south_axis_retains_nonzero_components_and_projectors() {
    let rounded = near_south(&0.0_f64);
    assert_ne!(rounded.spatial.px, 0.0);
    assert_eq!(rounded.spatial.pz, -1.0);
    near_south_controls(&0.0_f64);
    for bits in [113, 137, 256] {
        near_south_controls(&Float::with_val(bits, 0));
    }
}

fn axis_controls<T: KinematicScalar>(prototype: &T) {
    let tolerance = prototype.from_usize(2).pow(40).inv();
    for (z, plus, minus) in [
        (2, [1, 0], [0, 1]),
        (-2, [0, -1], [1, 0]),
        (0, [0, -1], [1, 0]),
    ] {
        let momentum = FourMomentum::from_args(
            prototype.from_usize(2),
            prototype.zero(),
            prototype.zero(),
            prototype.from_i64(z),
        );
        for (sign, target) in [(Sign::Positive, plus), (Sign::Negative, minus)] {
            for (actual, expected) in momentum.helicity_spinor(sign).unwrap().iter().zip(target) {
                complex_close(
                    actual,
                    &Complex::new(prototype.from_i64(expected), prototype.zero()),
                    &tolerance,
                );
            }
        }
        let [a, b] = momentum.transverse_basis().unwrap();
        for (actual, expected) in a.iter().zip([0, 1, 0, 0]) {
            close(actual, &prototype.from_i64(expected), &tolerance);
        }
        for (actual, expected) in b.iter().zip([0, 0, if z > 0 { 1 } else { -1 }, 0]) {
            close(actual, &prototype.from_i64(expected), &tolerance);
        }
        for kind in [
            WavefunctionKind::U,
            WavefunctionKind::V,
            WavefunctionKind::Epsilon,
        ] {
            for helicity in [Helicity::PLUS, Helicity::MINUS] {
                let state = momentum.wavefunction(kind, helicity).unwrap();
                assert_eq!(
                    momentum.wavefunction(kind.bar(), helicity).unwrap(),
                    state.bar()
                );
            }
        }
    }
}

#[test]
fn true_axes_and_rest_keep_their_madgraph_phases() {
    axis_controls(&0.0_f64);
    for bits in [113, 137, 256] {
        axis_controls(&Float::with_val(bits, 0));
    }
}

fn scale_controls<T: KinematicScalar>(prototype: &T) {
    let base = FourMomentum::from_args(
        prototype.from_usize(5),
        prototype.from_usize(3),
        prototype.from_usize(4),
        prototype.zero(),
    );
    let tolerance = prototype
        .from_usize(2)
        .pow(u64::from(prototype.get_precision() - 16))
        .inv();
    let huge = prototype.from_usize(2).pow(600);
    for scale in [huge.clone(), huge.inv()] {
        let momentum = base.map_ref(|component| component.clone() * scale.clone());
        let spinor_scale = scale.sqrt();
        let [a, b] = momentum.transverse_basis().unwrap();
        // Known geometric basis for the (3,4,0) direction, independent of scale.
        let target_a = [
            prototype.zero(),
            prototype.zero(),
            prototype.zero(),
            -prototype.one(),
        ];
        let target_b = [
            prototype.zero(),
            -prototype.from_usize(4) / prototype.from_usize(5),
            prototype.from_usize(3) / prototype.from_usize(5),
            prototype.zero(),
        ];
        for (actual, expected) in a.iter().zip(&target_a).chain(b.iter().zip(&target_b)) {
            close(actual, expected, &tolerance);
        }
        for kind in [
            WavefunctionKind::U,
            WavefunctionKind::UBar,
            WavefunctionKind::V,
            WavefunctionKind::VBar,
            WavefunctionKind::Epsilon,
            WavefunctionKind::EpsilonBar,
        ] {
            for helicity in [Helicity::PLUS, Helicity::MINUS] {
                let original = base.wavefunction(kind, helicity).unwrap();
                let scaled = momentum.wavefunction(kind, helicity).unwrap();
                let divisor = if kind.is_spinor() {
                    spinor_scale.clone()
                } else {
                    prototype.one()
                };
                for (actual, expected) in scaled.components().iter().zip(original.components()) {
                    let normalized = Complex::new(
                        actual.re.clone() / divisor.clone(),
                        actual.im.clone() / divisor.clone(),
                    );
                    complex_close(&normalized, expected, &tolerance);
                }
                if !kind.is_spinor() {
                    // Divide p by its scale before contraction: the test must
                    // not reintroduce an overflowing squared norm of its own.
                    let mut contraction = Complex::new(prototype.zero(), prototype.zero());
                    let mut norm = prototype.zero();
                    for (index, (component, p)) in scaled
                        .components()
                        .iter()
                        .zip(base.components())
                        .enumerate()
                    {
                        let metric = if index == 0 {
                            prototype.one()
                        } else {
                            -prototype.one()
                        };
                        contraction += component.clone()
                            * Complex::new(p.clone() * metric.clone(), prototype.zero());
                        norm += metric
                            * (component.re.clone() * component.re.clone()
                                + component.im.clone() * component.im.clone());
                    }
                    complex_close(
                        &contraction,
                        &Complex::new(prototype.zero(), prototype.zero()),
                        &tolerance,
                    );
                    close(&norm, &-prototype.one(), &tolerance);
                }
            }
        }
    }
}

#[test]
fn transverse_and_spinor_states_survive_extreme_uniform_rescaling() {
    scale_controls(&0.0_f64);
    for bits in [113, 137, 256] {
        scale_controls(&Float::with_val(bits, 0));
    }
}

fn hemisphere_scale_controls<T: KinematicScalar>(prototype: &T) {
    let normalization = prototype.from_usize(26).sqrt().inv();
    let tolerance = prototype
        .from_usize(2)
        .pow(u64::from(prototype.get_precision() - 16))
        .inv();
    let huge = prototype.from_usize(2).pow(600);
    for scale in [huge.clone(), huge.inv()] {
        for z in [-12, 12] {
            // Pythagorean direction (3,4,±12) with norm 13 gives an
            // independent exact component reference for both hemisphere charts.
            let momentum = FourMomentum::from_args(
                prototype.from_usize(13) * &scale,
                prototype.from_usize(3) * &scale,
                prototype.from_usize(4) * &scale,
                prototype.from_i64(z) * &scale,
            );
            let (upper, lower) = if z < 0 {
                (normalization.clone(), normalization.clone())
            } else {
                (
                    prototype.from_usize(5) * &normalization,
                    normalization.clone() / prototype.from_usize(5),
                )
            };
            let positive = [
                Complex::new(upper, prototype.zero()),
                Complex::new(
                    prototype.from_usize(3) * &lower,
                    prototype.from_usize(4) * lower,
                ),
            ];
            let negative = [-positive[1].conj(), positive[0].clone()];
            for (sign, expected) in [(Sign::Positive, positive), (Sign::Negative, negative)] {
                let xi = momentum.helicity_spinor(sign).unwrap();
                for (actual, expected) in xi.iter().zip(&expected) {
                    complex_close(actual, expected, &tolerance);
                }
            }
        }
    }
}

#[test]
fn both_hemisphere_charts_survive_extreme_uniform_rescaling() {
    hemisphere_scale_controls(&0.0_f64);
    for bits in [113, 137, 256] {
        hemisphere_scale_controls(&Float::with_val(bits, 0));
    }
}

#[test]
fn mixed_scales_preserve_nonzero_direction_components_or_require_rescue() {
    // Each source component is a normal binary64 number, but px/|p| is 2^-1200.
    // Losing it before multiplying a spinor by sqrt(E+|p|) would erase a final
    // component of order 2^-901, which binary64 can represent.
    let huge = 2.0_f64.powi(600);
    let rounded = FourMomentum::from_args(huge, huge.recip(), 0.0, huge);
    assert!(matches!(
        rounded.transverse_basis(),
        Err(WavefunctionError::InsufficientPrecision { bits: 53, .. })
    ));
    for sign in [Sign::Positive, Sign::Negative] {
        assert!(matches!(
            rounded.helicity_spinor(sign),
            Err(WavefunctionError::InsufficientPrecision { bits: 53, .. })
        ));
    }
    let kinds = [
        WavefunctionKind::U,
        WavefunctionKind::UBar,
        WavefunctionKind::V,
        WavefunctionKind::VBar,
        WavefunctionKind::Epsilon,
        WavefunctionKind::EpsilonBar,
    ];
    for kind in kinds {
        for helicity in [Helicity::PLUS, Helicity::MINUS] {
            assert!(matches!(
                rounded.wavefunction(kind, helicity),
                Err(WavefunctionError::InsufficientPrecision { bits: 53, .. })
            ));
        }
    }
    for bits in [113, 137, 256] {
        let prototype = Float::with_val(bits, 0);
        let huge = prototype.from_usize(2).pow(600);
        let tiny = huge.inv();
        let ratio = tiny.clone() / huge.clone();
        let momentum = FourMomentum::from_args(huge.clone(), tiny, prototype.zero(), huge);
        let tolerance = prototype.from_usize(2).pow(u64::from(bits - 16)).inv();
        let [a, _] = momentum.transverse_basis().unwrap();
        assert!(a[3] < prototype.zero());
        close(
            &(a[3].clone() / ratio.clone()),
            &-prototype.one(),
            &tolerance,
        );
        for sign in [Sign::Positive, Sign::Negative] {
            let xi = momentum.helicity_spinor(sign).unwrap();
            let index = if sign == Sign::Positive { 1 } else { 0 };
            let expected = if sign == Sign::Positive {
                prototype.one()
            } else {
                -prototype.one()
            } / prototype.from_usize(2);
            assert_ne!(xi[index].re, prototype.zero());
            close(
                &(xi[index].re.clone() / ratio.clone()),
                &expected,
                &tolerance,
            );
        }
        for kind in kinds {
            for helicity in [Helicity::PLUS, Helicity::MINUS] {
                let state = momentum.wavefunction(kind, helicity).unwrap();
                if !kind.is_spinor() {
                    assert_ne!(state.components()[3].re, prototype.zero());
                }
            }
        }
        // The small right-handed component is restored from the exact source.
        let state = momentum
            .wavefunction(WavefunctionKind::U, Helicity::PLUS)
            .unwrap();
        assert_ne!(state.components()[3].re, prototype.zero());
    }
}
#[test]
fn southern_mixed_scale_component_loss_still_requires_rescue() {
    let huge = 2.0_f64.powi(600);
    let rounded = FourMomentum::from_args(huge, huge.recip(), 0.0, -huge);
    for sign in [Sign::Positive, Sign::Negative] {
        assert!(matches!(
            rounded.helicity_spinor(sign),
            Err(WavefunctionError::InsufficientPrecision { bits: 53, .. })
        ));
    }
    assert!(matches!(
        rounded.wavefunction(WavefunctionKind::U, Helicity::PLUS),
        Err(WavefunctionError::InsufficientPrecision { bits: 53, .. })
    ));
    for bits in [113, 137, 256] {
        let prototype = Float::with_val(bits, 0);
        let huge = prototype.from_usize(2).pow(600);
        let tiny = huge.inv();
        let ratio = tiny.clone() / &huge;
        let momentum = FourMomentum::from_args(huge.clone(), tiny, prototype.zero(), -huge);
        let tolerance = prototype.from_usize(2).pow(u64::from(bits - 16)).inv();
        for (sign, index) in [(Sign::Positive, 0), (Sign::Negative, 1)] {
            let xi = momentum.helicity_spinor(sign).unwrap();
            assert_ne!(xi[index].re, prototype.zero());
            close(
                &(xi[index].re.clone() / &ratio),
                &(prototype.one() / prototype.from_usize(2)),
                &tolerance,
            );
        }
        let state = momentum
            .wavefunction(WavefunctionKind::U, Helicity::PLUS)
            .unwrap();
        assert_ne!(state.components()[2].re, prototype.zero());
    }
}

fn rotated_projectors<T: KinematicScalar>(prototype: &T) {
    let tolerance = prototype
        .from_usize(2)
        .pow(u64::from(prototype.get_precision() - 20))
        .inv();
    let c = |re: i64, im: i64| Complex::new(prototype.from_i64(re), prototype.from_i64(im));
    // Exact timelike source E=5, |p|=3, m=4. Permutations and spatial
    // reversal exercise both helicities and particle/antiparticle channels
    // with genuinely complex spinors, independently of chart formulas.
    for [x, y, z] in [[1, 2, 2], [2, 1, -2], [-2, -2, -1]] {
        let p = FourMomentum::from_args(
            prototype.from_i64(5),
            prototype.from_i64(x),
            prototype.from_i64(y),
            prototype.from_i64(z),
        );
        let sigma = [[c(z, 0), c(x, -y)], [c(x, y), c(-z, 0)]];
        let mut slash: [[Complex<T>; 4]; 4] = array::from_fn(|_| array::from_fn(|_| c(0, 0)));
        for i in 0..2 {
            for j in 0..2 {
                let energy = c(if i == j { 5 } else { 0 }, 0);
                slash[i][j + 2] = energy.clone() - sigma[i][j].clone();
                slash[i + 2][j] = energy + sigma[i][j].clone();
            }
        }
        for (kind, mass_sign) in [(WavefunctionKind::U, 1), (WavefunctionKind::V, -1)] {
            let mut projector: [[Complex<T>; 4]; 4] =
                array::from_fn(|_| array::from_fn(|_| c(0, 0)));
            for helicity in [Helicity::PLUS, Helicity::MINUS] {
                let state = p.wavefunction(kind, helicity).unwrap();
                let adjoint = p.wavefunction(kind.bar(), helicity).unwrap();
                assert_eq!(adjoint, state.bar());
                for (i, row) in projector.iter_mut().enumerate() {
                    let dirac = (0..4).fold(c(0, 0), |sum, j| {
                        sum + slash[i][j].clone() * state.components()[j].clone()
                    });
                    complex_close(
                        &dirac,
                        &(state.components()[i].clone() * c(4 * mass_sign, 0)),
                        &tolerance,
                    );
                    for (entry, adjoint_component) in row.iter_mut().zip(adjoint.components()) {
                        *entry += state.components()[i].clone() * adjoint_component.clone();
                    }
                }
            }
            for i in 0..4 {
                for j in 0..4 {
                    complex_close(
                        &projector[i][j],
                        &(slash[i][j].clone() + c(if i == j { 4 * mass_sign } else { 0 }, 0)),
                        &tolerance,
                    );
                }
            }
        }
    }
}

#[test]
fn rotated_spinor_projectors_preserve_particle_and_antiparticle_conventions() {
    rotated_projectors(&0.0_f64);
    for bits in [113, 137, 256] {
        rotated_projectors(&Float::with_val(bits, 0));
    }
}
