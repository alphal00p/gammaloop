//! Range and chart regressions, with each precision reconstructed independently.

use numerica::domains::float::{Float, FloatLike, Real};

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

#[test]
fn rounded_south_axis_requires_fresh_wider_reconstruction() {
    let rounded = near_south(&0.0_f64);
    assert_ne!(rounded.spatial.px, 0.0);
    assert_eq!(rounded.spatial.pz, -1.0);
    assert_ne!(rounded.transverse_basis().unwrap()[0][3], 0.0);
    for helicity in [Helicity::PLUS, Helicity::MINUS] {
        assert!(
            rounded
                .wavefunction(WavefunctionKind::Epsilon, helicity)
                .is_ok()
        );
    }
    for sign in [Sign::Positive, Sign::Negative] {
        assert!(matches!(
            rounded.helicity_spinor(sign),
            Err(WavefunctionError::InsufficientPrecision { bits: 53, .. })
        ));
    }
    for kind in [
        WavefunctionKind::U,
        WavefunctionKind::UBar,
        WavefunctionKind::V,
        WavefunctionKind::VBar,
    ] {
        for helicity in [Helicity::PLUS, Helicity::MINUS] {
            assert!(matches!(
                rounded.wavefunction(kind, helicity),
                Err(WavefunctionError::InsufficientPrecision { bits: 53, .. })
            ));
        }
    }
    for bits in [113, 137, 256] {
        let prototype = Float::with_val(bits, 0);
        let momentum = near_south(&prototype);
        let delta = prototype.from_usize(2).pow(30).inv();
        let normalization = (prototype.one() + delta.clone() * delta.clone())
            .sqrt()
            .inv();
        // This chart subtracts quantities separated at the 60th binary place.
        // Its expected cancellation loss is allowed, but the small nonzero
        // component and the nonaxis phase must both be resolved.
        let tolerance = prototype.from_usize(2).pow(u64::from(bits - 72)).inv();
        let positive = momentum.helicity_spinor(Sign::Positive).unwrap();
        let negative = momentum.helicity_spinor(Sign::Negative).unwrap();
        for (actual, expected) in [
            (&positive[0], delta.clone() * normalization.clone()),
            (&positive[1], normalization.clone()),
            (&negative[0], -normalization.clone()),
            (&negative[1], delta * normalization),
        ] {
            complex_close(
                actual,
                &Complex::new(expected, prototype.zero()),
                &tolerance,
            );
        }
        assert!(positive[0].re > prototype.zero());
        assert!(positive[1].re > prototype.zero());
        for kind in [
            WavefunctionKind::U,
            WavefunctionKind::UBar,
            WavefunctionKind::V,
            WavefunctionKind::VBar,
        ] {
            for helicity in [Helicity::PLUS, Helicity::MINUS] {
                assert!(momentum.wavefunction(kind, helicity).is_ok());
            }
        }
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
                        norm = norm
                            + metric
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
        let mut slash: [[Complex<T>; 4]; 4] =
            std::array::from_fn(|_| std::array::from_fn(|_| c(0, 0)));
        for i in 0..2 {
            for j in 0..2 {
                let energy = c(if i == j { 5 } else { 0 }, 0);
                slash[i][j + 2] = energy.clone() - sigma[i][j].clone();
                slash[i + 2][j] = energy + sigma[i][j].clone();
            }
        }
        for (kind, mass_sign) in [(WavefunctionKind::U, 1), (WavefunctionKind::V, -1)] {
            let mut projector: [[Complex<T>; 4]; 4] =
                std::array::from_fn(|_| std::array::from_fn(|_| c(0, 0)));
            for helicity in [Helicity::PLUS, Helicity::MINUS] {
                let state = p.wavefunction(kind, helicity).unwrap();
                let adjoint = p.wavefunction(kind.bar(), helicity).unwrap();
                assert_eq!(adjoint, state.bar());
                for i in 0..4 {
                    let dirac = (0..4).fold(c(0, 0), |sum, j| {
                        sum + slash[i][j].clone() * state.components()[j].clone()
                    });
                    complex_close(
                        &dirac,
                        &(state.components()[i].clone() * c(4 * mass_sign, 0)),
                        &tolerance,
                    );
                    for j in 0..4 {
                        projector[i][j] +=
                            state.components()[i].clone() * adjoint.components()[j].clone();
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
