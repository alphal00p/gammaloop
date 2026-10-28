use numerica::domains::float::{Float, FloatLike, Real};

use super::*;

fn close<T: KinematicScalar>(actual: &T, expected: &T, tolerance: &T) {
    assert!((actual.clone() - expected.clone()).norm() < *tolerance);
}

fn complex_close<T: KinematicScalar>(actual: &Complex<T>, expected: &Complex<T>, tolerance: &T) {
    close(&actual.re, &expected.re, tolerance);
    close(&actual.im, &expected.im, tolerance);
}

fn rational<T: KinematicScalar>(prototype: &T, numerator: i64, denominator: usize) -> T {
    prototype.from_i64(numerator) / prototype.from_usize(denominator)
}

fn check_dirac_weyl_and_vector_projectors<T: KinematicScalar>(
    owner: &FutureNullMomentum<T>,
    prototype: &T,
) {
    let momentum = owner.as_momentum();
    let energy = &momentum.temporal.value;
    let zero = Complex::new(prototype.zero(), prototype.zero());
    let one = Complex::new(prototype.one(), prototype.zero());
    let tolerance = prototype
        .from_usize(2)
        .pow(u64::from(prototype.get_precision() - 24))
        .inv();
    // Construct the Pauli and chiral gamma matrices directly from p/E, with
    // no helicity-chart formula or spinor-energy helper in the oracle.
    let [x, y, z] = momentum
        .spatial
        .components()
        .map(|value| value.clone() / energy.clone());
    let sigma = [
        [
            Complex::new(z.clone(), prototype.zero()),
            Complex::new(x.clone(), -y.clone()),
        ],
        [Complex::new(x, y), Complex::new(-z, prototype.zero())],
    ];
    let mut slash = vec![vec![zero.clone(); 4]; 4];
    for i in 0..2 {
        for j in 0..2 {
            let diagonal = if i == j { one.clone() } else { zero.clone() };
            slash[i][j + 2] = diagonal.clone() - sigma[i][j].clone();
            slash[i + 2][j] = diagonal + sigma[i][j].clone();
        }
    }
    let normalization = energy.sqrt();
    for kind in [WavefunctionKind::U, WavefunctionKind::V] {
        let mut completeness = vec![vec![zero.clone(); 4]; 4];
        for helicity in [Helicity::PLUS, Helicity::MINUS] {
            let state = owner.wavefunction(kind, helicity).unwrap();
            let adjoint = owner.wavefunction(kind.bar(), helicity).unwrap();
            assert_eq!(adjoint, state.bar());
            assert_eq!(adjoint.bar(), state);
            let normalize = |component: &Complex<T>| {
                Complex::new(
                    component.re.clone() / normalization.clone(),
                    component.im.clone() / normalization.clone(),
                )
            };
            let spinor = state.components().iter().map(normalize).collect::<Vec<_>>();
            let bar = adjoint
                .components()
                .iter()
                .map(normalize)
                .collect::<Vec<_>>();
            let active_right = (kind == WavefunctionKind::U && helicity.is_plus())
                || (kind == WavefunctionKind::V && helicity.is_minus());
            let active = if active_right { 2 } else { 0 };
            let inactive = if active_right { 0 } else { 2 };
            // These are exact source zeros, not tolerance-based residual tests.
            for component in &state.components()[inactive..inactive + 2] {
                assert_eq!(component.re, prototype.zero());
                assert_eq!(component.im, prototype.zero());
            }
            // Each active Weyl rank-one projector is I ± sigma.n.
            for i in 0..2 {
                for j in 0..2 {
                    let diagonal = if i == j { one.clone() } else { zero.clone() };
                    let expected = if active_right {
                        diagonal + sigma[i][j].clone()
                    } else {
                        diagonal - sigma[i][j].clone()
                    };
                    complex_close(
                        &(spinor[active + i].clone() * spinor[active + j].conj()),
                        &expected,
                        &tolerance,
                    );
                }
            }
            for i in 0..4 {
                let column = (0..4).fold(zero.clone(), |sum, j| {
                    sum + slash[i][j].clone() * spinor[j].clone()
                });
                let row = (0..4).fold(zero.clone(), |sum, j| {
                    sum + bar[j].clone() * slash[j][i].clone()
                });
                complex_close(&column, &zero, &tolerance);
                complex_close(&row, &zero, &tolerance);
                for (entry, bar_component) in completeness[i].iter_mut().zip(&bar) {
                    *entry += &spinor[i] * bar_component;
                }
            }
        }
        for i in 0..4 {
            for j in 0..4 {
                complex_close(&completeness[i][j], &slash[i][j], &tolerance);
            }
        }
    }
    let vectors = [Helicity::PLUS, Helicity::MINUS].map(|helicity| {
        let state = owner
            .wavefunction(WavefunctionKind::Epsilon, helicity)
            .unwrap();
        assert_eq!(
            owner
                .wavefunction(WavefunctionKind::EpsilonBar, helicity)
                .unwrap(),
            state.bar(),
        );
        state
    });
    let direction = momentum
        .spatial
        .components()
        .map(|value| value.clone() / energy.clone());
    for (index, state) in vectors.iter().enumerate() {
        assert_eq!(state.components()[0], zero);
        let ward = (0..3).fold(zero.clone(), |sum, i| {
            sum + state.components()[i + 1].clone()
                * Complex::new(direction[i].clone(), prototype.zero())
        });
        complex_close(&ward, &zero, &tolerance);
        for (other_index, other) in vectors.iter().enumerate() {
            let inner = (1..4).fold(zero.clone(), |sum, i| {
                sum + state.components()[i].clone() * other.components()[i].conj()
            });
            complex_close(
                &inner,
                if index == other_index { &one } else { &zero },
                &tolerance,
            );
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            let projector = vectors.iter().fold(zero.clone(), |sum, state| {
                sum + state.components()[i + 1].clone() * state.components()[j + 1].conj()
            });
            let expected = (if i == j {
                prototype.one()
            } else {
                prototype.zero()
            }) - direction[i].clone() * direction[j].clone();
            complex_close(
                &projector,
                &Complex::new(expected, prototype.zero()),
                &tolerance,
            );
        }
    }
}

fn rational_family<T: KinematicScalar>(prototype: &T) {
    // Independent rational parametrization of the future light cone. Cyclic
    // axis permutations and reflections exercise different chart orientations.
    for (a, b) in [
        (rational(prototype, 1, 3), rational(prototype, 2, 5)),
        (rational(prototype, 4, 3), rational(prototype, -3, 2)),
    ] {
        let square = a.clone() * a.clone() + b.clone() * b.clone();
        let expected_energy = prototype.one() + square.clone();
        let coordinates = [
            prototype.from_usize(2) * a,
            prototype.from_usize(2) * b,
            prototype.one() - square,
        ];
        for shift in 0..3 {
            for sign in [-1, 1] {
                let owner = FutureNullMomentum::from_spatial(ThreeMomentum::new(
                    coordinates[shift].clone(),
                    coordinates[(shift + 1) % 3].clone() * prototype.from_i64(sign),
                    coordinates[(shift + 2) % 3].clone(),
                ))
                .unwrap();
                let tolerance = prototype
                    .from_usize(2)
                    .pow(u64::from(prototype.get_precision() - 20))
                    .inv();
                close(
                    &owner.as_momentum().temporal.value,
                    &expected_energy,
                    &tolerance,
                );
                check_dirac_weyl_and_vector_projectors(&owner, prototype);
            }
        }
    }
}

#[test]
fn rational_null_families_obey_independent_dirac_weyl_and_vector_identities() {
    rational_family(&0.0_f64);
    for bits in [113, 137, 233, 256, 512, 1024] {
        rational_family(&Float::with_val(bits, 0));
    }
}

fn axes<T: KinematicScalar>(prototype: &T) {
    let tolerance = prototype
        .from_usize(2)
        .pow(u64::from(prototype.get_precision() - 16))
        .inv();
    for (z, targets) in [
        (
            2,
            [[0, 0, 2, 0], [0, 2, 0, 0], [0, -2, 0, 0], [0, 0, -2, 0]],
        ),
        (
            -2,
            [[0, 0, 0, -2], [2, 0, 0, 0], [-2, 0, 0, 0], [0, 0, 0, 2]],
        ),
    ] {
        let owner = FutureNullMomentum::from_spatial(ThreeMomentum::new(
            prototype.zero(),
            prototype.zero(),
            prototype.from_i64(z),
        ))
        .unwrap();
        for ((kind, helicity), expected) in [
            (WavefunctionKind::U, Helicity::PLUS),
            (WavefunctionKind::U, Helicity::MINUS),
            (WavefunctionKind::V, Helicity::PLUS),
            (WavefunctionKind::V, Helicity::MINUS),
        ]
        .into_iter()
        .zip(targets)
        {
            let state = owner.wavefunction(kind, helicity).unwrap();
            for (actual, expected) in state.components().iter().zip(expected) {
                complex_close(
                    actual,
                    &Complex::new(prototype.from_i64(expected), prototype.zero()),
                    &tolerance,
                );
            }
            assert_eq!(
                owner.wavefunction(kind.bar(), helicity).unwrap(),
                state.bar()
            );
        }
        check_dirac_weyl_and_vector_projectors(&owner, prototype);
    }
}

#[test]
fn exact_axes_keep_all_fermion_roles_and_madgraph_phases() {
    axes(&0.0_f64);
    for bits in [113, 137, 256, 512, 1024] {
        axes(&Float::with_val(bits, 0));
    }
}

fn scaled_sources<T: KinematicScalar>(prototype: &T) {
    let huge = prototype.from_usize(2).pow(1000);
    for scale in [huge.clone(), huge.inv()] {
        let owner = FutureNullMomentum::from_spatial(ThreeMomentum::new(
            prototype.from_usize(2) * scale.clone(),
            prototype.from_usize(3) * scale.clone(),
            -prototype.from_usize(4) * scale,
        ))
        .unwrap();
        check_dirac_weyl_and_vector_projectors(&owner, prototype);
    }
}

#[test]
fn tiny_and_huge_spatial_sources_keep_normalized_physical_identities() {
    scaled_sources(&0.0_f64);
    for bits in [113, 137, 256, 512, 1024] {
        scaled_sources(&Float::with_val(bits, 0));
    }
    // 2E overflows binary64, while all physical massless spinors remain finite.
    let energy = f64::MAX * 0.75;
    assert!((2.0 * energy).is_infinite());
    let owner = FutureNullMomentum::from_spatial(ThreeMomentum::new(0.0, 0.0, energy)).unwrap();
    check_dirac_weyl_and_vector_projectors(&owner, &0.0_f64);
    assert!(matches!(
        owner
            .as_momentum()
            .wavefunction(WavefunctionKind::U, Helicity::PLUS),
        Err(WavefunctionError::NonFinite)
    ));
}

fn near_south<T: KinematicScalar>(prototype: &T) -> FutureNullMomentum<T> {
    let delta = prototype.from_usize(2).pow(30).inv();
    FutureNullMomentum::from_spatial(ThreeMomentum::new(
        prototype.from_usize(2) * delta.clone(),
        prototype.zero(),
        -prototype.one() + delta.clone() * delta,
    ))
    .unwrap()
}

fn near_south_controls<T: KinematicScalar>(prototype: &T) {
    let owner = near_south(prototype);
    let actual = owner
        .wavefunction(WavefunctionKind::U, Helicity::PLUS)
        .unwrap();
    let delta = prototype.from_usize(2).pow(30).inv();
    let expected = prototype.from_usize(2).sqrt();
    let tolerance = prototype
        .from_usize(2)
        .pow(u64::from(prototype.get_precision() - 16))
        .inv();
    // The rationalized southern chart retains the small component even when
    // binary64 rounds the source's energy and pz to the exact-axis values.
    assert_ne!(actual.components()[2].re, prototype.zero());
    complex_close(
        &actual.components()[2],
        &Complex::new(expected.clone() * &delta, prototype.zero()),
        &(tolerance.clone() * delta),
    );
    complex_close(
        &actual.components()[3],
        &Complex::new(expected, prototype.zero()),
        &tolerance,
    );
    check_dirac_weyl_and_vector_projectors(&owner, prototype);
}

#[test]
fn near_south_states_retain_nonzero_components_and_independent_projectors() {
    let owner = near_south(&0.0_f64);
    assert_ne!(owner.as_momentum().spatial.px, 0.0);
    assert_eq!(owner.as_momentum().spatial.pz, -1.0);
    near_south_controls(&0.0_f64);
    for bits in [113, 137, 256, 512, 1024] {
        near_south_controls(&Float::with_val(bits, 0));
    }
}

fn mixed_source<T: KinematicScalar>(prototype: &T) -> FutureNullMomentum<T> {
    let huge = prototype.from_usize(2).pow(600);
    FutureNullMomentum::from_spatial(ThreeMomentum::new(
        huge.clone().inv(),
        prototype.zero(),
        huge,
    ))
    .unwrap()
}

#[test]
fn mixed_scale_direction_loss_is_retryable_before_a_physical_component_is_erased() {
    let owner = mixed_source(&0.0_f64);
    assert!(matches!(
        owner.wavefunction(WavefunctionKind::U, Helicity::PLUS),
        Err(WavefunctionError::InsufficientPrecision { bits: 53, .. })
    ));
    for bits in [113, 137, 256, 512, 1024] {
        let prototype = Float::with_val(bits, 0);
        let actual = mixed_source(&prototype)
            .wavefunction(WavefunctionKind::U, Helicity::PLUS)
            .unwrap();
        assert!(actual.components()[3].re > prototype.zero());
        let scaled = actual.components()[3].re.clone() * prototype.from_usize(2).pow(900);
        let expected = prototype.from_usize(2).sqrt().inv();
        let tolerance = prototype.from_usize(2).pow(u64::from(bits - 20)).inv();
        close(&scaled, &expected, &tolerance);
    }
}

#[test]
fn final_product_loss_guard_is_an_explicit_arithmetic_boundary() {
    // This directly tests the final multiplication guard, not a claim that
    // these chart values arise from representable binary64 spatial sources.
    // A real mixed-scale source can fail the preceding chart check instead.
    let nonzero = f64::from_bits(1);
    for component in [Complex::new(nonzero, 0.0), Complex::new(0.0, nonzero)] {
        assert!(matches!(
            FutureNullMomentum::<f64>::scale_component(&component, &0.25),
            Err(WavefunctionError::InsufficientPrecision { bits: 53, .. })
        ));
    }
    assert_eq!(
        FutureNullMomentum::<f64>::scale_component(&Complex::new(0.0, 0.0), &0.25).unwrap(),
        Complex::new(0.0, 0.0)
    );
}

#[test]
fn invalid_sources_and_unsupported_helicities_have_explicit_errors() {
    assert_eq!(
        FutureNullMomentum::from_spatial(ThreeMomentum::new(0.0, -0.0, 0.0)),
        Err(FutureNullMomentumError::ZeroSpatialMomentum)
    );
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            FutureNullMomentum::from_spatial(ThreeMomentum::new(invalid, 1.0, 2.0)),
            Err(FutureNullMomentumError::Wavefunction(
                WavefunctionError::NonFinite
            ))
        );
    }
    assert_eq!(
        FutureNullMomentum::from_spatial(ThreeMomentum::new(f64::MAX, f64::MAX, 0.0)),
        Err(FutureNullMomentumError::Wavefunction(
            WavefunctionError::NonFinite
        ))
    );
    let owner = FutureNullMomentum::from_spatial(ThreeMomentum::new(-2.0, 3.0, -4.0)).unwrap();
    assert!(owner.as_momentum().temporal.value > 0.0);
    for kind in [
        WavefunctionKind::U,
        WavefunctionKind::UBar,
        WavefunctionKind::V,
        WavefunctionKind::VBar,
        WavefunctionKind::Epsilon,
        WavefunctionKind::EpsilonBar,
    ] {
        assert!(owner.wavefunction(kind, Helicity::ZERO).is_err());
    }
    assert_eq!(
        owner
            .wavefunction(WavefunctionKind::Scalar, Helicity::ZERO)
            .unwrap()
            .components(),
        &[Complex::new(1.0, 0.0)]
    );
    assert!(
        owner
            .wavefunction(WavefunctionKind::Scalar, Helicity::PLUS)
            .is_err()
    );
    let mut copied_raw = *owner.as_momentum();
    copied_raw.temporal.value = -3.0;
    assert_ne!(copied_raw, *owner.as_momentum());
    assert!(owner.as_momentum().temporal.value > 0.0);
    assert_eq!(owner.clone().into_momentum(), *owner.as_momentum());
}

#[test]
fn generic_off_shell_energy_roots_keep_both_real_and_imaginary_residuals() {
    let delta = 2.0_f64.powi(-40);
    for (offset, imaginary) in [(delta, false), (-delta, true)] {
        let generic = FourMomentum::from_args(1.0 + offset, 0.0, 0.0, 1.0);
        let residual = generic.spinor_energy_factor(Sign::Negative);
        if imaginary {
            assert!(residual.im > 0.0);
            assert_eq!(residual.re, 0.0);
        } else {
            assert!(residual.re > 0.0);
            assert_eq!(residual.im, 0.0);
        }
        let spinor = generic
            .wavefunction(WavefunctionKind::U, Helicity::PLUS)
            .unwrap();
        assert_ne!(spinor.components()[0], Complex::new(0.0, 0.0));
    }
}
