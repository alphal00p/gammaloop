use super::*;
use numerica::domains::float::{DoubleFloat, FloatLike, Real, RealLike};

fn c(re: f64, im: f64) -> Complex<f64> {
    Complex::new(re, im)
}
fn close(a: Complex<f64>, b: Complex<f64>) {
    assert!(
        (a.re - b.re).abs() < 2e-12 && (a.im - b.im).abs() < 2e-12,
        "{a:?} != {b:?}"
    );
}
fn state(p: &FourMomentum<f64>, kind: WavefunctionKind, h: Helicity) -> Vec<Complex<f64>> {
    p.wavefunction(kind, h).unwrap().into_components()
}

#[test]
fn madgraph_massless_axis_and_complex_spinor_fixtures() {
    // Independent reference values retained from GammaLoop's MadGraph controls.
    let minus_z = FourMomentum::from_args(2., 0., 0., -2.);
    for (kind, h, target) in [
        (WavefunctionKind::U, Helicity::PLUS, [0., 0., 0., -2.]),
        (WavefunctionKind::U, Helicity::MINUS, [2., 0., 0., 0.]),
        (WavefunctionKind::V, Helicity::PLUS, [-2., 0., 0., 0.]),
        (WavefunctionKind::V, Helicity::MINUS, [0., 0., 0., 2.]),
    ] {
        for (a, b) in state(&minus_z, kind, h).into_iter().zip(target) {
            close(a, c(b, 0.));
        }
    }
    let y = FourMomentum::from_args(4., 0., 4., 0.);
    for (kind, h, target) in [
        (
            WavefunctionKind::U,
            Helicity::PLUS,
            [c(0., 0.), c(0., 0.), c(2., 0.), c(0., 2.)],
        ),
        (
            WavefunctionKind::U,
            Helicity::MINUS,
            [c(0., 2.), c(2., 0.), c(0., 0.), c(0., 0.)],
        ),
        (
            WavefunctionKind::V,
            Helicity::PLUS,
            [c(0., -2.), c(-2., 0.), c(0., 0.), c(0., 0.)],
        ),
        (
            WavefunctionKind::V,
            Helicity::MINUS,
            [c(0., 0.), c(0., 0.), c(-2., 0.), c(0., -2.)],
        ),
    ] {
        for (a, b) in state(&y, kind, h).into_iter().zip(target) {
            close(a, b);
        }
    }
    let q = std::f64::consts::FRAC_1_SQRT_2;
    for (h, x) in [(Helicity::PLUS, -q), (Helicity::MINUS, q)] {
        let e = state(&minus_z, WavefunctionKind::Epsilon, h);
        for (a, b) in e
            .into_iter()
            .zip([c(0., 0.), c(x, 0.), c(0., q), c(0., 0.)])
        {
            close(a, b);
        }
    }
}

#[test]
fn vector_transversality_norm_and_longitudinal_completeness() {
    for p in [
        FourMomentum::from_args(5., 0., 0., 3.),
        FourMomentum::from_args(5., 1., 2., 2.),
        // Rest frame: the natural frame for a massive vector decay.
        FourMomentum::from_args(4., 0., 0., 0.),
    ] {
        let pv = [p.temporal.value, p.spatial.px, p.spatial.py, p.spatial.pz];
        let states = [Helicity::MINUS, Helicity::ZERO, Helicity::PLUS]
            .map(|h| state(&p, WavefunctionKind::Epsilon, h));
        for (i, a) in states.iter().enumerate() {
            let dot = (0..4).fold(c(0., 0.), |s, k| {
                s + a[k] * c(if k == 0 { pv[k] } else { -pv[k] }, 0.)
            });
            close(dot, c(0., 0.));
            for (j, b) in states.iter().enumerate() {
                let dot = (0..4).fold(c(0., 0.), |s, k| {
                    s + a[k] * b[k].conj() * c(if k == 0 { 1. } else { -1. }, 0.)
                });
                close(dot, c(if i == j { -1. } else { 0. }, 0.));
            }
        }
        for i in 0..4 {
            for j in 0..4 {
                let sum = states.iter().fold(c(0., 0.), |s, e| s + e[i] * e[j].conj());
                let metric = if i != j {
                    0.
                } else if i == 0 {
                    1.
                } else {
                    -1.
                };
                close(sum, c(-metric + pv[i] * pv[j] / p.mass_squared(), 0.));
            }
        }
    }
}

#[test]
fn rest_frame_longitudinal_is_small_momentum_limit_along_plus_z() {
    let rest = state(
        &FourMomentum::from_args(4., 0., 0., 0.),
        WavefunctionKind::Epsilon,
        Helicity::ZERO,
    );
    assert_eq!(rest, [c(0., 0.), c(0., 0.), c(0., 0.), c(1., 0.)]);
    let pz: f64 = 1e-13;
    let moving = state(
        &FourMomentum::from_args((16. + pz * pz).sqrt(), 0., 0., pz),
        WavefunctionKind::Epsilon,
        Helicity::ZERO,
    );
    for (a, b) in moving.into_iter().zip(rest) {
        close(a, b);
    }
}

#[test]
fn longitudinal_polarization_retains_tiny_momentum_direction() {
    for magnitude in [1.0e-200_f64, f64::MIN_POSITIVE, f64::from_bits(1)] {
        assert_eq!(magnitude * magnitude, 0.0);
        for (px, pz, expected_x, expected_z) in [
            (magnitude, 0.0, 1.0, 0.0),
            (-magnitude, 0.0, -1.0, 0.0),
            (0.0, magnitude, 0.0, 1.0),
            (0.0, -magnitude, 0.0, -1.0),
        ] {
            let momentum = FourMomentum::from_args(4.0, px, 0.0, pz);
            let polarization = state(&momentum, WavefunctionKind::Epsilon, Helicity::ZERO);
            assert_eq!(
                polarization,
                [
                    c(magnitude / 4.0, 0.0),
                    c(expected_x, 0.0),
                    c(0.0, 0.0),
                    c(expected_z, 0.0),
                ],
                "momentum {momentum:?}"
            );
            assert_eq!(
                state(&momentum, WavefunctionKind::EpsilonBar, Helicity::ZERO),
                polarization
            );
        }
    }
}

#[test]
fn massive_spinor_completeness_and_dirac_adjoint() {
    // In the chiral basis slash(p) has E-sigma.p and E+sigma.p off-diagonal.
    let p = FourMomentum::from_args(5., 1., 2., 2.);
    let mass = 4.;
    let mut slash = [[c(0., 0.); 4]; 4];
    let sigma = [[c(2., 0.), c(1., -2.)], [c(1., 2.), c(-2., 0.)]];
    for i in 0..2 {
        for j in 0..2 {
            let e = c(if i == j { 5. } else { 0. }, 0.);
            slash[i][j + 2] = e - sigma[i][j];
            slash[i + 2][j] = e + sigma[i][j];
        }
    }
    for (kind, sgn) in [(WavefunctionKind::U, 1.), (WavefunctionKind::V, -1.)] {
        let states = [Helicity::PLUS, Helicity::MINUS].map(|h| p.wavefunction(kind, h).unwrap());
        let mut sum = [[c(0., 0.); 4]; 4];
        for s in states {
            let adj = s.bar();
            assert_eq!(adj.bar(), s);
            for i in 0..4 {
                let lhs = (0..4).fold(c(0., 0.), |a, j| a + slash[i][j] * s.components()[j]);
                close(lhs, s.components()[i] * c(sgn * mass, 0.));
                for (j, entry) in sum[i].iter_mut().enumerate() {
                    *entry += s.components()[i] * adj.components()[j];
                }
            }
        }
        for i in 0..4 {
            for j in 0..4 {
                close(
                    sum[i][j],
                    slash[i][j] + c(if i == j { sgn * mass } else { 0. }, 0.),
                );
            }
        }
    }
}

#[test]
fn scalar_adjoint_and_all_direct_barred_states() {
    let p = FourMomentum::from_args(5., 1., 2., 2.);
    assert_eq!(
        state(&p, WavefunctionKind::Scalar, Helicity::ZERO),
        vec![c(1., 0.)]
    );
    let scalar = Wavefunction::from_components(WavefunctionKind::Scalar, vec![c(2., 3.)]).unwrap();
    assert_eq!(scalar.bar().components(), &[c(2., -3.)]);
    for kind in [
        WavefunctionKind::Epsilon,
        WavefunctionKind::U,
        WavefunctionKind::V,
    ] {
        for h in [Helicity::MINUS, Helicity::PLUS] {
            let s = p.wavefunction(kind, h).unwrap();
            assert_eq!(p.wavefunction(kind.bar(), h).unwrap(), s.bar());
            assert_eq!(s.bar().bar(), s);
        }
    }
    // Preserve GammaLoop's zero-spatial-momentum spinor branch and imaginary root.
    let rest = FourMomentum::from_args(4., 0., 0., 0.);
    assert_eq!(
        rest.helicity_spinor(Sign::Positive),
        [c(0., 0.), c(-1., 0.)]
    );
    assert!(
        rest.wavefunction(WavefunctionKind::U, Helicity::PLUS)
            .is_ok()
    );
    close(
        FourMomentum::from_args(1., 0., 0., 2.).spinor_energy_factor(Sign::Negative),
        c(0., 1.),
    );
}

#[test]
fn invalid_states_have_explicit_errors() {
    let p = FourMomentum::from_args(5., 0., 0., 3.);
    assert_eq!(
        p.wavefunction(WavefunctionKind::U, Helicity::ZERO),
        Err(WavefunctionError::Helicity(WavefunctionKind::U))
    );
    assert!(
        p.wavefunction(WavefunctionKind::Scalar, Helicity::PLUS)
            .is_err()
    );
    for p in [
        FourMomentum::from_args(2., 0., 0., 2.),
        FourMomentum::from_args(0., 0., 0., 0.),
        FourMomentum::from_args(1., 0., 0., 2.),
    ] {
        assert_eq!(
            p.wavefunction(WavefunctionKind::Epsilon, Helicity::ZERO),
            Err(WavefunctionError::UndefinedLongitudinal)
        );
    }
    assert_eq!(
        FourMomentum::from_args(f64::NAN, 0., 0., 1.)
            .wavefunction(WavefunctionKind::Epsilon, Helicity::PLUS),
        Err(WavefunctionError::NonFinite)
    );
    assert!(Wavefunction::from_components(WavefunctionKind::U, vec![c(0., 0.)]).is_err());
    assert!(
        Wavefunction::from_components(WavefunctionKind::Scalar, vec![c(f64::INFINITY, 0.)])
            .is_err()
    );
}

#[test]
fn shared_numerica_precision_stays_generic() {
    let zero = DoubleFloat::default();
    let p = FourMomentum::from_args(
        zero.from_usize(5),
        zero.from_usize(1),
        zero.from_usize(2),
        zero.from_usize(2),
    );
    let p64 = FourMomentum::from_args(5., 1., 2., 2.);
    for kind in [
        WavefunctionKind::Epsilon,
        WavefunctionKind::U,
        WavefunctionKind::V,
    ] {
        let state = p.wavefunction(kind, Helicity::PLUS).unwrap();
        let target = p64.wavefunction(kind, Helicity::PLUS).unwrap();
        for (a, b) in state.components().iter().zip(target.components()) {
            close(c(a.re.to_f64(), a.im.to_f64()), *b);
        }
    }
}

fn test_rotations() -> [Rotation<f64>; 4] {
    [
        Rotation::quarter_turn(Axis::X),
        Rotation::quarter_turn(Axis::Y),
        Rotation::quarter_turn(Axis::Z),
        Rotation::euler(0.3, -0.2, 0.7),
    ]
}

/// Rotate the spatial real and imaginary parts of a vector state.
fn rotated_vector(rotation: &Rotation<f64>, state: &[Complex<f64>]) -> Vec<Complex<f64>> {
    let part = |f: fn(&Complex<f64>) -> f64| {
        rotation.rotate_three(&ThreeMomentum::new(
            f(&state[1]),
            f(&state[2]),
            f(&state[3]),
        ))
    };
    let (re, im) = (part(|c| c.re), part(|c| c.im));
    vec![state[0], c(re.px, im.px), c(re.py, im.py), c(re.pz, im.pz)]
}

#[test]
fn rest_vector_states_follow_the_rotated_quantization_axis() {
    let rest = FourMomentum::from_args(4., 0., 0., 0.);
    for rotation in test_rotations() {
        for h in [Helicity::MINUS, Helicity::ZERO, Helicity::PLUS] {
            let unrotated = state(&rest, WavefunctionKind::Epsilon, h);
            let expected = rotated_vector(&rotation, &unrotated);
            for kind in [WavefunctionKind::Epsilon, WavefunctionKind::EpsilonBar] {
                let transported = rest
                    .wavefunction_after_rotations(kind, h, &[rotation])
                    .unwrap()
                    .into_components();
                let expected = if kind == WavefunctionKind::Epsilon {
                    expected.clone()
                } else {
                    expected.iter().map(Complex::conj).collect()
                };
                for (a, b) in transported.into_iter().zip(expected) {
                    close(a, b);
                }
            }
        }
    }
    // The longitudinal axis +z follows each quarter turn explicitly.
    for (axis, target) in [
        (Axis::X, [0., 0., -1., 0.]),
        (Axis::Y, [0., 1., 0., 0.]),
        (Axis::Z, [0., 0., 0., 1.]),
    ] {
        let transported = rest
            .wavefunction_after_rotations(
                WavefunctionKind::Epsilon,
                Helicity::ZERO,
                &[Rotation::quarter_turn(axis)],
            )
            .unwrap()
            .into_components();
        for (a, b) in transported.into_iter().zip(target) {
            close(a, c(b, 0.));
        }
    }
    assert_eq!(
        rest.wavefunction_after_rotations(WavefunctionKind::Epsilon, Helicity::ZERO, &[]),
        rest.wavefunction(WavefunctionKind::Epsilon, Helicity::ZERO)
    );
}

#[test]
fn transported_moving_states_differ_from_rotated_momentum_states_by_a_phase() {
    let p = FourMomentum::from_args(5., 1., 2., 2.);
    for rotation in test_rotations() {
        let rotated = rotation.rotate_four(&p);
        for (kind, helicities) in [
            (
                WavefunctionKind::Epsilon,
                &[Helicity::MINUS, Helicity::ZERO, Helicity::PLUS][..],
            ),
            (WavefunctionKind::U, &[Helicity::MINUS, Helicity::PLUS][..]),
            (WavefunctionKind::V, &[Helicity::MINUS, Helicity::PLUS][..]),
            (
                WavefunctionKind::UBar,
                &[Helicity::MINUS, Helicity::PLUS][..],
            ),
            (
                WavefunctionKind::VBar,
                &[Helicity::MINUS, Helicity::PLUS][..],
            ),
        ] {
            for &h in helicities {
                let transported = rotated
                    .wavefunction_after_rotations(kind, h, &[rotation])
                    .unwrap()
                    .into_components();
                let direct = state(&rotated, kind, h);
                let norm: f64 = direct.iter().map(|z| z.norm_squared()).sum();
                let overlap = direct
                    .iter()
                    .zip(&transported)
                    .fold(c(0., 0.), |sum, (d, t)| sum + d.conj() * *t);
                let phase = c(overlap.re / norm, overlap.im / norm);
                close(c(phase.norm_squared(), 0.), c(1., 0.));
                for (t, d) in transported.into_iter().zip(direct) {
                    close(t, phase * d);
                }
            }
        }
    }
}

#[test]
fn transported_states_compose_rotations_in_order() {
    let [x, y, z, euler] = test_rotations();
    for p in [
        FourMomentum::from_args(4., 0., 0., 0.),
        FourMomentum::from_args(5., 1., 2., 2.),
    ] {
        let rotated = euler.rotate_four(&z.rotate_four(&x.rotate_four(&p)));
        for (kind, h) in [
            (WavefunctionKind::Epsilon, Helicity::ZERO),
            (WavefunctionKind::Epsilon, Helicity::PLUS),
            (WavefunctionKind::U, Helicity::PLUS),
            (WavefunctionKind::VBar, Helicity::MINUS),
        ] {
            let stepwise = p
                .wavefunction(kind, h)
                .unwrap()
                .rotate(&x)
                .unwrap()
                .rotate(&z)
                .unwrap()
                .rotate(&euler)
                .unwrap();
            let transported = rotated
                .wavefunction_after_rotations(kind, h, &[x, z, euler])
                .unwrap();
            for (a, b) in transported.components().iter().zip(stepwise.components()) {
                close(*a, *b);
            }
        }
    }
    // A full turn about y is the identity on vectors and minus one on spinors.
    let rest = FourMomentum::from_args(4., 0., 0., 0.);
    let turns = [y; 4];
    let spinor = rest
        .wavefunction(WavefunctionKind::U, Helicity::PLUS)
        .unwrap();
    let turned = rest
        .wavefunction_after_rotations(WavefunctionKind::U, Helicity::PLUS, &turns)
        .unwrap();
    for (a, b) in turned.components().iter().zip(spinor.components()) {
        close(*a, -*b);
    }
}

#[test]
fn transported_states_reject_nonfinite_rotation_components() {
    let rest = FourMomentum::from_args(4.0, 0.0, 0.0, 0.0);
    for angle in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let rotation = Rotation::euler(angle, 0.0, 0.0);
        for (kind, helicity) in [
            (WavefunctionKind::Epsilon, Helicity::ZERO),
            (WavefunctionKind::EpsilonBar, Helicity::PLUS),
            (WavefunctionKind::U, Helicity::PLUS),
            (WavefunctionKind::VBar, Helicity::MINUS),
        ] {
            let state = rest.wavefunction(kind, helicity).unwrap();
            assert_eq!(state.rotate(&rotation), Err(WavefunctionError::NonFinite));
            assert_eq!(
                rest.wavefunction_after_rotations(kind, helicity, &[rotation]),
                Err(WavefunctionError::NonFinite),
            );
        }
    }
}

#[test]
fn transported_spinors_retain_the_working_precision() {
    let zero = DoubleFloat::default();
    let momentum = FourMomentum::from_args(zero.from_usize(4), zero, zero, zero);
    let quarter = Rotation::quarter_turn(Axis::X);
    let euler = Rotation::euler(zero.pi() / zero.from_usize(2), zero, zero);
    let tolerance = zero
        .from_usize(2)
        .pow(u64::from(zero.get_precision()))
        .inv()
        * zero.from_usize(256);
    for kind in [WavefunctionKind::U, WavefunctionKind::VBar] {
        let initial = momentum.wavefunction(kind, Helicity::PLUS).unwrap();
        let exact = initial.rotate(&quarter).unwrap();
        let trigonometric = initial.rotate(&euler).unwrap();
        for (a, b) in exact.components().iter().zip(trigonometric.components()) {
            assert!((a.re - b.re).norm() <= tolerance);
            assert!((a.im - b.im).norm() <= tolerance);
        }
        let full_turn = momentum
            .wavefunction_after_rotations(kind, Helicity::PLUS, &[quarter; 4])
            .unwrap();
        for (a, b) in full_turn.components().iter().zip(initial.components()) {
            assert!((a.re + b.re).norm() <= tolerance);
            assert!((a.im + b.im).norm() <= tolerance);
        }
    }
}
