use super::*;

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
        FourMomentum::from_args(4., 0., 0., 0.),
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
    use numerica::domains::float::{DoubleFloat, FloatLike, RealLike};
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
