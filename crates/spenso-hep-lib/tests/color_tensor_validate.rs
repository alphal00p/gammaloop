use idenso::{color::CS, representations::initialize};
use spenso::{
    algebra::complex::Complex,
    structure::abstract_index::AbstractIndex,
    tensors::data::{GetTensorData, SparseTensor},
};
use spenso_hep_lib::{su3_generator_data, su3_structure_f_data};

fn color_t()
-> SparseTensor<Complex<f64>, spenso::network::library::symbolic::ExplicitKey<AbstractIndex>> {
    initialize();
    su3_generator_data(CS.t_strct::<AbstractIndex>(3, 8))
}

fn color_f() -> SparseTensor<f64, spenso::network::library::symbolic::ExplicitKey<AbstractIndex>> {
    initialize();
    su3_structure_f_data(CS.f_strct::<AbstractIndex>(8))
}

fn t_component(
    tensor: &SparseTensor<
        Complex<f64>,
        spenso::network::library::symbolic::ExplicitKey<AbstractIndex>,
    >,
    a: usize,
    i: usize,
    j: usize,
) -> Complex<f64> {
    tensor
        .get_owned([a, i, j])
        .unwrap_or_else(|_| Complex::new(0., 0.))
}

fn f_component(
    tensor: &SparseTensor<f64, spenso::network::library::symbolic::ExplicitKey<AbstractIndex>>,
    a: usize,
    b: usize,
    c: usize,
) -> f64 {
    tensor.get_owned([a, b, c]).unwrap_or(0.)
}

fn assert_complex_close(actual: Complex<f64>, expected: Complex<f64>) {
    assert!(
        (actual.re - expected.re).abs() < 1e-12 && (actual.im - expected.im).abs() < 1e-12,
        "expected {expected}, got {actual}"
    );
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-12,
        "expected {expected}, got {actual}"
    );
}

#[test]
fn su3_generators_have_standard_trace_normalization() {
    let t = color_t();

    for a in 0..8 {
        for b in 0..8 {
            let mut trace = Complex::new(0., 0.);
            for i in 0..3 {
                for j in 0..3 {
                    trace += t_component(&t, a, i, j) * t_component(&t, b, j, i);
                }
            }

            assert_complex_close(trace, Complex::new(if a == b { 0.5 } else { 0. }, 0.));
        }
    }
}

#[test]
fn su3_generators_satisfy_fierz_identity() {
    let t = color_t();

    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let mut contraction = Complex::new(0., 0.);
                    for a in 0..8 {
                        contraction += t_component(&t, a, i, j) * t_component(&t, a, k, l);
                    }

                    let expected = 0.5 * if i == l && k == j { 1. } else { 0. }
                        - (1. / 6.) * if i == j && k == l { 1. } else { 0. };
                    assert_complex_close(contraction, Complex::new(expected, 0.));
                }
            }
        }
    }
}

#[test]
fn su3_structure_constants_have_adjoint_casimir() {
    let f = color_f();

    for a in 0..8 {
        for b in 0..8 {
            let mut contraction = 0.;
            for c in 0..8 {
                for d in 0..8 {
                    contraction += f_component(&f, a, c, d) * f_component(&f, b, c, d);
                }
            }

            assert_close(contraction, if a == b { 3. } else { 0. });
        }
    }
}

#[test]
fn su3_generators_have_fundamental_casimir() {
    let t = color_t();

    for i in 0..3 {
        for k in 0..3 {
            let mut contraction = Complex::new(0., 0.);
            for a in 0..8 {
                for j in 0..3 {
                    contraction += t_component(&t, a, i, j) * t_component(&t, a, j, k);
                }
            }

            assert_complex_close(
                contraction,
                Complex::new(if i == k { 4. / 3. } else { 0. }, 0.),
            );
        }
    }
}
