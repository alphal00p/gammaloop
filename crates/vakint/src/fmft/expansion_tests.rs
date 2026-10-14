use super::*;
use crate::matad::MATAD;
use symbolica::{
    coefficient::CoefficientView,
    domains::float::{Complex, Float, Real},
};

#[test]
fn pr6_keeps_finite_and_linear_coefficients_at_their_declared_orders() {
    Vakint::initialize_vakint_symbols();
    let fmft = FMFT::with_settings(VakintSettings::default());
    let expanded = fmft
        .expand_masters(vk_parse!("PR6").unwrap().as_view())
        .unwrap();
    let series = expanded.series(vk_symbol!("ep"), 0, 2).unwrap();
    for (power, expected) in [(0, "PR6ep0"), (1, "PR6ep1"), (2, "Oep(2,I_PR6)")] {
        assert_eq!(
            series.coefficient(power.into()).unwrap(),
            vk_parse!(expected).unwrap()
        );
    }
    for name in ["PR6ep0", "PR6ep1"] {
        let literal = &MASTERS_NUMERIC_SUBSTITUTIONS[&vk_parse!(name).unwrap()].0;
        let AtomView::Num(number) = literal.as_view() else {
            panic!("expected stored coefficient")
        };
        let CoefficientView::Float(real, _) = number.get_coeff_view() else {
            panic!("expected finite source data")
        };
        let decimal_digits = f64::from(real.to_float().prec()) / std::f64::consts::LOG2_10;
        assert!((49.0..51.0).contains(&decimal_digits));
    }
}

#[test]
fn pr6_finite_and_linear_terms_match_independent_tadpole_times_five_line_results() {
    Vakint::initialize_vakint_symbols();
    let settings = VakintSettings {
        run_time_decimal_precision: 80,
        ..Default::default()
    };
    let fmft = FMFT::with_settings(settings.clone());
    let matad = MATAD::with_settings(settings);
    let expanded = fmft
        .expand_masters(vk_parse!("PR6").unwrap().as_view())
        .unwrap();
    let numeric = fmft.substitute_masters(expanded.as_view()).unwrap();
    let fmft_series = numeric.series(vk_symbol!("ep"), 0, 1).unwrap();

    // Independent lower-topology FORM reduction, in unit-mass FMFT/MATAD
    // normalization: PR6 = K1 * I3L_pinch_6(1,1,1,1,1,0).
    // Expand the exact MATAD masters before numerical HPL substitution.
    let tadpole = vk_parse!("-4*Gam(1,1)/(d^2-6*d+8)").unwrap();
    let five_line = vk_parse!(
        "-(-16*Gam(1,1)^3/(d^6-20*d^5+165*d^4-718*d^3+1736*d^2-2208*d+1152)\
        +(12-3*d)/(2*d-6)*miD5-4*Gam(1,1)/(d^2-7*d+12)*miT111\
        +(3*d-8)/(8*d-24)*miBN)"
    )
    .unwrap();
    let product = (tadpole * five_line)
        .replace(vk_parse!("d").unwrap().to_pattern())
        .with(vk_parse!("4-2*ep").unwrap().to_pattern());
    let exact_product = matad.expand_matad_masters(product.as_view()).unwrap();
    let product_series = exact_product.series(vk_symbol!("ep"), 0, 1).unwrap();
    let bits = 320;
    let tolerance = Float::parse("1e-50", Some(bits)).unwrap();
    // Separate numerical radial-K1 and seven-state five-line flows, each
    // independently refined before multiplication (70 internal digits).
    // The FMFT literals carry only ~50 digits; this test makes no higher claim.
    for (power, independent) in [
        (0, "284.47597235659825421953653349254165446876076021105"),
        (1, "836.36692108160463119910749181123692194992731952507"),
    ] {
        let exact = product_series.coefficient(power.into()).unwrap();
        let exact = matad.substitute_masters(exact.as_view()).unwrap();
        exact.visitor(&mut |part| {
            if let AtomView::Num(number) = part {
                assert!(!matches!(
                    number.get_coeff_view(),
                    CoefficientView::Float(_, _)
                ));
            }
            true
        });
        assert!(!exact.contains_symbol(vk_symbol!("Oep")));
        let hpls = matad.substitute_hpls(exact.as_view()).unwrap();
        let lower = matad
            .substitute_additional_constants(hpls.as_view())
            .unwrap();
        let lower = lower
            .evaluate_with_prec(&HashMap::<Atom, Complex<Float>>::new(), bits)
            .unwrap();
        let actual = fmft_series
            .coefficient(power.into())
            .unwrap()
            .evaluate_with_prec(&HashMap::<Atom, Complex<Float>>::new(), bits)
            .unwrap();
        let independent = Float::parse(independent, Some(bits)).unwrap();
        assert!(actual.im.as_raw().is_zero() && lower.im.as_raw().is_zero());
        assert!((actual.re.clone() - &independent).norm() < independent.norm() * &tolerance);
        assert!((actual.re - &lower.re).norm() < lower.re.norm() * &tolerance);
    }
}
