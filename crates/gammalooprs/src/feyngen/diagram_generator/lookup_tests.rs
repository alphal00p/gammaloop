use super::*;
use crate::initialisation::test_initialise;

#[test]
fn numerator_comparisons_preserve_exact_ratios_and_failure_semantics() {
    use crate::utils::symbolica_ext::Q_I;

    test_initialise().unwrap();
    // These are already sampled scalar coefficients, not graph numerators.
    // Populate the same exact coefficient cache as the production comparison owner.
    let record = |canonical: Option<Atom>, samples: &[&str]| {
        let samples = samples
            .iter()
            .map(|sample| symbolica::parse!(*sample))
            .collect::<Vec<_>>();
        let polynomials = samples
            .iter()
            .map(|sample| sample.to_polynomial(&Q_I.clone(), None))
            .collect::<Vec<_>>();
        ProcessedNumeratorForComparison {
            diagram_id: 0,
            canonized_numerator: canonical,
            sample_evaluations_are_zero: polynomials.iter().map(|p| p.is_zero()).collect(),
            sample_evaluations: samples,
            sample_evaluations_as_polynomial: polynomials,
        }
    };

    for (name, left, right, sign, scalar) in [
        ("empty", &[][..], &[][..], None, None),
        (
            "identical",
            &["2", "3", "5"][..],
            &["2", "3", "5"][..],
            Some("1"),
            Some("1"),
        ),
        (
            "opposite",
            &["-2", "-3", "-5"][..],
            &["2", "3", "5"][..],
            Some("-1"),
            Some("-1"),
        ),
        (
            "numeric rescaling",
            &["4", "6", "10"][..],
            &["2", "3", "5"][..],
            None,
            Some("2"),
        ),
        (
            "complex rescaling",
            &["2+4𝑖", "3+6𝑖", "5+10𝑖"][..],
            &["2", "3", "5"][..],
            None,
            Some("1+2𝑖"),
        ),
        (
            "coupling rescaling",
            &["2*g^2", "3*g^2", "5*g^2"][..],
            &["2", "3", "5"][..],
            None,
            Some("g^2"),
        ),
        (
            "coupling ratio",
            &["2*g", "3*g", "5*g"][..],
            &["2*h", "3*h", "5*h"][..],
            None,
            Some("g/h"),
        ),
        (
            "coefficient polynomial fallback",
            &["g*(g+1)", "2*g*(g+1)"][..],
            &["g^2+g", "2*g^2+2*g"][..],
            Some("1"),
            Some("1"),
        ),
        (
            "negative coefficient polynomial fallback",
            &["-g*(g+1)", "-2*g*(g+1)"][..],
            &["g^2+g", "2*g^2+2*g"][..],
            Some("-1"),
            Some("-1"),
        ),
        (
            "first sample differs",
            &["7", "3", "5"][..],
            &["2", "3", "5"][..],
            None,
            None,
        ),
        (
            "second sample differs",
            &["2", "7", "5"][..],
            &["2", "3", "5"][..],
            None,
            None,
        ),
        (
            "last sample differs",
            &["2", "3", "7"][..],
            &["2", "3", "5"][..],
            None,
            None,
        ),
        (
            "inconsistent signs",
            &["2", "-3", "5"][..],
            &["2", "3", "5"][..],
            None,
            None,
        ),
        (
            "zero numerator",
            &["0", "3", "5"][..],
            &["2", "3", "5"][..],
            None,
            None,
        ),
        (
            "zero reference",
            &["2", "3", "5"][..],
            &["0", "3", "5"][..],
            None,
            None,
        ),
        (
            "all zero",
            &["0", "0"][..],
            &["0", "0"][..],
            Some("1"),
            Some("1"),
        ),
        (
            "shared zero and unit ratio",
            &["0", "3", "5"][..],
            &["0", "3", "5"][..],
            Some("1"),
            Some("1"),
        ),
        // The existing contract assigns 0/0 the ratio one, rather than skipping
        // it. A nonunit ratio at another sample must therefore reject grouping.
        (
            "shared zero and nonunit ratio",
            &["0", "6", "10"][..],
            &["0", "3", "5"][..],
            None,
            None,
        ),
        (
            "late shared zero and nonunit ratio",
            &["4", "6", "0"][..],
            &["2", "3", "0"][..],
            None,
            None,
        ),
    ] {
        let left = record(None, left);
        let right = record(None, right);
        assert_eq!(
            left.compare_with_sign_only(&right),
            sign.map(|expected| symbolica::parse!(expected)),
            "sign comparison: {name}"
        );
        assert_eq!(
            left.compare_with_scalar_rescaling(&right),
            scalar.map(|expected| symbolica::parse!(expected)),
            "scalar comparison: {name}"
        );
    }

    // A canonical symbolic match precedes sample rejection. Deliberately give
    // these records different masks to certify that the new cache prefilter
    // cannot override an already established exact match.
    let canonical = symbolica::parse_lit!(spenso::gamma(
        spenso::bis(4, 1),
        spenso::bis(4, 2),
        spenso::mink(4, 3)
    ));
    for (ratio, expected_sign) in [
        (Atom::one(), Some(Atom::one())),
        (-Atom::one(), Some(-Atom::one())),
        (symbolica::parse_lit!(g ^ 2), None),
    ] {
        let left = record(Some(&ratio * &canonical), &["0", "2"]);
        let right = record(Some(canonical.clone()), &["1", "0"]);
        assert_ne!(
            left.sample_evaluations_are_zero,
            right.sample_evaluations_are_zero
        );
        assert_eq!(left.compare_with_sign_only(&right), expected_sign);
        assert_eq!(left.compare_with_scalar_rescaling(&right), Some(ratio));
    }

    // A residual Lorentz index is not a scalar coupling ratio, so an attempted
    // symbolic comparison must still reject it when no samples are available.
    let indexed = record(Some(canonical), &[]);
    let scalar = record(Some(Atom::one()), &[]);
    assert_eq!(indexed.compare_with_sign_only(&scalar), None);
    assert_eq!(indexed.compare_with_scalar_rescaling(&scalar), None);
}
