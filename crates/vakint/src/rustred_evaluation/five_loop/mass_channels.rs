//! Combine scalar spectators with their integer common-mass scaling before IBP.

use super::SavedApplication;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::coefficient::CoefficientView;

impl SavedApplication {
    pub(super) fn mass_channel(
        spectator: &Atom,
        mass_squared: &Atom,
        exponent: i64,
        anchor: i64,
        diagnostic_max_degree: Option<usize>,
    ) -> (Atom, i64) {
        let Some(shift) = exponent.checked_sub(anchor) else {
            return (spectator.clone(), exponent);
        };
        if diagnostic_max_degree.is_some_and(|max_degree| {
            !usize::try_from(shift.unsigned_abs()).is_ok_and(|shift| shift <= max_degree)
        }) {
            return (spectator.clone(), exponent);
        }
        let mut exact = true;
        for expression in [spectator, mass_squared] {
            expression.visitor(&mut |part| {
                if let AtomView::Num(number) = part
                    && matches!(number.get_coeff_view(), CoefficientView::Float(..))
                {
                    exact = false;
                }
                exact
            });
        }
        if !exact {
            // Preserve raw floating inputs and their original precision. Staged
            // source identities remain symbolic and need no reevaluation here.
            return (spectator.clone(), exponent);
        }
        // Use one existing key's exponent as origin: raising a numeric mass to
        // -sum(a) could create an enormous rational for large dotted inputs.
        // An optional diagnostic can bound this difference; ordinary grouping
        // relies on exact integer scaling and the shared RAM/time budget.
        // Atom normalization combines powers without expanding compound masses;
        // (m²)^(Ld/2+anchor) stays outside the unchanged d-only native coefficients.
        (spectator * mass_squared.clone().pow(shift), anchor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;
    use symbolica::{domains::rational::Rational, parse};

    fn channel_sum(mass: &Atom, rows: &[(i64, Atom, Atom)]) -> (Atom, usize) {
        let mut channels = BTreeMap::<(Atom, i64), Atom>::new();
        let anchor = -rows[0].0;
        for (degree, spectator, value) in rows {
            let channel =
                SavedApplication::mass_channel(spectator, mass, -*degree, anchor, Some(512));
            *channels.entry(channel).or_insert(Atom::Zero) += value;
        }
        let n = channels.len();
        let result = channels
            .into_iter()
            .map(|((spectator, exponent), value)| {
                spectator * mass.clone().pow(parse!("10-5*ep") + exponent) * value
            })
            .sum();
        (result, n)
    }

    #[test]
    fn homogeneous_channels_preserve_the_full_five_loop_laurent_scaling() {
        let mass = parse!("mu");
        let rows = [
            (10, Atom::one(), parse!("2/ep^2+3/ep+5")),
            (11, mass.clone(), parse!("-3*(-1/ep^2+7)")),
            (12, mass.clone().pow(2), parse!("22/ep")),
        ];
        let (result, channels) = channel_sum(&mass, &rows);
        assert_eq!(channels, 1);
        let ep = symbolica::symbol!("ep");
        let expanded = result
            .series(ep, Atom::Zero.as_view(), Rational::from(0))
            .unwrap()
            .to_atom();
        let expected = parse!("5/ep^2+(25-25*log(mu))/ep-16-125*log(mu)+125/2*log(mu)^2");
        assert!((expanded - expected).expand().together().is_zero());
        // The diagnostic has engineering degree nine instead of ten, hence
        // one additional mass squared. No rank inference enters grouping.
        let diagnostic = rows.map(|(d, s, f)| (d - 1, s, f));
        let (other, channels) = channel_sum(&mass, &diagnostic);
        assert_eq!(channels, 1);
        assert!((other - &mass * result).expand().together().is_zero());
    }

    #[test]
    fn compound_and_nonunit_masses_agree_with_independent_per_integral_scaling() {
        for mass in [parse!("7/3"), parse!("x^2+y^2"), parse!("4*x^2")] {
            let rows = [
                (10, Atom::one(), parse!("f1")),
                (11, mass.clone(), parse!("-3*f2")),
                (12, mass.clone().pow(2), parse!("2*f3")),
            ];
            let before = rows
                .iter()
                .map(|(d, s, f)| s * mass.clone().pow(parse!("10-5*ep") - *d) * f)
                .sum::<Atom>();
            let (after, _) = channel_sum(&mass, &rows);
            // An exact integer dimension specialization also catches incorrect
            // scaling of nonunit numeric and compound masses without requiring
            // branch-sensitive symbolic power simplifications from the CAS.
            for epsilon in [parse!("0"), parse!("1/5"), parse!("-2/5")] {
                let residual = (before.clone() - &after)
                    .replace(parse!("ep"))
                    .with(epsilon);
                assert!(residual.expand().together().is_zero());
            }
        }
    }

    #[test]
    fn canceled_unknowns_share_a_channel_with_exact_rational_spectators() {
        let mass = parse!("mu");
        let rows = [
            (11, &mass / parse!("d-5"), parse!("unknown")),
            (10, parse!("1/(d-5)"), parse!("-unknown")),
        ];
        let (result, channels) = channel_sum(&mass, &rows);
        assert_eq!(channels, 1);
        assert!(result.is_zero());
    }

    #[test]
    fn compound_powers_stay_compact_and_staged_sources_are_not_evaluated() {
        let mass = parse!("x^2+y^2");
        let channel =
            SavedApplication::mass_channel(&mass.clone().pow(500), &mass, -510, 0, Some(512));
        assert_eq!(channel, (mass.pow(-10), 0));
        let source = parse!("source_float(17)");
        let mass = parse!("mu");
        assert_eq!(
            SavedApplication::mass_channel(
                &(&source * mass.clone().pow(2)),
                &mass,
                -12,
                -10,
                Some(2)
            ),
            (source.clone(), -10)
        );
        // No huge rational is materialized for a large absolute index; unsafe
        // differences and shifts beyond an explicit diagnostic grant stay put.
        let numeric_mass = parse!("7/3");
        for (exponent, anchor, limit) in [
            (i64::MIN, i64::MIN, 6),
            (i64::MIN, i64::MAX, 6),
            (-1_000_000_000_000, 0, 6),
            (-12, -10, 1),
        ] {
            assert_eq!(
                SavedApplication::mass_channel(
                    &source,
                    &numeric_mass,
                    exponent,
                    anchor,
                    Some(limit)
                ),
                (source.clone(), exponent)
            );
        }
    }

    #[test]
    fn ordinary_symbolic_mass_channels_ignore_dormant_degree_thresholds() {
        let source = parse!("source_float(17)");
        let mass = parse!("mu");
        let exponent = -1_000_000_000_000_i64;
        let expected = (&source * mass.clone().pow(exponent), 0);
        assert_eq!(
            SavedApplication::mass_channel(&source, &mass, exponent, 0, None),
            expected
        );
        assert_eq!(
            SavedApplication::mass_channel(&source, &mass, exponent, 0, Some(6)),
            (source.clone(), exponent)
        );
        assert_eq!(
            SavedApplication::mass_channel(&source, &mass, i64::MIN, i64::MAX, None),
            (source, i64::MIN)
        );
    }

    #[test]
    fn raw_float_components_keep_the_original_mass_factorization_and_precision() {
        use symbolica::domains::float::Float;
        let value = Atom::num(Float::parse("1.234567890123456789", Some(73)).unwrap());
        for (spectator, mass) in [
            (value.clone(), parse!("mu")),
            (parse!("x"), &value + parse!("mu")),
            (
                symbolica::symbol!("wrapped").call(value.clone()),
                parse!("mu"),
            ),
        ] {
            let original = (spectator.clone(), -12);
            for diagnostic_max_degree in [None, Some(6)] {
                assert_eq!(
                    SavedApplication::mass_channel(
                        &spectator,
                        &mass,
                        -12,
                        -10,
                        diagnostic_max_degree
                    ),
                    original
                );
            }
        }
    }
}
