//! Extract authenticated zero-flow physical rows with existing native Atom
//! powers. No relations among loop-dependent denominators are searched.
//! The pinned native integer-power primitive accepts at most u32::MAX.
//! This is a representation check, not a bound on rational output bytes:
//! numeric scalar materialization still consumes ordinary native CAS resources.
use feynkit_graph::IntegralFamily;
use symbolica::{
    atom::{Atom, AtomCore},
    domains::integer::Integer,
};

use crate::VakintError;

pub(super) fn extract(
    family: &IntegralFamily,
    powers: &[Integer],
    constant_slots: &[usize],
) -> Result<(Atom, Vec<Integer>), VakintError> {
    let invalid = |detail| VakintError::InvalidIntegralFormat(detail);
    let mut coefficient = Atom::one();
    let mut remaining = powers.to_vec();
    for &slot in constant_slots {
        if powers[slot] > Integer::from(0) {
            let denominator = family.denominators()[slot].expand();
            // Check before multiplication or repeated-denominator coalescing:
            // a numerator cannot legitimize an originally undefined pole.
            if denominator.is_zero() {
                return Err(invalid(
                    "a positive-power zero-flow denominator is identically zero".into(),
                ));
            }
            // Atom::pow may materialize numeric factors even in a symbolic
            // denominator, e.g. the minus sign in -m^2. Numerica::Integer::pow
            // otherwise panics here instead of returning a recoverable error.
            if powers[slot] > Integer::from(u32::MAX) {
                return Err(invalid(format!(
                    "positive zero-flow power {} exceeds the native scalar exponent representation limit {}",
                    powers[slot],
                    u32::MAX,
                )));
            }
            coefficient *= denominator.pow(Atom::num(-powers[slot].clone()));
            remaining[slot] = Integer::from(0);
        }
    }
    Ok((coefficient, remaining))
}

#[cfg(test)]
mod tests {
    use super::*;
    use feynkit_kinematics::Kinematics;
    use symbolica::symbol;

    fn family(mass_squared: i64) -> IntegralFamily {
        let q = Atom::var(symbol!("bridge_normalization::q"));
        let kinematics = Kinematics::new().with_momenta([q.clone()]).unwrap();
        let square = kinematics.scalar_product(&q, &q).unwrap();
        IntegralFamily::new(
            vec![q],
            vec![],
            vec![square - Atom::num(mass_squared), Atom::num(-mass_squared)],
            &kinematics,
        )
        .unwrap()
    }

    #[test]
    fn native_bridge_powers_are_scalar_without_changing_loop_powers() {
        for mass_squared in [1, 3] {
            let family = family(mass_squared);
            for bridge_power in [1, 2, 0, -1] {
                let powers = [Integer::from(5), Integer::from(bridge_power)];
                let (coefficient, remaining) = extract(&family, &powers, &[1]).unwrap();
                let expected = if bridge_power > 0 {
                    Atom::num(-mass_squared).pow(-bridge_power)
                } else {
                    Atom::one()
                };
                assert!((coefficient - expected).together().is_zero());
                assert_eq!(remaining[0], powers[0]);
                assert_eq!(remaining[1], Integer::from(bridge_power.min(0)));
                assert_eq!(family.loop_momenta().len(), 1);
            }
        }
    }

    #[test]
    fn native_bridge_extraction_rejects_massless_poles_and_unsupported_native_exponents() {
        let massive = family(3);
        let massless = family(0);
        assert!(extract(&massless, &[1.into(), 1.into()], &[1]).is_err());
        let unit = family(1);
        for exponent in [
            Integer::from(i64::MAX),
            Integer::from(i64::MAX) + Integer::from(1),
        ] {
            let error = extract(&unit, &[Integer::from(1), exponent.clone()], &[1]).unwrap_err();
            assert!(
                matches!(error, VakintError::InvalidIntegralFormat(message)
                if message == format!(
                    "positive zero-flow power {exponent} exceeds the native scalar exponent representation limit {}", u32::MAX)),
                "unsupported native exponents must return the specific typed ingress rejection"
            );
        }
        let (coefficient, powers) = extract(&massive, &[1.into(), 0.into()], &[1]).unwrap();
        assert_eq!(coefficient, Atom::one());
        assert_eq!(powers, [Integer::from(1), Integer::from(0)]);
    }
}
