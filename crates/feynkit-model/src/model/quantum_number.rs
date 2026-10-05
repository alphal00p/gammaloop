//! Exact electric charges and optional hypercharges at the model JSON boundary.

use serde::de::Error;
use serde::{Deserialize, Deserializer, Serializer};
use symbolica::domains::integer::Integer;
use symbolica::domains::rational::Rational;

use super::quantum_number_token::RawScalar;

/// Recover an exact quantum number from a floating-point model input.
///
/// Finite integer-valued floats represent their integer value. For fractional
/// magnitudes below 4096, fractions with denominator at most one million are
/// farther apart than the binary64 spacing, so at most one rounds back to the
/// input. Other values need an explicit fraction string.
pub fn quantum_number_from_float(value: f64) -> Option<Rational> {
    if value.fract() == 0.0 {
        return Rational::try_from(value).ok();
    }
    let magnitude = value.abs();
    if !value.is_finite() || magnitude >= 4096.0 {
        return None;
    }
    let exact = Rational::try_from(magnitude).ok()?;
    let closest = exact.truncate_denominator(&Integer::from(1_000_000));
    (closest.to_f64() == magnitude).then(|| if value < 0.0 { -closest } else { closest })
}

/// Read one serialized model charge or hypercharge using the model JSON rules.
pub fn quantum_number_from_json(token: &str) -> Result<Rational, serde_json::Error> {
    let mut reader = serde_json::Deserializer::from_str(token);
    let number = deserialize(&mut reader)?;
    reader.end()?;
    Ok(number)
}

fn rational<E: Error>(raw: RawScalar) -> Result<Rational, E> {
    let token = raw.0.get();
    if token.starts_with('"') {
        let text: String = serde_json::from_str(token).map_err(E::custom)?;
        let (n, d) = text.split_once('/').unwrap_or((&text, "1"));
        let n: Integer = n.trim().parse().map_err(E::custom)?;
        let d: Integer = d.trim().parse().map_err(E::custom)?;
        if d.is_zero() {
            return Err(E::custom("quantum number denominator cannot be zero"));
        }
        return Ok(Rational::from((n, d)));
    }
    if let Ok(n) = token.parse::<i64>() {
        return Ok(n.into());
    }
    if let Ok(n) = token.parse::<u64>() {
        return Ok(n.into());
    }
    // Parse just this token with correctly rounded binary64 conversion. Other
    // numerical fields retain the surrounding deserializer's precision policy.
    let value = token.parse::<f64>().map_err(E::custom)?;
    if value.abs() <= 9_007_199_254_740_991.0
        && let Some(exact) = quantum_number_from_float(value)
    {
        return Ok(exact);
    }
    Err(E::custom(format!(
        "unsupported floating-point quantum number {token}. Supply an exact fraction string such as '2/3'"
    )))
}

pub(super) fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<Rational, D::Error> {
    rational(RawScalar::deserialize(d)?)
}

pub(super) fn serialize<S: Serializer>(value: &Rational, s: S) -> Result<S::Ok, S::Error> {
    s.collect_str(value)
}

pub(super) mod optional {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    use symbolica::domains::rational::Rational;

    use super::{RawScalar, rational};

    pub(crate) fn deserialize<'de, D: Deserializer<'de>>(
        d: D,
    ) -> Result<Option<Rational>, D::Error> {
        Option::<RawScalar>::deserialize(d)?
            .map(rational)
            .transpose()
    }

    pub(crate) fn serialize<S: Serializer>(
        value: &Option<Rational>,
        s: S,
    ) -> Result<S::Ok, S::Error> {
        value.as_ref().map(ToString::to_string).serialize(s)
    }
}

#[cfg(test)]
mod tests {
    use serde::Deserialize;
    use serde_json::{from_reader, from_slice, from_str, from_value, json};
    use symbolica::domains::rational::Rational;

    #[derive(Deserialize)]
    struct Record {
        #[serde(deserialize_with = "super::deserialize")]
        charge: Rational,
        #[serde(deserialize_with = "super::optional::deserialize")]
        hypercharge: Option<Rational>,
        ordinary: f64,
    }

    #[test]
    fn decimal_recovery_is_local_to_quantum_number_fields() {
        let token = "0.09090909090909091";
        let text = format!(r#"{{"charge":{token},"hypercharge":null,"ordinary":{token}}}"#);
        let ambient = from_str::<f64>(token).unwrap();
        let expected = Rational::from((1, 11));
        for record in [
            from_str::<Record>(&text).unwrap(),
            from_slice(text.as_bytes()).unwrap(),
            from_reader(text.as_bytes()).unwrap(),
        ] {
            assert_eq!(record.charge, expected);
            assert!(record.hypercharge.is_none());
            assert_eq!(record.ordinary.to_bits(), ambient.to_bits());
        }

        let record: Record = from_value(json!({
            "charge": 1.0 / 11.0,
            "hypercharge": -2.0 / 3.0,
            "ordinary": ambient,
        }))
        .unwrap();
        assert_eq!(record.charge, expected);
        assert_eq!(record.hypercharge, Some(Rational::from((-2, 3))));
        assert_eq!(record.ordinary.to_bits(), ambient.to_bits());
    }
}
